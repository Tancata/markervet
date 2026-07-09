#!/usr/bin/env python3
"""Aggregate all screens into a keep/drop scorecard (MarkerVet stage 8).

This is where MarkerVet's **aggressive** drop policy lives.  For concatenation,
even a little incongruent data biases the supermatrix tree, so the policy is
conservative-keep / aggressive-drop:

    DROP a marker if it fails ANY enabled criterion.

The criteria, each read from the per-marker outputs of earlier stages and each
switched on/off (with its threshold) in the regime YAML:

    compositional   -- chi-square homogeneity fail            (artefact screen)
    saturation      -- regression slope below threshold       (artefact screen)
    congruence      -- high-support RF / quartet conflict vs reference
    monophyly       -- clade_recovery_score below threshold
    phylter         -- gene-level outlier, or > X% taxon cells removed
    reconciliation  -- inferred transfers above threshold

For robustness the policy is *fail-closed*: if an enabled criterion's result
file is missing or unreadable, the marker is dropped with a `missing:<crit>`
reason rather than silently kept -- a marker we could not vet is not a marker we
trust for deep phylogeny.

Output: scorecard.tsv, one row per marker, with the decision, the number of
failed criteria, a human-readable reason list, and the key statistic behind
each criterion.
"""
from __future__ import print_function

import argparse
import csv
import os
import sys
from collections import OrderedDict

import yaml


def read_one_row_tsv(path):
    """Read a one-data-row TSV into a dict; return None if absent/empty."""
    if not path or not os.path.exists(path):
        return None
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return rows[0] if rows else None


def as_bool(v):
    """Interpret assorted truthy/falsey string spellings as bool, else None."""
    if v is None:
        return None
    s = str(v).strip().lower()
    if s in ("true", "1", "yes", "pass"):
        return True
    if s in ("false", "0", "no", "fail"):
        return False
    return None


def as_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def load_phylter_summary(path):
    """Return {marker: {'gene_outlier': bool, 'frac_removed': float}}."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["marker"]] = {
                "gene_outlier": as_bool(row.get("gene_outlier")),
                "frac_removed": as_float(row.get("frac_removed")),
                "cells_removed": row.get("cells_removed"),
            }
    return out


def evaluate_marker(marker, cfg, paths, phylter):
    """Return (decision, reasons, metrics) for one marker under the policy."""
    reasons = []
    metrics = OrderedDict()

    # ---- 1. compositional homogeneity --------------------------------------
    if cfg["artefact"]["compositional"]["enabled"]:
        row = read_one_row_tsv(paths["comp"])
        metrics["comp_p"] = row.get("p_value") if row else "NA"
        passed = as_bool(row.get("pass")) if row else None
        if passed is None:
            reasons.append("missing:compositional")
        elif not passed:
            reasons.append("compositional_heterogeneous(p=%s)" % metrics["comp_p"])

    # ---- 2. saturation ------------------------------------------------------
    if cfg["artefact"]["saturation"]["enabled"]:
        row = read_one_row_tsv(paths["sat"])
        metrics["sat_slope"] = row.get("slope") if row else "NA"
        passed = as_bool(row.get("pass")) if row else None
        if passed is None:
            reasons.append("missing:saturation")
        elif not passed:
            reasons.append("saturated(slope=%s)" % metrics["sat_slope"])

    # ---- 3. congruence (RF + quartet vs reference) -------------------------
    if cfg["congruence"]["enabled"]:
        row = read_one_row_tsv(paths["cong"])
        metrics["rf_norm"] = row.get("rf_norm") if row else "NA"
        metrics["quartet_norm"] = row.get("quartet_norm") if row else "NA"
        passed = as_bool(row.get("pass")) if row else None
        if passed is None:
            reasons.append("missing:congruence")
        elif not passed:
            rf_ok = as_bool(row.get("rf_pass"))
            q_ok = as_bool(row.get("quartet_pass"))
            if rf_ok is False:
                reasons.append("rf_conflict(%s)" % metrics["rf_norm"])
            if q_ok is False:
                reasons.append("quartet_conflict(%s)" % metrics["quartet_norm"])
            if rf_ok is None and q_ok is None:
                reasons.append("congruence_fail")

    # ---- 4. monophyly / clade recovery -------------------------------------
    if cfg["monophyly"].get("clades"):
        row = read_one_row_tsv(paths["mono"])
        metrics["clade_recovery_score"] = row.get("clade_recovery_score") if row else "NA"
        metrics["intruders"] = row.get("total_intruders") if row else "NA"
        metrics["outliers"] = row.get("total_outliers") if row else "NA"
        passed = as_bool(row.get("pass")) if row else None
        if passed is None:
            reasons.append("missing:monophyly")
        elif not passed:
            reasons.append("low_clade_recovery(%s)" % metrics["clade_recovery_score"])

    # ---- 5. PhylteR outlier vote -------------------------------------------
    if cfg["outliers_hgt"]["phylter"]["enabled"]:
        pr = phylter.get(marker)
        frac_max = cfg["outliers_hgt"]["phylter"]["taxon_cell_frac_max"]
        if pr is None:
            reasons.append("missing:phylter")
        else:
            metrics["phylter_frac_removed"] = pr["frac_removed"]
            if pr["gene_outlier"]:
                reasons.append("phylter_gene_outlier")
            elif pr["frac_removed"] is not None and pr["frac_removed"] > frac_max:
                reasons.append("phylter_cells_removed(%.2f>%.2f)"
                               % (pr["frac_removed"], frac_max))

    # ---- 6. reconciliation transfer vote -----------------------------------
    if cfg["outliers_hgt"]["reconciliation"]["enabled"]:
        row = read_one_row_tsv(paths["recon"])
        t_max = cfg["outliers_hgt"]["reconciliation"]["transfers_max"]
        transfers = as_float(row.get("transfers")) if row else None
        metrics["transfers"] = transfers if transfers is not None else "NA"
        if transfers is None:
            reasons.append("missing:reconciliation")
        elif transfers > t_max:
            reasons.append("excess_transfers(%g>%g)" % (transfers, t_max))

    # ---- decision -----------------------------------------------------------
    # Aggressive policy: any reason at all => DROP.
    decision = "DROP" if reasons else "KEEP"
    return decision, reasons, metrics


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True, help="regime YAML")
    ap.add_argument("--markers", required=True,
                    help="comma-separated marker names, or a file (one per line)")
    ap.add_argument("--comp-dir", required=True)
    ap.add_argument("--sat-dir", required=True)
    ap.add_argument("--cong-dir", required=True)
    ap.add_argument("--mono-dir", required=True)
    ap.add_argument("--recon-dir", required=True)
    ap.add_argument("--phylter-summary", required=True)
    ap.add_argument("--out", required=True, help="scorecard.tsv")
    args = ap.parse_args(argv)

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)

    if os.path.exists(args.markers):
        with open(args.markers) as fh:
            markers = [ln.strip() for ln in fh if ln.strip()]
    else:
        markers = [m for m in args.markers.split(",") if m]

    phylter = load_phylter_summary(args.phylter_summary)

    # Assemble scorecard rows.
    all_metric_keys = []
    scored = []
    for marker in markers:
        paths = {
            "comp": os.path.join(args.comp_dir, marker + ".tsv"),
            "sat": os.path.join(args.sat_dir, marker + ".tsv"),
            "cong": os.path.join(args.cong_dir, marker + ".tsv"),
            "mono": os.path.join(args.mono_dir, marker + ".summary.tsv"),
            "recon": os.path.join(args.recon_dir, marker + ".tsv"),
        }
        decision, reasons, metrics = evaluate_marker(marker, cfg, paths, phylter)
        for k in metrics:
            if k not in all_metric_keys:
                all_metric_keys.append(k)
        scored.append((marker, decision, reasons, metrics))

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    header = ["marker", "decision", "n_failed", "reasons"] + all_metric_keys
    n_keep = 0
    with open(args.out, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for marker, decision, reasons, metrics in scored:
            if decision == "KEEP":
                n_keep += 1
            row = [marker, decision, str(len(reasons)),
                   ";".join(reasons) if reasons else "-"]
            row += [str(metrics.get(k, "NA")) for k in all_metric_keys]
            fh.write("\t".join(row) + "\n")

    print("[aggregate_scorecard] {k}/{n} markers KEEP (aggressive policy)".format(
        k=n_keep, n=len(markers)), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
