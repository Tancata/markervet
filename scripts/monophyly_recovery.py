#!/usr/bin/env python3
"""MonoPhy-style clade-recovery scoring for a gene tree (MarkerVet stage 6).

Whether a marker recovers the clades we are confident about (d__Archaea,
d__Bacteria, Asgard, DPANN, Eukaryota, ...) is simultaneously a *filter* (a
marker that scatters known clades is untrustworthy for deep phylogeny) and a
*scientific result* (the clade x marker recovery matrix is a headline output).

For each named clade we reproduce the logic of the MonoPhy R package
(Schwery & O'Meara 2016):

  * members present  -- leaves assigned to the clade;
  * intruders        -- foreign leaves nested inside the clade's MRCA span;
  * outliers         -- clade members that fall outside the clade's largest
                        monophyletic block (i.e. placed among other groups);
  * recovery fraction = (size of the largest all-member block) / members present.

The comparison is **support-gated**: before scoring, internal branches whose
support (UFBoot, parsed from the IQ-TREE "aLRT/UFBoot" label) falls below
`--support-gate` are collapsed to polytomies, so we never penalise a marker for
non-monophyly that its own data do not strongly support.  (Stage 4 already
collapses by the tree-building gate; this re-applies the -- possibly different --
monophyly gate.)

Two outputs:
  * a long clade x marker table (one row per clade) -> assembled into the
    clade x marker matrix by the aggregate/report stage;
  * a one-row per-marker summary carrying `clade_recovery_score`
    (mean recovery over clades present) and the keep/drop `pass` flag.
"""
from __future__ import print_function

import argparse
import json
import os
import re
import sys
from collections import OrderedDict

from ete3 import Tree

_SUPPORT_RE = re.compile(r"^\s*([0-9]*\.?[0-9]+)\s*(?:/\s*([0-9]*\.?[0-9]+)\s*)?$")


def parse_ufboot(label):
    """Extract UFBoot (second field of 'aLRT/UFBoot', else the lone value)."""
    if not label:
        return None
    m = _SUPPORT_RE.match(str(label))
    if not m:
        return None
    return float(m.group(2)) if m.group(2) is not None else float(m.group(1))


def collapse_below(tree, gate):
    """Collapse internal branches whose parsed UFBoot support < gate.

    Mirrors scripts/collapse_support.py but applies the monophyly-specific gate,
    so weakly held groupings become polytomies and cannot create spurious
    intruders/outliers.
    """
    to_collapse = []
    for node in tree.traverse("postorder"):
        if node.is_leaf() or node.is_root():
            continue
        sup = parse_ufboot(node.name)
        if sup is not None and sup < gate:
            to_collapse.append(node)
    for node in to_collapse:
        node.delete(prevent_nondicotomic=False, preserve_branch_length=True)


def load_membership(taxon_table, clades):
    """Return a function leaf_name -> clade_name (or None).

    Each clade entry is {name, pattern, column?}.  `pattern` is a regex matched
    against the leaf name, or against `column` of the taxon_table when given.
    The first matching clade in list order wins (so list more specific clades
    first if they overlap).
    """
    # Optional per-taxon annotation columns.
    table = {}
    header = []
    if taxon_table and os.path.exists(taxon_table):
        with open(taxon_table) as fh:
            for i, line in enumerate(fh):
                parts = line.rstrip("\n").split("\t")
                if i == 0:
                    header = parts
                    continue
                if parts and parts[0]:
                    table[parts[0]] = dict(zip(header[1:], parts[1:]))

    compiled = []
    for c in clades:
        compiled.append((c["name"], re.compile(c["pattern"]), c.get("column")))

    def assign(leaf):
        for name, rx, column in compiled:
            target = leaf
            if column and leaf in table:
                target = table[leaf].get(column, "")
            if rx.search(target):
                return name
        return None

    return assign


def largest_member_block(tree, clade):
    """Size of the biggest maximal subtree whose leaves are all `clade` members.

    ete3.get_monophyletic returns every maximal node all of whose leaves carry
    the target attribute value; the largest of these is the clade's recovered
    core.  Returns 0 if the clade is absent.
    """
    best = 0
    for node in tree.get_monophyletic(values=[clade], target_attr="clade"):
        n = len(node.get_leaves())
        if n > best:
            best = n
    return best


def score_clade(tree, clade, members_present):
    """Compute recovery/intruder/outlier statistics for one clade."""
    m = len(members_present)
    if m == 0:
        return OrderedDict([("present", 0), ("members", 0), ("recovered", 0),
                            ("recovery", "NA"), ("intruders", "NA"),
                            ("outliers", "NA"), ("status", "absent")])
    if m == 1:
        # A singleton is trivially "monophyletic".
        return OrderedDict([("present", 1), ("members", 1), ("recovered", 1),
                            ("recovery", 1.0), ("intruders", 0),
                            ("outliers", 0), ("status", "singleton")])

    is_mono, _, _ = tree.check_monophyly(values=[clade], target_attr="clade",
                                         unrooted=True)
    if is_mono:
        return OrderedDict([("present", m), ("members", m), ("recovered", m),
                            ("recovery", 1.0), ("intruders", 0),
                            ("outliers", 0), ("status", "monophyletic")])

    # Non-monophyletic: intruders live inside the members' MRCA; outliers are
    # members outside the largest all-member block.
    mrca = tree.get_common_ancestor(list(members_present))
    span_leaves = mrca.get_leaf_names()
    intruders = sum(1 for l in span_leaves if l not in members_present)
    recovered = largest_member_block(tree, clade)
    outliers = m - recovered
    recovery = recovered / float(m)
    status = "intruded" if intruders and outliers == 0 else "non-monophyletic"
    return OrderedDict([("present", m), ("members", m), ("recovered", recovered),
                        ("recovery", round(recovery, 4)),
                        ("intruders", intruders), ("outliers", outliers),
                        ("status", status)])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree", required=True, help="collapsed gene tree (Newick)")
    ap.add_argument("--marker", required=True)
    ap.add_argument("--clades-json", required=True,
                    help="JSON list of {name, pattern, column?} clade definitions")
    ap.add_argument("--taxon-table", default=None,
                    help="optional TSV of per-taxon annotations")
    ap.add_argument("--support-gate", type=float, default=80.0,
                    help="collapse branches with UFBoot < this before scoring")
    ap.add_argument("--recovery-min", type=float, default=0.80,
                    help="drop marker if clade_recovery_score < this")
    ap.add_argument("--out-matrix", required=True,
                    help="long clade x marker TSV (one row per clade)")
    ap.add_argument("--out-summary", required=True,
                    help="one-row per-marker summary TSV")
    ap.add_argument("--out-json", default=None)
    args = ap.parse_args(argv)

    with open(args.clades_json) as fh:
        clades = json.load(fh)

    tree = Tree(args.tree, format=1)
    collapse_below(tree, args.support_gate)

    assign = load_membership(args.taxon_table, clades)
    # Label every leaf with its clade (or "none") for ete3's monophyly checks.
    present_by_clade = {c["name"]: set() for c in clades}
    for leaf in tree:
        cl = assign(leaf.name)
        leaf.add_feature("clade", cl if cl is not None else "none")
        if cl in present_by_clade:
            present_by_clade[cl].add(leaf.name)

    # Score every clade and accumulate the per-marker headline.
    matrix_rows = []
    recoveries = []
    total_intruders = 0
    total_outliers = 0
    for c in clades:
        name = c["name"]
        stats = score_clade(tree, name, present_by_clade[name])
        row = OrderedDict([("marker", args.marker), ("clade", name)])
        row.update(stats)
        matrix_rows.append(row)
        if stats["recovery"] != "NA":
            recoveries.append(float(stats["recovery"]))
            total_intruders += int(stats["intruders"])
            total_outliers += int(stats["outliers"])

    # clade_recovery_score: mean recovery over the clades actually present.
    score = sum(recoveries) / len(recoveries) if recoveries else float("nan")
    passed = bool(recoveries and score >= args.recovery_min)

    # --- write the long clade x marker table --------------------------------
    os.makedirs(os.path.dirname(os.path.abspath(args.out_matrix)), exist_ok=True)
    cols = ["marker", "clade", "present", "members", "recovered", "recovery",
            "intruders", "outliers", "status"]
    with open(args.out_matrix, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in matrix_rows:
            fh.write("\t".join(str(row[c]) for c in cols) + "\n")

    # --- write the one-row per-marker summary -------------------------------
    summary = OrderedDict([
        ("marker", args.marker),
        ("n_clades_present", len(recoveries)),
        ("clade_recovery_score", round(score, 4) if score == score else "NA"),
        ("total_intruders", total_intruders),
        ("total_outliers", total_outliers),
        ("recovery_min", args.recovery_min),
        ("pass", passed),
    ])
    with open(args.out_summary, "w") as fh:
        fh.write("\t".join(summary.keys()) + "\n")
        fh.write("\t".join(str(v) for v in summary.values()) + "\n")
    if args.out_json:
        with open(args.out_json, "w") as fh:
            json.dump({"summary": summary, "clades": matrix_rows}, fh, indent=2)

    print("[monophyly_recovery] {m}: score={s} intruders={i} outliers={o} -> {v}"
          .format(m=args.marker, s=summary["clade_recovery_score"],
                  i=total_intruders, o=total_outliers,
                  v="PASS" if passed else "FAIL"), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
