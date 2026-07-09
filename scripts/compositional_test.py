#!/usr/bin/env python3
"""Compositional homogeneity test for a protein alignment (MarkerVet stage 3).

Deep phylogenetic signal is easily corrupted by *compositional attraction*:
lineages that independently drift to similar amino-acid compositions (e.g.
thermophiles, or fast-evolving eukaryotes) are pulled together regardless of
true relationships.  A marker whose taxa differ significantly in amino-acid
composition is therefore a liability for concatenation.

This script implements a chi-square test of homogeneity across taxa, the same
idea used by PhyloBayes / IQ-TREE's composition test: build the taxon x
amino-acid count matrix and ask whether all rows (taxa) are drawn from the same
underlying composition.

    H0 : every taxon shares the alignment-wide amino-acid composition.
    fail (reject H0, composition is heterogeneous) if p < alpha.

The statistic is

    X^2 = sum_{taxon t, residue a}  (O_ta - E_ta)^2 / E_ta

with expected counts E_ta = (row_total_t * col_total_a) / grand_total, i.e. the
standard contingency-table expectation under independence of taxon and residue.
Degrees of freedom = (n_taxa - 1) * (n_states - 1), counting only the amino-acid
columns that actually occur.

Output is a one-row TSV (plus a JSON sidecar) so the aggregator can consume it:

    marker  n_taxa  n_sites  chi2  df  p_value  min_expected  pass

`pass` is True when the composition is acceptably homogeneous (p >= alpha).
"""
from __future__ import print_function

import argparse
import json
import os
import sys
from collections import OrderedDict

import numpy as np
from scipy import stats

# The 20 canonical amino acids.  Ambiguity codes (BZJX*), gaps (-.) and stops
# are ignored: they carry no clean compositional signal and would inflate df.
AA = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {a: i for i, a in enumerate(AA)}


def read_fasta(path):
    """Minimal FASTA reader -> OrderedDict{name: sequence}.

    Kept dependency-free (no Biopython) so the artefact screen can run in a
    bare Python environment.
    """
    seqs = OrderedDict()
    name = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def composition_matrix(seqs):
    """Return (taxa, counts) where counts[t, a] is the residue count.

    Rows are taxa, columns are the 20 amino acids.  Non-standard characters are
    silently skipped.
    """
    taxa = list(seqs.keys())
    counts = np.zeros((len(taxa), len(AA)), dtype=float)
    for i, t in enumerate(taxa):
        for ch in seqs[t].upper():
            j = AA_INDEX.get(ch)
            if j is not None:
                counts[i, j] += 1.0
    return taxa, counts


def chi_square_homogeneity(counts):
    """Chi-square test of taxon x residue homogeneity.

    Drops all-zero taxa (no observed residues) and all-zero residue columns
    (absent amino acids) before computing expectations, so the degrees of
    freedom reflect only informative cells.  Returns
    (chi2, df, p_value, min_expected, n_taxa_used, n_states_used).
    """
    # Remove taxa / columns that are entirely empty -- they contribute zero
    # expected counts and would create 0/0 cells.
    row_ok = counts.sum(axis=1) > 0
    counts = counts[row_ok]
    col_ok = counts.sum(axis=0) > 0
    counts = counts[:, col_ok]

    n_taxa, n_states = counts.shape
    if n_taxa < 2 or n_states < 2:
        # Not enough variation to test -- treat as "no evidence of heterogeneity".
        return 0.0, 0, 1.0, float("nan"), n_taxa, n_states

    grand = counts.sum()
    row_tot = counts.sum(axis=1, keepdims=True)
    col_tot = counts.sum(axis=0, keepdims=True)
    expected = row_tot * col_tot / grand

    chi2 = float(((counts - expected) ** 2 / expected).sum())
    df = (n_taxa - 1) * (n_states - 1)
    p_value = float(stats.chi2.sf(chi2, df))
    return chi2, df, p_value, float(expected.min()), n_taxa, n_states


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alignment", required=True, help="aligned protein FASTA")
    ap.add_argument("--marker", required=True, help="marker name (for the report)")
    ap.add_argument("--alpha", type=float, default=0.05,
                    help="p-value threshold; FAIL if p < alpha (default 0.05)")
    ap.add_argument("--out-tsv", required=True, help="one-row TSV result")
    ap.add_argument("--out-json", default=None, help="optional JSON sidecar")
    args = ap.parse_args(argv)

    seqs = read_fasta(args.alignment)
    _, counts = composition_matrix(seqs)
    chi2, df, p_value, min_exp, n_taxa, n_states = chi_square_homogeneity(counts)

    # p >= alpha => homogeneous enough => pass.
    passed = bool(p_value >= args.alpha)

    n_sites = max((len(s) for s in seqs.values()), default=0)
    row = OrderedDict([
        ("marker", args.marker),
        ("n_taxa", n_taxa),
        ("n_sites", n_sites),
        ("chi2", round(chi2, 4)),
        ("df", df),
        ("p_value", p_value),
        ("min_expected", min_exp),
        ("alpha", args.alpha),
        ("pass", passed),
    ])

    os.makedirs(os.path.dirname(os.path.abspath(args.out_tsv)), exist_ok=True)
    with open(args.out_tsv, "w") as fh:
        fh.write("\t".join(row.keys()) + "\n")
        fh.write("\t".join(str(v) for v in row.values()) + "\n")

    if args.out_json:
        with open(args.out_json, "w") as fh:
            json.dump(row, fh, indent=2)

    print("[compositional_test] {marker}: chi2={chi2:.1f} df={df} "
          "p={p:.3g} -> {verdict}".format(marker=args.marker, chi2=chi2, df=df,
                                          p=p_value,
                                          verdict="PASS" if passed else "FAIL"),
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
