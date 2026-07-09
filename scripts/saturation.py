#!/usr/bin/env python3
"""Substitution-saturation test for a marker (MarkerVet stage 3).

Saturation -- multiple substitutions at the same site erasing historical signal
-- is the second classic failure mode for deep phylogeny (alongside
compositional attraction).  A saturated marker records recent noise, not ancient
divergence, and misleads both concatenation and coalescence.

We measure saturation the way Philippe and colleagues do: compare *observed*
distances (uncorrected p-distances read straight off the alignment) with
*inferred* distances (patristic distances -- summed branch lengths -- from the
model-corrected gene tree).  For each pair of taxa:

    x = patristic distance (tree)      # what the model thinks really happened
    y = p-distance        (alignment)  # what we can still see

Fit y = slope * x  through the cloud of pairwise points.  As divergence grows,
observed p-distances plateau (they cannot exceed ~1) while patristic distances
keep climbing, so the slope falls toward 0.  A steep slope (near 1) means little
saturation; a shallow slope means heavy saturation.

    fail (too saturated) if slope < min_slope.

We report the regression slope, its R^2, and the number of pairs used.  Using
the slope (rather than a raw correlation) makes the statistic comparable across
markers with different overall divergence.

Output: one-row TSV (+ optional JSON):

    marker  n_pairs  slope  r_squared  mean_pdist  mean_patristic  pass
"""
from __future__ import print_function

import argparse
import itertools
import json
import os
import sys
from collections import OrderedDict

import numpy as np

try:
    import dendropy
except ImportError:  # pragma: no cover - dendropy is an env dependency
    dendropy = None


def read_fasta(path):
    """Minimal FASTA reader -> OrderedDict{name: sequence}."""
    seqs = OrderedDict()
    name, chunks = None, []
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


GAP = set("-.?")


def p_distance(a, b):
    """Uncorrected p-distance between two aligned sequences.

    Only columns where both taxa have a non-gap, non-ambiguous residue are
    counted (pairwise deletion).  Returns None if there is no overlap.
    """
    diff = 0
    comp = 0
    for ca, cb in zip(a, b):
        if ca in GAP or cb in GAP or ca == "X" or cb == "X":
            continue
        comp += 1
        if ca != cb:
            diff += 1
    if comp == 0:
        return None
    return diff / float(comp)


def patristic_distances(tree_path, taxa):
    """Return {frozenset({t1, t2}): patristic_distance} for taxa pairs.

    Uses dendropy's phylogenetic distance matrix over the shared taxon set.
    """
    tree = dendropy.Tree.get(path=tree_path, schema="newick",
                             preserve_underscores=True)
    pdm = tree.phylogenetic_distance_matrix()
    label_to_taxon = {t.label: t for t in tree.taxon_namespace}
    out = {}
    present = [t for t in taxa if t in label_to_taxon]
    for t1, t2 in itertools.combinations(present, 2):
        d = pdm.distance(label_to_taxon[t1], label_to_taxon[t2])
        out[frozenset((t1, t2))] = d
    return out


def slope_through_origin(x, y):
    """Least-squares slope of y = slope * x forced through the origin, plus R^2.

    A through-origin model is appropriate here because zero evolutionary
    distance must give zero observed difference.  R^2 is computed relative to
    that same origin-anchored model.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    denom = float((x * x).sum())
    if denom == 0:
        return float("nan"), float("nan")
    slope = float((x * y).sum() / denom)
    resid = y - slope * x
    ss_res = float((resid ** 2).sum())
    ss_tot = float((y ** 2).sum())
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return slope, r2


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alignment", required=True, help="aligned protein FASTA")
    ap.add_argument("--tree", required=True,
                    help="gene tree with branch lengths (Newick)")
    ap.add_argument("--marker", required=True)
    ap.add_argument("--min-slope", type=float, default=0.55,
                    help="FAIL if regression slope < this (default 0.55)")
    ap.add_argument("--max-pairs", type=int, default=20000,
                    help="cap on taxon pairs sampled for speed (default 20000)")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-json", default=None)
    args = ap.parse_args(argv)

    if dendropy is None:
        sys.exit("saturation.py requires dendropy (see envs/markervet.yaml)")

    seqs = read_fasta(args.alignment)
    taxa = list(seqs.keys())
    patr = patristic_distances(args.tree, taxa)

    xs, ys = [], []
    pairs = list(patr.keys())
    # Deterministically subsample pairs if the alignment is huge (evenly spaced
    # rather than random so results are reproducible across runs).
    if len(pairs) > args.max_pairs:
        step = len(pairs) / float(args.max_pairs)
        pairs = [pairs[int(i * step)] for i in range(args.max_pairs)]

    for pair in pairs:
        t1, t2 = tuple(pair)
        pd = p_distance(seqs[t1], seqs[t2])
        if pd is None:
            continue
        xs.append(patr[pair])
        ys.append(pd)

    if len(xs) < 3:
        slope, r2 = float("nan"), float("nan")
        passed = False
    else:
        slope, r2 = slope_through_origin(xs, ys)
        passed = bool(slope >= args.min_slope)

    row = OrderedDict([
        ("marker", args.marker),
        ("n_pairs", len(xs)),
        ("slope", round(slope, 4) if slope == slope else "NA"),
        ("r_squared", round(r2, 4) if r2 == r2 else "NA"),
        ("mean_pdist", round(float(np.mean(ys)), 4) if ys else "NA"),
        ("mean_patristic", round(float(np.mean(xs)), 4) if xs else "NA"),
        ("min_slope", args.min_slope),
        ("pass", passed),
    ])

    os.makedirs(os.path.dirname(os.path.abspath(args.out_tsv)), exist_ok=True)
    with open(args.out_tsv, "w") as fh:
        fh.write("\t".join(row.keys()) + "\n")
        fh.write("\t".join(str(v) for v in row.values()) + "\n")
    if args.out_json:
        with open(args.out_json, "w") as fh:
            json.dump(row, fh, indent=2)

    print("[saturation] {m}: slope={s} (min {mn}) -> {v}".format(
        m=args.marker, s=row["slope"], mn=args.min_slope,
        v="PASS" if passed else "FAIL"), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
