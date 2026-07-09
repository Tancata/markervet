#!/usr/bin/env python3
"""Topological distance of a gene tree from the reference (MarkerVet stage 5).

After weak branches have been collapsed (stage 4), a *remaining* disagreement
with the trusted reference topology is real conflict, not noise -- exactly the
kind of signal that flags a horizontally transferred or otherwise misleading
marker.  This script quantifies that disagreement two complementary ways, both
via dendropy on the shared taxon set:

  * Normalised Robinson-Foulds (RF) -- symmetric difference of bipartitions,
    divided by the maximum possible (2 * (n_taxa - 3) for two binary trees).
    Sensitive to any misplaced clade.

  * Normalised quartet distance -- fraction of the C(n,4) taxon quartets whose
    resolved topology differs between the two trees.  More robust than RF to a
    single rogue tip and interpretable across trees of different sizes.

Both are computed on the induced topologies restricted to taxa shared by the
gene tree and the reference (deep marker sets are patchy, so this pruning is
essential).  Collapsed polytomies are handled by dendropy's unrooted RF; the
quartet distance treats unresolved quartets as "not conflicting".

Output: one-row TSV (+ optional JSON):

    marker  n_shared  rf  rf_norm  quartet  quartet_norm  rf_pass  quartet_pass  pass

`pass` is the AND of the enabled sub-tests (a marker passes congruence only if
it is within BOTH the RF and quartet ceilings).
"""
from __future__ import print_function

import argparse
import itertools
import json
import os
import sys
from collections import OrderedDict

import dendropy
from dendropy.calculate import treecompare


def load_common(gene_path, ref_path):
    """Load gene + reference trees onto a shared taxon namespace, pruned to the
    intersection of their leaf sets.  Returns (gene, ref, shared_labels)."""
    tns = dendropy.TaxonNamespace()
    gene = dendropy.Tree.get(path=gene_path, schema="newick",
                             taxon_namespace=tns, preserve_underscores=True)
    ref = dendropy.Tree.get(path=ref_path, schema="newick",
                            taxon_namespace=tns, preserve_underscores=True)

    gene_labels = {l.taxon.label for l in gene.leaf_node_iter() if l.taxon}
    ref_labels = {l.taxon.label for l in ref.leaf_node_iter() if l.taxon}
    shared = gene_labels & ref_labels

    keep = [t for t in tns if t.label in shared]
    gene.retain_taxa(keep)
    ref.retain_taxa(keep)
    # Recompute bipartitions after pruning so RF is defined on the induced trees.
    gene.encode_bipartitions()
    ref.encode_bipartitions()
    return gene, ref, sorted(shared)


def normalised_rf(gene, ref, n_shared):
    """Symmetric-difference RF and its normalisation by 2*(n-3)."""
    rf = treecompare.symmetric_difference(gene, ref)
    max_rf = 2 * (n_shared - 3)
    rf_norm = rf / float(max_rf) if max_rf > 0 else 0.0
    return rf, rf_norm


def _quartet_resolution(tree, labels):
    """Map each taxon quartet to its resolved partner-pairing, or None if the
    quartet is unresolved on this tree.

    For four leaves {a,b,c,d} an unrooted binary tree induces exactly one of
    three pairings: ab|cd, ac|bd, ad|bc.  We identify the pairing by comparing
    pairwise path lengths (in edges): the two leaves on the same side of the
    central branch are closer to each other than to the other pair.  A polytomy
    yields ties, which we report as unresolved (None).
    """
    pdm = tree.phylogenetic_distance_matrix()
    taxon = {t.label: t for t in tree.taxon_namespace if t.label in set(labels)}
    res = {}
    for quad in itertools.combinations(labels, 4):
        a, b, c, d = (taxon[x] for x in quad)
        # steps=True -> topological (edge-count) distance, ignoring branch length.
        d_ab = pdm.path_edge_count(a, b) + pdm.path_edge_count(c, d)
        d_ac = pdm.path_edge_count(a, c) + pdm.path_edge_count(b, d)
        d_ad = pdm.path_edge_count(a, d) + pdm.path_edge_count(b, c)
        m = min(d_ab, d_ac, d_ad)
        # Ambiguous (tie) => unresolved on this tree.
        if [d_ab, d_ac, d_ad].count(m) > 1:
            res[quad] = None
        elif m == d_ab:
            res[quad] = "ab|cd"
        elif m == d_ac:
            res[quad] = "ac|bd"
        else:
            res[quad] = "ad|bc"
    return res


def quartet_distance(gene, ref, labels, max_quartets=200000):
    """Fraction of resolved quartets that conflict between the two trees.

    Quartets unresolved on *either* tree (a collapsed polytomy) are not counted
    as conflict -- consistent with the "collapse then compare" principle.  For
    large taxon sets the quartet space is deterministically subsampled.
    """
    quads = list(itertools.combinations(labels, 4))
    if len(quads) > max_quartets:
        step = len(quads) / float(max_quartets)
        quads = [quads[int(i * step)] for i in range(max_quartets)]
        labels_used = sorted({x for q in quads for x in q})
    else:
        labels_used = labels

    g_res = _quartet_resolution(gene, labels_used)
    r_res = _quartet_resolution(ref, labels_used)

    compared = 0
    conflict = 0
    for q in quads:
        gr, rr = g_res.get(q), r_res.get(q)
        if gr is None or rr is None:
            continue  # unresolved somewhere -> not counted
        compared += 1
        if gr != rr:
            conflict += 1
    frac = conflict / float(compared) if compared else 0.0
    return conflict, compared, frac


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene-tree", required=True,
                    help="collapsed gene tree (Newick)")
    ap.add_argument("--reference", required=True, help="reference tree (Newick)")
    ap.add_argument("--marker", required=True)
    ap.add_argument("--rf-max", type=float, default=0.35,
                    help="normalised RF ceiling (default 0.35)")
    ap.add_argument("--quartet-max", type=float, default=0.30,
                    help="normalised quartet ceiling (default 0.30)")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-json", default=None)
    args = ap.parse_args(argv)

    gene, ref, shared = load_common(args.gene_tree, args.reference)
    n = len(shared)

    if n < 4:
        # Too little overlap to compare -- report as failing so the aggressive
        # policy can drop it (a marker we cannot vet is a marker we cannot trust).
        row = OrderedDict([
            ("marker", args.marker), ("n_shared", n),
            ("rf", "NA"), ("rf_norm", "NA"),
            ("quartet", "NA"), ("quartet_norm", "NA"),
            ("rf_pass", False), ("quartet_pass", False), ("pass", False),
        ])
    else:
        rf, rf_norm = normalised_rf(gene, ref, n)
        qconf, qcomp, q_norm = quartet_distance(gene, ref, shared)
        rf_pass = bool(rf_norm <= args.rf_max)
        q_pass = bool(q_norm <= args.quartet_max)
        row = OrderedDict([
            ("marker", args.marker), ("n_shared", n),
            ("rf", rf), ("rf_norm", round(rf_norm, 4)),
            ("quartet", qconf), ("quartet_norm", round(q_norm, 4)),
            ("rf_pass", rf_pass), ("quartet_pass", q_pass),
            ("pass", bool(rf_pass and q_pass)),
        ])

    os.makedirs(os.path.dirname(os.path.abspath(args.out_tsv)), exist_ok=True)
    with open(args.out_tsv, "w") as fh:
        fh.write("\t".join(row.keys()) + "\n")
        fh.write("\t".join(str(v) for v in row.values()) + "\n")
    if args.out_json:
        with open(args.out_json, "w") as fh:
            json.dump(row, fh, indent=2)

    print("[tree_distances] {m}: RFnorm={rf} qnorm={q} -> {v}".format(
        m=args.marker, rf=row["rf_norm"], q=row["quartet_norm"],
        v="PASS" if row["pass"] else "FAIL"), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
