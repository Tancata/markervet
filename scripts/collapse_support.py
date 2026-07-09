#!/usr/bin/env python3
"""Collapse weakly-supported internal branches of a gene tree (MarkerVet stage 4).

*The* central design principle of MarkerVet: at domain-level depths a single
gene tree is mostly estimation error, so comparing raw gene trees measures noise,
not conflict.  Before any tree-vs-tree comparison (congruence, monophyly,
PhylteR, reconciliation) we therefore collapse every internal branch whose
support is too low into a polytomy.  What survives is only the signal the data
actually support.

IQ-TREE writes combined support as an internal-node label of the form

    "aLRT/UFBoot"     e.g.  "83.4/97"

(SH-aLRT first, UFBoot second) when run with `-alrt ... -bb ...`.  This script
parses that label and collapses a branch according to `--logic`:

    either : collapse if aLRT < alrt_min  OR  UFBoot < ufboot_min   (strict)
    both   : collapse only if aLRT < alrt_min AND UFBoot < ufboot_min (lenient)

Collapsing sets the branch length to zero and splices the node's children into
its parent (ete3's `delete`), producing a multifurcating tree.  Terminal
branches and the root are never touched.

Also emits a tiny JSON summary (branches examined / collapsed / retained,
mean retained support) that feeds the report.
"""
from __future__ import print_function

import argparse
import json
import os
import re
import sys

from ete3 import Tree

# Matches "aLRT/UFBoot", "aLRT", or a bare number, tolerating surrounding
# whitespace.  Groups: (first_value, optional_second_value).
_SUPPORT_RE = re.compile(r"^\s*([0-9]*\.?[0-9]+)\s*(?:/\s*([0-9]*\.?[0-9]+)\s*)?$")


def parse_support(label):
    """Parse an IQ-TREE internal label into (alrt, ufboot).

    Returns (None, None) when the label is empty or unparseable (e.g. a named
    internal node), signalling "no support info -> do not collapse".
    """
    if label is None:
        return None, None
    m = _SUPPORT_RE.match(str(label))
    if not m:
        return None, None
    first = float(m.group(1))
    second = float(m.group(2)) if m.group(2) is not None else None
    if second is None:
        # Only one value present.  By convention we treat a lone value as
        # UFBoot (the more commonly reported single support); aLRT unknown.
        return None, first
    return first, second


def should_collapse(alrt, ufboot, alrt_min, ufboot_min, logic):
    """Decide whether a branch fails the support gate.

    A gate is only applied when its value is available; if support could not be
    parsed at all the branch is conservatively retained.
    """
    fail_alrt = alrt is not None and alrt < alrt_min
    fail_uf = ufboot is not None and ufboot < ufboot_min
    if alrt is None and ufboot is None:
        return False
    if logic == "both":
        # Lenient: collapse only when *both* available gates fail.  If only one
        # gate is available, fall back to that single gate.
        gates = [g for g in (fail_alrt if alrt is not None else None,
                             fail_uf if ufboot is not None else None)
                 if g is not None]
        return all(gates) if gates else False
    # Default "either": collapse if any available gate fails (strict).
    return fail_alrt or fail_uf


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree", required=True, help="input gene tree (Newick)")
    ap.add_argument("--out", required=True, help="collapsed tree (Newick)")
    ap.add_argument("--marker", default="?")
    ap.add_argument("--ufboot-min", type=float, default=95.0)
    ap.add_argument("--alrt-min", type=float, default=80.0)
    ap.add_argument("--logic", choices=["either", "both"], default="either")
    ap.add_argument("--out-json", default=None)
    args = ap.parse_args(argv)

    # format=1 keeps internal-node names/labels (that is where IQ-TREE stores
    # the aLRT/UFBoot string).
    tree = Tree(args.tree, format=1)

    examined = collapsed = 0
    retained_support = []

    # Collect first, then delete: mutating the tree while traversing is unsafe.
    to_collapse = []
    for node in tree.traverse("postorder"):
        if node.is_leaf() or node.is_root():
            continue
        examined += 1
        alrt, ufboot = parse_support(node.name)
        if should_collapse(alrt, ufboot, args.alrt_min, args.ufboot_min, args.logic):
            to_collapse.append(node)
        else:
            # Record the UFBoot (or aLRT) we kept, for reporting.
            keep = ufboot if ufboot is not None else alrt
            if keep is not None:
                retained_support.append(keep)

    for node in to_collapse:
        node.delete(prevent_nondicotomic=False, preserve_branch_length=True)
        collapsed += 1

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    # format=1 so we round-trip any surviving internal labels.
    tree.write(outfile=args.out, format=1)

    summary = {
        "marker": args.marker,
        "internal_branches": examined,
        "collapsed": collapsed,
        "retained": examined - collapsed,
        "ufboot_min": args.ufboot_min,
        "alrt_min": args.alrt_min,
        "logic": args.logic,
        "mean_retained_support": (sum(retained_support) / len(retained_support))
        if retained_support else None,
    }
    if args.out_json:
        with open(args.out_json, "w") as fh:
            json.dump(summary, fh, indent=2)

    print("[collapse_support] {m}: collapsed {c}/{e} internal branches".format(
        m=args.marker, c=collapsed, e=examined), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
