# =============================================================================
# Stage 4 -- gene trees + support-based collapse
# -----------------------------------------------------------------------------
# Build one gene tree per marker with a site-heterogeneous mixture model
# (LG+C20/C60+G) -- the class of model that resists long-branch attraction at
# these depths -- then COLLAPSE every branch that fails the support gate, so that
# all downstream comparisons see signal, not estimation noise.
#
# PMSF (posterior mean site frequencies) is a fast approximation to the full
# CXX mixture: infer a guide tree under the mixture, then re-optimise under PMSF.
# When gene_trees.pmsf is true we run the two-step protocol; otherwise a single
# full-mixture search.  UFBoot + SH-aLRT give the two support values that the
# "aLRT/UFBoot" internal labels encode.
#
# Output:
#   {RESULTS}/genetrees/{marker}.treefile        ML gene tree (support-labelled)
#   {RESULTS}/genetrees/{marker}.collapsed.nwk    weak branches collapsed
# =============================================================================

GT_DIR = os.path.join(RESULTS, "genetrees")
_gt = config["gene_trees"]
_col = config["collapse"]

# Search / support flags for the single-step (non-PMSF) gene-tree rule.
#   fast: true  -> IQ-TREE quick search (-fast) with SH-aLRT support only.
#                  UFBoot is incompatible with -fast, so no -B; the lone SH-aLRT
#                  value is what the collapse / monophyly support gates read.
#   fast: false -> full search with UFBoot + SH-aLRT (labels "aLRT/UFBoot").
if _gt.get("fast", False):
    _SEARCH_FLAGS = "-fast -alrt %d" % _gt["alrt"]
else:
    _SEARCH_FLAGS = "-B %d -alrt %d" % (_gt["ufboot"], _gt["alrt"])


if _gt.get("pmsf", True):

    rule gene_tree_guide:
        # Step 1: a guide tree under the mixture model (no bootstrap needed).
        input:
            clean=os.path.join(CLEAN_DIR, "{marker}.clean.faa"),
        output:
            guide=os.path.join(GT_DIR, "guide", "{marker}.treefile"),
        params:
            model=_gt["model"],
            prefix=os.path.join(GT_DIR, "guide", "{marker}"),
        threads: _gt.get("threads", 4)
        log:
            os.path.join(GT_DIR, "guide", "{marker}.log"),
        shell:
            r"""
            mkdir -p $(dirname {output.guide})
            iqtree3 -s {input.clean} -m {params.model} -T {threads} \
                    --prefix {params.prefix} -redo > {log} 2>&1
            """

    rule gene_tree_pmsf:
        # Step 2: PMSF final tree, seeded by the guide tree (-ft), with UFBoot
        # + SH-aLRT support.
        input:
            clean=os.path.join(CLEAN_DIR, "{marker}.clean.faa"),
            guide=os.path.join(GT_DIR, "guide", "{marker}.treefile"),
        output:
            tree=os.path.join(GT_DIR, "{marker}.treefile"),
        params:
            model=_gt["model"],
            ufboot=_gt["ufboot"],
            alrt=_gt["alrt"],
            prefix=os.path.join(GT_DIR, "{marker}"),
        threads: _gt.get("threads", 4)
        log:
            os.path.join(GT_DIR, "{marker}.log"),
        shell:
            r"""
            mkdir -p $(dirname {output.tree})
            iqtree3 -s {input.clean} -m {params.model} -ft {input.guide} \
                    -B {params.ufboot} -alrt {params.alrt} -T {threads} \
                    --prefix {params.prefix} -redo > {log} 2>&1
            """

else:

    rule gene_tree_full:
        # Single-step full-mixture search (slower; no PMSF approximation).
        input:
            clean=os.path.join(CLEAN_DIR, "{marker}.clean.faa"),
        output:
            tree=os.path.join(GT_DIR, "{marker}.treefile"),
        params:
            model=_gt["model"],
            search=_SEARCH_FLAGS,
            prefix=os.path.join(GT_DIR, "{marker}"),
        threads: _gt.get("threads", 4)
        log:
            os.path.join(GT_DIR, "{marker}.log"),
        shell:
            r"""
            mkdir -p $(dirname {output.tree})
            iqtree3 -s {input.clean} -m {params.model} \
                    {params.search} -T {threads} \
                    --prefix {params.prefix} -redo > {log} 2>&1
            """


rule collapse_support:
    # Collapse internal branches failing the UFBoot/aLRT gate into polytomies.
    input:
        tree=os.path.join(GT_DIR, "{marker}.treefile"),
    output:
        collapsed=os.path.join(GT_DIR, "{marker}.collapsed.nwk"),
    params:
        ufboot_min=_col["ufboot_min"],
        alrt_min=_col["alrt_min"],
        logic=_col.get("logic", "either"),
        script=os.path.join(SCRIPTS, "collapse_support.py"),
    shell:
        r"""
        python {params.script} --tree {input.tree} --out {output.collapsed} \
            --marker {wildcards.marker} \
            --ufboot-min {params.ufboot_min} --alrt-min {params.alrt_min} \
            --logic {params.logic} --out-json {output.collapsed}.json
        """
