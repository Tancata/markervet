# =============================================================================
# Stage 5 -- congruence vs the reference
# -----------------------------------------------------------------------------
# Two complementary comparisons of each collapsed gene tree against the trusted
# reference topology:
#
#   * IQ-TREE gCF/sCF -- gene- and site-concordance factors: how often the
#     reference's branches are recovered by this marker's tree/sites.  Reported
#     for the record (and an optional gcf_min gate).
#   * scripts/tree_distances.py -- normalised Robinson-Foulds and quartet
#     distance on the shared taxon set (dendropy).  These drive the drop policy.
#
# Because weak branches were already collapsed (stage 4), any residual conflict
# is real, not noise.
#
# Output:
#   {RESULTS}/congruence/{marker}.tsv            RF + quartet + pass/fail
#   {RESULTS}/congruence/{marker}.cf.stat        IQ-TREE concordance factors
# =============================================================================

CONG_DIR = os.path.join(RESULTS, "congruence")
_cong = config["congruence"]
REFERENCE = config["reference_tree"]


rule concordance_factors:
    # gCF/sCF of the reference branches given this marker's gene tree + alignment.
    # (Single-gene gCF is coarse but sCF from the alignment is informative.)
    input:
        ref=REFERENCE,
        gene=os.path.join(GT_DIR, "{marker}.collapsed.nwk"),
        aln=os.path.join(CLEAN_DIR, "{marker}.clean.faa"),
    output:
        stat=os.path.join(CONG_DIR, "{marker}.cf.stat"),
    params:
        prefix=os.path.join(CONG_DIR, "{marker}.cf"),
    threads: 2
    log:
        os.path.join(CONG_DIR, "{marker}.cf.log"),
    shell:
        r"""
        mkdir -p $(dirname {output.stat})
        # --gcf takes the gene tree(s); --scf draws quartets from the alignment.
        iqtree3 -t {input.ref} --gcf {input.gene} -s {input.aln} \
                --scf 100 -T {threads} --prefix {params.prefix} -redo \
                > {log} 2>&1
        # Guarantee the declared output exists even for degenerate inputs.
        [ -f {params.prefix}.cf.stat ] && cp {params.prefix}.cf.stat {output.stat} \
            || echo "# no concordance factors computed" > {output.stat}
        """


rule tree_distances:
    # Normalised RF + quartet distance vs the reference (the drop-policy inputs).
    input:
        gene=os.path.join(GT_DIR, "{marker}.collapsed.nwk"),
        ref=REFERENCE,
    output:
        tsv=os.path.join(CONG_DIR, "{marker}.tsv"),
    params:
        rf_max=_cong["rf_max"],
        quartet_max=_cong["quartet_max"],
        script=os.path.join(SCRIPTS, "tree_distances.py"),
    shell:
        r"""
        python {params.script} --gene-tree {input.gene} --reference {input.ref} \
            --marker {wildcards.marker} --rf-max {params.rf_max} \
            --quartet-max {params.quartet_max} \
            --out-tsv {output.tsv} --out-json {output.tsv}.json
        """
