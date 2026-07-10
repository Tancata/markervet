# =============================================================================
# Stage 10 -- supertree (ASTER / ASTRAL)
# -----------------------------------------------------------------------------
# Summarise the passer gene trees into a single coalescent species tree with
# ASTER (https://github.com/chaoszhang/ASTER).  This is the coalescence
# counterpart to the concatenation supermatrix: it takes exactly the KEEP
# markers' collapsed gene trees (the checkpoint-selected set assembled by
# outputs.smk) and estimates the species tree under the multi-species coalescent.
#
# Tool is chosen in the regime YAML (supertree.tool):
#   astral      -- ASTRAL-IV, single-copy gene trees (the default here)
#   wastral     -- weighted ASTRAL (uses branch support/length weighting)
#   astral-pro  -- ASTRAL-Pro, for multi-copy (paralogous) gene trees
#
# Collapsed gene trees carry polytomies (weak branches were collapsed upstream);
# ASTER treats these as soft, and tolerates the differing taxon sets across
# markers, so no extra harmonisation is needed.
#
# Output:
#   {RESULTS}/outputs/supertree.nwk         coalescent species tree
#   {RESULTS}/outputs/supertree.input.nwk   the concatenated gene-tree input
# =============================================================================

_st = config.get("supertree", {})
ST_TOOL = _st.get("tool", "astral")


rule supertree:
    # keep_gene_trees (defined in outputs.smk) fans in over the KEEP set via the
    # aggregate_scorecard checkpoint, so the supertree is built from passers only.
    input:
        trees=keep_gene_trees,
    output:
        nwk=os.path.join(OUT_DIR, "supertree.nwk"),
        combined=os.path.join(OUT_DIR, "supertree.input.nwk"),
    params:
        tool=ST_TOOL,
        extra=_st.get("extra_args", ""),
    threads: _st.get("threads", 4)
    log:
        os.path.join(OUT_DIR, "supertree.log"),
    shell:
        r"""
        mkdir -p $(dirname {output.nwk})
        # One gene tree per line -> ASTER input.
        cat {input.trees} > {output.combined}
        if [ ! -s {output.combined} ]; then
            echo "[supertree] no KEEP gene trees -- nothing to summarise" > {log}
            : > {output.nwk}
        else
            {params.tool} -t {threads} {params.extra} \
                -o {output.nwk} {output.combined} > {log} 2>&1
        fi
        """
