# =============================================================================
# Stage 7 -- outlier / HGT detection  (two independent votes)
# -----------------------------------------------------------------------------
# Only reached AFTER the artefact screen, so composition/saturation artefacts
# are not mistaken for transfer.  Two orthogonal signals:
#
#   PhylteR (scripts/phylter_run.R) -- ONE DISTATIS run over ALL collapsed gene
#     trees; flags gene x taxon outlier cells and whole-gene outliers.  Broad,
#     distance-based, topology-agnostic.
#
#   Reconciliation (GeneRax or ALE) -- per-marker gene-tree/species-tree
#     reconciliation; counts inferred transfers, localising HGT.  A marker with
#     excess transfers is a transfer-laden marker.
#
# Output:
#   {RESULTS}/outliers_hgt/phylter/summary.tsv    per-marker PhylteR summary
#   {RESULTS}/outliers_hgt/phylter/cells.tsv      removed (marker, taxon) cells
#   {RESULTS}/outliers_hgt/recon/{marker}.tsv     per-marker transfer count
# =============================================================================

HGT_DIR = os.path.join(RESULTS, "outliers_hgt")
PHYLTER_DIR = os.path.join(HGT_DIR, "phylter")
RECON_DIR = os.path.join(HGT_DIR, "recon")
_hgt = config["outliers_hgt"]
REFERENCE = config["reference_tree"]


rule phylter:
    # PhylteR runs once across the whole marker set (it needs the full
    # gene x taxon block to define the consensus).  Requesting the summary pulls
    # every marker's collapsed tree.
    input:
        trees=expand(os.path.join(GT_DIR, "{marker}.collapsed.nwk"),
                     marker=MARKERS),
    output:
        summary=os.path.join(PHYLTER_DIR, "summary.tsv"),
        cells=os.path.join(PHYLTER_DIR, "cells.tsv"),
    params:
        tree_dir=GT_DIR,
        script=os.path.join(SCRIPTS, "phylter_run.R"),
    log:
        os.path.join(PHYLTER_DIR, "phylter.log"),
    shell:
        r"""
        mkdir -p {PHYLTER_DIR}
        Rscript {params.script} --tree-dir {params.tree_dir} \
            --tree-ext .collapsed.nwk \
            --out-cells {output.cells} --out-summary {output.summary} \
            > {log} 2>&1
        """


if _hgt["reconciliation"]["tool"] == "ale":

    rule reconciliation_ale:
        # ALE undated reconciliation.  ALEobserve builds the .ale from a bootstrap
        # sample; ALEml_undated reconciles against the species tree and reports a
        # "# of Transfers" line we parse.
        input:
            gene=os.path.join(GT_DIR, "{marker}.treefile"),
            species=REFERENCE,
        output:
            tsv=os.path.join(RECON_DIR, "{marker}.tsv"),
        params:
            outdir=os.path.join(RECON_DIR, "{marker}"),
        log:
            os.path.join(RECON_DIR, "{marker}.log"),
        shell:
            r"""
            mkdir -p {params.outdir} $(dirname {output.tsv})
            cp {input.gene} {params.outdir}/gene.treefile
            ( cd {params.outdir} && \
              ALEobserve gene.treefile && \
              ALEml_undated {input.species} gene.treefile.ale ) > {log} 2>&1
            # ALE writes "<...>.uml_rec" carrying a "# of Transfers :" line.
            transfers=$(grep -m1 -i "Transfers" {params.outdir}/*.uml_rec \
                        | grep -oE "[0-9]+(\.[0-9]+)?" | head -n1)
            [ -z "$transfers" ] && transfers=NA
            printf "marker\ttransfers\n{wildcards.marker}\t$transfers\n" \
                > {output.tsv}
            """

else:

    rule reconciliation_generax:
        # GeneRax UndatedDTL reconciliation.  Leaves are already taxon-labelled,
        # so the gene->species mapping is the identity.  Transfers are the H
        # (horizontal) events annotated in the reconciled NHX gene tree.
        input:
            gene=os.path.join(GT_DIR, "{marker}.collapsed.nwk"),
            aln=os.path.join(CLEAN_DIR, "{marker}.clean.faa"),
            species=REFERENCE,
        output:
            tsv=os.path.join(RECON_DIR, "{marker}.tsv"),
        params:
            outdir=os.path.join(RECON_DIR, "{marker}"),
            model=config["gene_trees"]["model"],
        log:
            os.path.join(RECON_DIR, "{marker}.log"),
        shell:
            r"""
            mkdir -p {params.outdir} $(dirname {output.tsv})
            # Identity gene->species mapping (leaf name == taxon name).
            grep "^>" {input.aln} | sed 's/^>//' | awk '{{print $1" "$1}}' \
                > {params.outdir}/mapping.link
            # GeneRax family description file.
            cat > {params.outdir}/families.txt <<EOF
[FAMILIES]
- {wildcards.marker}
starting_gene_tree = {input.gene}
alignment = {input.aln}
mapping = {params.outdir}/mapping.link
subst_model = {params.model}
EOF
            generax -f {params.outdir}/families.txt -s {input.species} \
                    -r UndatedDTL -p {params.outdir}/run > {log} 2>&1 || true
            # Count transfer (H) events in the reconciled NHX gene tree.
            rec=$(find {params.outdir}/run -name "*.newick" 2>/dev/null | head -n1)
            if [ -n "$rec" ]; then
                transfers=$(grep -o "H=Y" "$rec" | wc -l | tr -d ' ')
            else
                transfers=NA
            fi
            printf "marker\ttransfers\n{wildcards.marker}\t$transfers\n" \
                > {output.tsv}
            """
