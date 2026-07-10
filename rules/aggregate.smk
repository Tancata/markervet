# =============================================================================
# Stage 8 -- aggregate scorecard
# -----------------------------------------------------------------------------
# Collect every per-marker screen and apply the AGGRESSIVE drop policy
# (scripts/aggregate_scorecard.py): a marker is dropped if it fails ANY enabled
# criterion.  Requesting the scorecard therefore pulls the entire per-marker DAG.
#
# Output:
#   {RESULTS}/aggregate/scorecard.tsv          decision + reasons per marker
#   {RESULTS}/aggregate/effective_config.yaml  the fully-resolved config used
#   {RESULTS}/aggregate/markers.txt            the marker list screened
# =============================================================================

AGG_DIR = os.path.join(RESULTS, "aggregate")


rule dump_effective_config:
    # Freeze the resolved config so the scorecard is reproducible and the
    # aggregator reads thresholds from exactly what this run used.
    output:
        yaml=os.path.join(AGG_DIR, "effective_config.yaml"),
    run:
        import yaml as _yaml
        os.makedirs(os.path.dirname(output.yaml), exist_ok=True)
        with open(output.yaml, "w") as fh:
            _yaml.safe_dump(dict(config), fh, sort_keys=False)


rule marker_list:
    output:
        txt=os.path.join(AGG_DIR, "markers.txt"),
    run:
        os.makedirs(os.path.dirname(output.txt), exist_ok=True)
        with open(output.txt, "w") as fh:
            fh.write("\n".join(MARKERS) + "\n")


# A checkpoint (not a plain rule): the set of KEEP markers is only known after
# the scorecard is written, and stage 9 fans out over exactly that set.
checkpoint aggregate_scorecard:
    input:
        cfg=os.path.join(AGG_DIR, "effective_config.yaml"),
        markers=os.path.join(AGG_DIR, "markers.txt"),
        # Fan-in of every per-marker screen (forces the whole pipeline to run).
        comp=expand(os.path.join(COMP_DIR, "{marker}.tsv"), marker=MARKERS),
        sat=expand(os.path.join(SAT_DIR, "{marker}.tsv"), marker=MARKERS),
        cong=expand(os.path.join(CONG_DIR, "{marker}.tsv"), marker=MARKERS),
        mono=expand(os.path.join(MONO_DIR, "{marker}.summary.tsv"), marker=MARKERS),
        # Reconciliation is an optional vote: only require its per-marker outputs
        # (and thus run the reconciliation rule) when it is enabled in the regime
        # YAML.  Disabled => empty input, so the vote drops out of the DAG cleanly.
        recon=(expand(os.path.join(RECON_DIR, "{marker}.tsv"), marker=MARKERS)
               if config["outliers_hgt"]["reconciliation"]["enabled"] else []),
        phylter=os.path.join(PHYLTER_DIR, "summary.tsv"),
    output:
        scorecard=os.path.join(AGG_DIR, "scorecard.tsv"),
    params:
        script=os.path.join(SCRIPTS, "aggregate_scorecard.py"),
    shell:
        r"""
        python {params.script} --config {input.cfg} --markers {input.markers} \
            --comp-dir {COMP_DIR} --sat-dir {SAT_DIR} --cong-dir {CONG_DIR} \
            --mono-dir {MONO_DIR} --recon-dir {RECON_DIR} \
            --phylter-summary {input.phylter} --out {output.scorecard}
        """
