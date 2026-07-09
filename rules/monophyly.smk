# =============================================================================
# Stage 6 -- monophyly recovery  (headline metric + filter)
# -----------------------------------------------------------------------------
# For each named clade in the regime config (d__Archaea, d__Bacteria, Asgard,
# DPANN, Eukaryota, ...), score how faithfully this marker recovers it:
# intruders, outliers, and a per-clade recovery fraction, support-gated
# (scripts/monophyly_recovery.py, MonoPhy-style).  The per-marker mean is the
# clade_recovery_score used by the drop policy; the long per-clade table is
# assembled into the clade x marker matrix by the report stage.
#
# Output:
#   {RESULTS}/monophyly/{marker}.matrix.tsv     one row per clade (long form)
#   {RESULTS}/monophyly/{marker}.summary.tsv     one-row per-marker headline
#   {RESULTS}/monophyly/clades.json              clade definitions (from config)
# =============================================================================

MONO_DIR = os.path.join(RESULTS, "monophyly")
_mono = config["monophyly"]


rule monophyly_clades_json:
    # Serialise the clade definitions from config so the Python script can load
    # them without re-parsing YAML.
    output:
        json=os.path.join(MONO_DIR, "clades.json"),
    run:
        import json as _json
        os.makedirs(os.path.dirname(output.json), exist_ok=True)
        with open(output.json, "w") as fh:
            _json.dump(_mono["clades"], fh, indent=2)


rule monophyly_recovery:
    input:
        tree=os.path.join(GT_DIR, "{marker}.collapsed.nwk"),
        clades=os.path.join(MONO_DIR, "clades.json"),
    output:
        matrix=os.path.join(MONO_DIR, "{marker}.matrix.tsv"),
        summary=os.path.join(MONO_DIR, "{marker}.summary.tsv"),
    params:
        support_gate=_mono["support_gate"],
        recovery_min=_mono["recovery_min"],
        # Emit the --taxon-table flag only when a table is configured.
        tt_flag=lambda wc: ("--taxon-table " + config["taxon_table"])
        if config.get("taxon_table") else "",
        script=os.path.join(SCRIPTS, "monophyly_recovery.py"),
    shell:
        r"""
        python {params.script} --tree {input.tree} --marker {wildcards.marker} \
            --clades-json {input.clades} {params.tt_flag} \
            --support-gate {params.support_gate} \
            --recovery-min {params.recovery_min} \
            --out-matrix {output.matrix} --out-summary {output.summary} \
            --out-json {output.summary}.json
        """
