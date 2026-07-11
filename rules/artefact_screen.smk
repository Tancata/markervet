# =============================================================================
# Stage 3 -- artefact screen  (runs BEFORE HGT calling)
# -----------------------------------------------------------------------------
# The two eukaryote-placement failure modes -- compositional attraction and
# substitution saturation -- are screened here, together with rogue-tip removal,
# so that they are never mislabelled as horizontal transfer in stage 7.
#
#   quicktree     : a fast ML tree from the trimmed alignment, used only as
#                   scaffolding for TreeShrink and the saturation slope (the
#                   final, model-rich gene tree is built in stage 4).
#   treeshrink    : remove rogue tips (abnormally long root-to-tip paths).
#   composition   : chi-square amino-acid homogeneity     (scripts/compositional_test.py)
#   saturation    : p-distance ~ patristic slope          (scripts/saturation.py)
#   recode        : optional AA recoding track (Dayhoff-6 / SR-4 / SR-6),
#                   enabled for the three-domain regime to blunt composition.
#
# The cleaned (post-TreeShrink) trimmed alignment,
#   {RESULTS}/artefact/clean/{marker}.clean.faa,
# is the canonical alignment consumed by every later stage.
# =============================================================================

ART_DIR = os.path.join(RESULTS, "artefact")
CLEAN_DIR = os.path.join(ART_DIR, "clean")
COMP_DIR = os.path.join(ART_DIR, "composition")
SAT_DIR = os.path.join(ART_DIR, "saturation")

_art = config["artefact"]


rule quicktree:
    # Fast ML tree (IQ-TREE -fast, LG+G) purely as scaffolding for TreeShrink
    # and the saturation regression; not used for any topological inference.
    input:
        trimmed=os.path.join(TRIM_DIR, "{marker}.trimmed.faa"),
    output:
        tree=os.path.join(ART_DIR, "quicktree", "{marker}.treefile"),
    params:
        prefix=os.path.join(ART_DIR, "quicktree", "{marker}"),
    threads: 2
    log:
        os.path.join(ART_DIR, "quicktree", "{marker}.log"),
    shell:
        r"""
        mkdir -p $(dirname {output.tree})
        iqtree3 -s {input.trimmed} -m LG+G -fast -T {threads} \
                --prefix {params.prefix} -redo > {log} 2>&1
        """


if _art["treeshrink"]["enabled"]:

    rule treeshrink:
        # TreeShrink flags tips whose removal disproportionately shortens the
        # tree (rogue / mislabelled / contaminant sequences) and drops them from
        # both the tree and the alignment.
        input:
            trimmed=os.path.join(TRIM_DIR, "{marker}.trimmed.faa"),
            tree=os.path.join(ART_DIR, "quicktree", "{marker}.treefile"),
        output:
            clean=os.path.join(CLEAN_DIR, "{marker}.clean.faa"),
            removed=os.path.join(CLEAN_DIR, "{marker}.removed.txt"),
        params:
            quantile=_art["treeshrink"]["quantile"],
            outdir=os.path.join(ART_DIR, "treeshrink", "{marker}"),
        log:
            os.path.join(ART_DIR, "treeshrink", "{marker}.log"),
        shell:
            r"""
            mkdir -p {params.outdir} $(dirname {output.clean})
            run_treeshrink.py -t {input.tree} -a {input.trimmed} \
                -q {params.quantile} -O out -o {params.outdir} > {log} 2>&1
            # TreeShrink writes {{outdir}}/out.txt: whitespace-separated removed
            # tip labels (one line per gene).  In single-gene mode it does NOT
            # reliably emit the shrunk alignment, so we derive the cleaned FASTA
            # ourselves by dropping the removed tips from the trimmed alignment.
            ( [ -f {params.outdir}/out.txt ] && tr '\t' '\n' \
                < {params.outdir}/out.txt | sed '/^$/d' > {output.removed} ) \
                || : > {output.removed}
            # Read the removed-tip set in BEGIN via getline (robust when the
            # removed list is EMPTY -- the NR==FNR two-file idiom breaks there,
            # consuming the whole alignment as "removed" and emitting nothing).
            awk -v rmf={output.removed} \
                'BEGIN{{while((getline l < rmf)>0){{if(l!="") rm[l]=1}}}} \
                 /^>/{{k=substr($1,2); drop=(k in rm)}} !drop' \
                {input.trimmed} > {output.clean}
            """

else:

    rule treeshrink_skip:
        # TreeShrink disabled: the cleaned alignment is just the trimmed one.
        input:
            trimmed=os.path.join(TRIM_DIR, "{marker}.trimmed.faa"),
        output:
            clean=os.path.join(CLEAN_DIR, "{marker}.clean.faa"),
            removed=os.path.join(CLEAN_DIR, "{marker}.removed.txt"),
        shell:
            "mkdir -p $(dirname {output.clean}) && "
            "cp {input.trimmed} {output.clean} && : > {output.removed}"


rule composition_test:
    # Chi-square test of amino-acid compositional homogeneity across taxa.
    input:
        clean=os.path.join(CLEAN_DIR, "{marker}.clean.faa"),
    output:
        tsv=os.path.join(COMP_DIR, "{marker}.tsv"),
    params:
        alpha=_art["compositional"]["alpha"],
        script=os.path.join(SCRIPTS, "compositional_test.py"),
    shell:
        r"""
        python {params.script} --alignment {input.clean} --marker {wildcards.marker} \
            --alpha {params.alpha} --out-tsv {output.tsv} \
            --out-json {output.tsv}.json
        """


rule saturation_test:
    # Saturation via the slope of observed p-distance on patristic distance.
    input:
        clean=os.path.join(CLEAN_DIR, "{marker}.clean.faa"),
        tree=os.path.join(ART_DIR, "quicktree", "{marker}.treefile"),
    output:
        tsv=os.path.join(SAT_DIR, "{marker}.tsv"),
    params:
        min_slope=_art["saturation"]["min_slope"],
        script=os.path.join(SCRIPTS, "saturation.py"),
    shell:
        r"""
        python {params.script} --alignment {input.clean} --tree {input.tree} \
            --marker {wildcards.marker} --min-slope {params.min_slope} \
            --out-tsv {output.tsv} --out-json {output.tsv}.json
        """


# ---- optional amino-acid recoding track ------------------------------------
# Recoding lumps the 20 amino acids into a handful of exchange-rate classes,
# which strongly reduces compositional/saturation artefacts at the cost of
# resolution -- the right trade for the three-domain regime.  The recoded
# alignment is emitted as an alternative concatenation track in stage 9.
_RECODE_MAPS = {
    # Dayhoff-6 classes (Hrdy et al.): STPAG / DENQ / HKR / MIVL / FYW / C
    "dayhoff6": {
        **{c: "0" for c in "STPAG"}, **{c: "1" for c in "DENQ"},
        **{c: "2" for c in "HKR"}, **{c: "3" for c in "MIVL"},
        **{c: "4" for c in "FYW"}, **{c: "5" for c in "C"}},
    # Susko-Roger 4-state (SR4): AGNPST / CHWY / DEKQR / FILMV
    "sr4": {
        **{c: "0" for c in "AGNPST"}, **{c: "1" for c in "CHWY"},
        **{c: "2" for c in "DEKQR"}, **{c: "3" for c in "FILMV"}},
    # Susko-Roger 6-state (SR6): APST / DENG / QKR / MIVL / WC / FYH
    "sr6": {
        **{c: "0" for c in "APST"}, **{c: "1" for c in "DENG"},
        **{c: "2" for c in "QKR"}, **{c: "3" for c in "MIVL"},
        **{c: "4" for c in "WC"}, **{c: "5" for c in "FYH"}},
}

if config.get("recoding", {}).get("enabled", False):

    rule recode_alignment:
        # Recode each residue to its exchange class; gaps/ambiguities preserved.
        input:
            clean=os.path.join(CLEAN_DIR, "{marker}.clean.faa"),
        output:
            recoded=os.path.join(ART_DIR, "recoded", "{marker}.recoded.faa"),
        params:
            scheme=config["recoding"]["scheme"],
        run:
            table = _RECODE_MAPS[params.scheme]
            os.makedirs(os.path.dirname(output.recoded), exist_ok=True)
            with open(input.clean) as ih, open(output.recoded, "w") as oh:
                for line in ih:
                    if line.startswith(">"):
                        oh.write(line)
                    else:
                        oh.write("".join(table.get(c.upper(), c if c in "-.?" else "-")
                                         for c in line.strip()) + "\n")
