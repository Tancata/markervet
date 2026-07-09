# =============================================================================
# Stage 2 -- alignment + trimming
# -----------------------------------------------------------------------------
#   align :  MAFFT   (--auto, or L-INS-i for the accurate/slow track)
#   trim  :  ClipKIT (smart-gap)  OR  BMGE (entropy)
#
# Output:
#   {RESULTS}/align/{marker}.aln.faa       aligned protein FASTA
#   {RESULTS}/trim/{marker}.trimmed.faa    trimmed alignment (input to stage 3)
#
# Both tools are external; the shell blocks are runnable templates.
# =============================================================================

ALIGN_DIR = os.path.join(RESULTS, "align")
TRIM_DIR = os.path.join(RESULTS, "trim")

_acfg = config["align"]
_tcfg = config["trim"]


rule mafft_align:
    input:
        faa=os.path.join(MARKERS_DIR, "{marker}.faa"),
    output:
        aln=os.path.join(ALIGN_DIR, "{marker}.aln.faa"),
    params:
        # auto  -> `mafft --auto` (picks a strategy by size)
        # linsi -> the accurate L-INS-i strategy, preferred for few/deep markers
        mode=_acfg.get("mode", "auto"),
    threads: _acfg.get("threads", 4)
    log:
        os.path.join(ALIGN_DIR, "{marker}.mafft.log"),
    shell:
        r"""
        mkdir -p $(dirname {output.aln})
        if [ "{params.mode}" = "linsi" ]; then
            mafft --thread {threads} --localpair --maxiterate 1000 \
                  {input.faa} > {output.aln} 2> {log}
        else
            mafft --thread {threads} --auto {input.faa} > {output.aln} 2> {log}
        fi
        """


if _tcfg.get("tool", "clipkit") == "bmge":

    rule trim_bmge:
        # BMGE: entropy-based trimming; -h sets the smoothed-entropy cutoff
        # (lower = stricter, removes more fast/noisy columns).
        input:
            aln=os.path.join(ALIGN_DIR, "{marker}.aln.faa"),
        output:
            trimmed=os.path.join(TRIM_DIR, "{marker}.trimmed.faa"),
        params:
            entropy=_tcfg.get("bmge_entropy", 0.5),
        log:
            os.path.join(TRIM_DIR, "{marker}.bmge.log"),
        shell:
            r"""
            mkdir -p $(dirname {output.trimmed})
            bmge -i {input.aln} -t AA -h {params.entropy} \
                 -of {output.trimmed} > {log} 2>&1
            """

else:

    rule trim_clipkit:
        # ClipKIT smart-gap: keeps parsimony/phylogenetically-informative sites
        # while removing gappy columns; a strong default for deep alignments.
        input:
            aln=os.path.join(ALIGN_DIR, "{marker}.aln.faa"),
        output:
            trimmed=os.path.join(TRIM_DIR, "{marker}.trimmed.faa"),
        params:
            mode=_tcfg.get("clipkit_mode", "smart-gap"),
        log:
            os.path.join(TRIM_DIR, "{marker}.clipkit.log"),
        shell:
            r"""
            mkdir -p $(dirname {output.trimmed})
            clipkit {input.aln} -m {params.mode} \
                    -o {output.trimmed} > {log} 2>&1
            """
