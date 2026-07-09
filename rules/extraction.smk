# =============================================================================
# Stage 1 -- extraction
# -----------------------------------------------------------------------------
# Produce one unaligned per-marker protein FASTA:  {RESULTS}/markers/{marker}.faa
#
# Two modes, chosen by config["markers"]["extract"]:
#
#   extract: false  -- the user already has per-marker FASTAs.  We simply stage
#                      (copy) them into the results tree so every downstream rule
#                      has a single, uniform input path.  The hmmsearch machinery
#                      below is never invoked.
#
#   extract: true   -- run hmmsearch with each marker's HMM against every
#                      proteome, enforce single-copy (one sequence per taxon,
#                      the best-scoring hit, dropped if a second hit scores
#                      within `paralog_margin` bits => ambiguous paralogy), and
#                      assemble the per-marker FASTA.
#
# The single-copy logic is real Python (in the `run:` block); hmmsearch itself is
# the external tool and carries a runnable shell template.
# =============================================================================

MARKERS_DIR = os.path.join(RESULTS, "markers")
_mcfg = config["markers"]

if not _mcfg.get("extract", False):

    # ---- supplied per-marker FASTAs: just stage them -----------------------
    rule stage_marker_fasta:
        input:
            faa=lambda wc: os.path.join(_mcfg["fasta_dir"],
                                        wc.marker + _mcfg.get("fasta_ext", ".faa")),
        output:
            faa=os.path.join(MARKERS_DIR, "{marker}.faa"),
        shell:
            # cp (not symlink) so the results dir is self-contained/portable.
            "mkdir -p $(dirname {output.faa}) && cp {input.faa} {output.faa}"

else:

    # ---- de novo extraction from proteomes via HMMER -----------------------
    HMM_DIR = _mcfg["hmm_dir"]
    PROTEOMES_DIR = _mcfg["proteomes_dir"]

    def proteome_files():
        return sorted(glob.glob(os.path.join(PROTEOMES_DIR, "*.faa")))

    rule hmmsearch_marker:
        # Search one marker HMM against the concatenation of all proteomes.
        # --domtblout gives per-domain hits we parse for scores/coordinates.
        input:
            hmm=lambda wc: os.path.join(HMM_DIR, wc.marker + ".hmm"),
            proteomes=proteome_files(),
        output:
            domtbl=os.path.join(RESULTS, "extraction", "hmmsearch",
                                "{marker}.domtblout"),
        params:
            evalue=config["markers"].get("hmmsearch_evalue", "1e-10"),
        threads: 2
        shell:
            r"""
            mkdir -p $(dirname {output.domtbl})
            # Concatenate proteomes on the fly so a single hmmsearch sees them all.
            cat {input.proteomes} | \
            hmmsearch --cpu {threads} -E {params.evalue} \
                      --domtblout {output.domtbl} {input.hmm} - > /dev/null
            """

    rule build_marker_fasta:
        # Enforce single-copy and extract the winning sequences.  Pure-Python
        # logic: best hit per taxon, dropped when a paralog scores within
        # `paralog_margin` bits of the best (ambiguous single-copy).
        input:
            domtbl=os.path.join(RESULTS, "extraction", "hmmsearch",
                                "{marker}.domtblout"),
            proteomes=proteome_files(),
        output:
            faa=os.path.join(MARKERS_DIR, "{marker}.faa"),
        params:
            margin=config["markers"].get("paralog_margin", 20.0),
        run:
            import re as _re
            from collections import defaultdict

            # --- map every target sequence id to its best bit-score ---------
            best = {}          # seqid -> best domain bitscore
            for line in open(input.domtbl):
                if line.startswith("#"):
                    continue
                f = line.split()
                if len(f) < 14:
                    continue
                seqid = f[0]
                score = float(f[13])          # domain-level bit score
                if seqid not in best or score > best[seqid]:
                    best[seqid] = score

            # --- taxon of each hit ------------------------------------------
            # Sequence headers are assumed prefixed "taxon|...": adjust the
            # split here to match your proteome header convention.
            def taxon_of(seqid):
                return seqid.split("|")[0]

            by_taxon = defaultdict(list)
            for seqid, score in best.items():
                by_taxon[taxon_of(seqid)].append((score, seqid))

            # --- single-copy decision per taxon -----------------------------
            keep = {}          # taxon -> winning seqid
            for taxon, hits in by_taxon.items():
                hits.sort(reverse=True)      # highest score first
                top_score, top_id = hits[0]
                if len(hits) > 1 and (top_score - hits[1][0]) < params.margin:
                    # Runner-up too close: ambiguous paralogy, skip this taxon.
                    continue
                keep[taxon] = top_id

            wanted = set(keep.values())

            # --- pull the sequences out of the proteomes --------------------
            seqs = {}
            for pf in input.proteomes:
                name, chunks = None, []
                for ln in open(pf):
                    ln = ln.rstrip("\n")
                    if ln.startswith(">"):
                        if name in wanted:
                            seqs[name] = "".join(chunks)
                        name = ln[1:].split()[0]
                        chunks = []
                    else:
                        chunks.append(ln.strip())
                if name in wanted:
                    seqs[name] = "".join(chunks)

            os.makedirs(os.path.dirname(output.faa), exist_ok=True)
            with open(output.faa, "w") as out:
                for taxon, seqid in sorted(keep.items()):
                    if seqid in seqs:
                        # Relabel to the taxon so all markers share leaf names.
                        out.write(">%s\n%s\n" % (taxon, seqs[seqid]))
