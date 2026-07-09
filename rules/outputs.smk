# =============================================================================
# Stage 9 -- outputs
# -----------------------------------------------------------------------------
# Assemble the deliverables for the markers the scorecard chose to KEEP:
#
#   * cleaned/trimmed alignments of passers  -> outputs/alignments/{marker}.faa
#   * per-marker gene trees (coalescence in) -> outputs/gene_trees/{marker}.nwk
#   * concatenated supermatrix + partitions  -> outputs/supermatrix.{fasta,partitions}
#   * HTML + markdown report                 -> report/markervet_report.{html,md}
#
# The KEEP set is read from the aggregate_scorecard *checkpoint*, so this stage
# fans out over exactly the markers that passed.
# =============================================================================

OUT_DIR = os.path.join(RESULTS, "outputs")
REPORT_DIR = os.path.join(RESULTS, "report")


def keep_markers(wildcards):
    """Parse the scorecard checkpoint and return the list of KEEP markers."""
    sc = checkpoints.aggregate_scorecard.get(**wildcards).output.scorecard
    keep = []
    with open(sc) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        di = header.index("decision")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > di and f[di] == "KEEP":
                keep.append(f[0])
    return keep


# ---- passer alignments ------------------------------------------------------
rule passer_alignment:
    input:
        clean=os.path.join(CLEAN_DIR, "{marker}.clean.faa"),
    output:
        faa=os.path.join(OUT_DIR, "alignments", "{marker}.faa"),
    shell:
        "mkdir -p $(dirname {output.faa}) && cp {input.clean} {output.faa}"


# ---- passer gene trees (coalescence input) ----------------------------------
rule passer_gene_tree:
    input:
        tree=os.path.join(GT_DIR, "{marker}.collapsed.nwk"),
    output:
        nwk=os.path.join(OUT_DIR, "gene_trees", "{marker}.nwk"),
    shell:
        "mkdir -p $(dirname {output.nwk}) && cp {input.tree} {output.nwk}"


def keep_alignments(wildcards):
    return expand(os.path.join(OUT_DIR, "alignments", "{marker}.faa"),
                  marker=keep_markers(wildcards))


def keep_gene_trees(wildcards):
    return expand(os.path.join(OUT_DIR, "gene_trees", "{marker}.nwk"),
                  marker=keep_markers(wildcards))


# A sentinel that pulls every passer's gene tree (the coalescence input set).
rule gene_trees_done:
    input:
        keep_gene_trees,
    output:
        done=os.path.join(OUT_DIR, "gene_trees.done"),
    run:
        with open(output.done, "w") as fh:
            fh.write("\n".join(os.path.basename(f) for f in input) + "\n")


# ---- concatenated supermatrix + partitions ----------------------------------
rule supermatrix:
    # Real concatenation logic: union the taxa across passer alignments, pad
    # each marker's missing taxa with gaps, concatenate in marker order, and
    # emit a RAxML/IQ-TREE-style partitions file recording each marker's span.
    input:
        alns=keep_alignments,
    output:
        fasta=os.path.join(OUT_DIR, "supermatrix.fasta"),
        parts=os.path.join(OUT_DIR, "supermatrix.partitions"),
    run:
        from collections import OrderedDict

        def read_fasta(path):
            seqs, name, chunks = OrderedDict(), None, []
            for ln in open(path):
                ln = ln.rstrip("\n")
                if ln.startswith(">"):
                    if name is not None:
                        seqs[name] = "".join(chunks)
                    name, chunks = ln[1:].split()[0], []
                else:
                    chunks.append(ln.strip())
            if name is not None:
                seqs[name] = "".join(chunks)
            return seqs

        alns = sorted(input.alns)
        markers = [os.path.basename(a)[:-4] for a in alns]  # strip ".faa"

        # Pass 1: load each alignment, record its width and the global taxon set.
        loaded, widths, all_taxa = [], [], OrderedDict()
        for a in alns:
            seqs = read_fasta(a)
            width = max((len(s) for s in seqs.values()), default=0)
            loaded.append(seqs)
            widths.append(width)
            for t in seqs:
                all_taxa[t] = True
        taxa = list(all_taxa.keys())

        # Pass 2: build each taxon's concatenated row, gap-padding absences.
        os.makedirs(os.path.dirname(output.fasta), exist_ok=True)
        with open(output.fasta, "w") as fh:
            for t in taxa:
                row = []
                for seqs, width in zip(loaded, widths):
                    row.append(seqs.get(t, "-" * width))
                fh.write(">%s\n%s\n" % (t, "".join(row)))

        # Partitions file (protein models, one charset per marker).
        with open(output.parts, "w") as fh:
            start = 1
            for marker, width in zip(markers, widths):
                end = start + width - 1
                if width > 0:
                    fh.write("LG, %s = %d-%d\n" % (marker, start, end))
                    start = end + 1

        print("[supermatrix] %d taxa x %d markers, %d columns"
              % (len(taxa), len(markers), sum(widths)))


# ---- report -----------------------------------------------------------------
rule report:
    # Human-readable summary: keep/drop counts, drop-reason tally, and the
    # clade x marker recovery matrix.  Emits both markdown and a standalone HTML.
    input:
        scorecard=os.path.join(AGG_DIR, "scorecard.tsv"),
        # Pull all monophyly matrices so the clade x marker table is complete.
        mono=expand(os.path.join(MONO_DIR, "{marker}.matrix.tsv"), marker=MARKERS),
    output:
        md=os.path.join(REPORT_DIR, "markervet_report.md"),
        html=os.path.join(REPORT_DIR, "markervet_report.html"),
    params:
        regime=config.get("regime", "?"),
    run:
        import csv as _csv
        from collections import Counter, defaultdict

        # --- scorecard summary ---------------------------------------------
        rows = list(_csv.DictReader(open(input.scorecard), delimiter="\t"))
        n = len(rows)
        keep = [r for r in rows if r["decision"] == "KEEP"]
        drop = [r for r in rows if r["decision"] == "DROP"]

        # Tally individual drop reasons (split the reason string, strip stats).
        reason_counts = Counter()
        for r in drop:
            for reason in r["reasons"].split(";"):
                if reason and reason != "-":
                    reason_counts[reason.split("(")[0]] += 1

        # --- clade x marker recovery matrix --------------------------------
        clade_recovery = defaultdict(dict)   # clade -> {marker: recovery}
        clades = []
        for mf in input.mono:
            for row in _csv.DictReader(open(mf), delimiter="\t"):
                clade_recovery[row["clade"]][row["marker"]] = row["recovery"]
                if row["clade"] not in clades:
                    clades.append(row["clade"])
        markers = [r["marker"] for r in rows]

        # --- markdown -------------------------------------------------------
        md = []
        md.append("# MarkerVet report -- regime: `%s`\n" % params.regime)
        md.append("**%d markers screened** -- **%d KEEP**, **%d DROP** "
                  "(aggressive policy: drop on any failure)\n"
                  % (n, len(keep), len(drop)))
        md.append("## Drop reasons\n")
        if reason_counts:
            md.append("| reason | markers |\n|---|---|")
            for reason, c in reason_counts.most_common():
                md.append("| %s | %d |" % (reason, c))
        else:
            md.append("_No markers dropped._")
        md.append("\n## Decisions\n")
        md.append("| marker | decision | reasons |\n|---|---|---|")
        for r in rows:
            md.append("| %s | %s | %s |" % (r["marker"], r["decision"], r["reasons"]))
        md.append("\n## Clade x marker recovery\n")
        header = "| clade | " + " | ".join(markers) + " |"
        md.append(header)
        md.append("|" + "---|" * (len(markers) + 1))
        for clade in clades:
            cells = [str(clade_recovery[clade].get(m, "NA")) for m in markers]
            md.append("| %s | %s |" % (clade, " | ".join(cells)))
        md_text = "\n".join(md) + "\n"

        os.makedirs(os.path.dirname(output.md), exist_ok=True)
        with open(output.md, "w") as fh:
            fh.write(md_text)

        # --- minimal standalone HTML (no external assets) -------------------
        def esc(s):
            return (str(s).replace("&", "&amp;").replace("<", "&lt;")
                    .replace(">", "&gt;"))

        html = ["<!doctype html><meta charset='utf-8'>",
                "<title>MarkerVet report -- %s</title>" % esc(params.regime),
                "<style>body{font-family:system-ui,sans-serif;margin:2rem;}"
                "table{border-collapse:collapse;margin:1rem 0;}"
                "td,th{border:1px solid #ccc;padding:3px 8px;font-size:13px;}"
                ".KEEP{color:#137333;font-weight:600;}"
                ".DROP{color:#c5221f;font-weight:600;}</style>",
                "<h1>MarkerVet report &mdash; regime: %s</h1>" % esc(params.regime),
                "<p><b>%d</b> markers screened &mdash; <b>%d KEEP</b>, "
                "<b>%d DROP</b> (aggressive policy).</p>" % (n, len(keep), len(drop)),
                "<h2>Drop reasons</h2><table><tr><th>reason</th><th>markers</th></tr>"]
        for reason, c in reason_counts.most_common():
            html.append("<tr><td>%s</td><td>%d</td></tr>" % (esc(reason), c))
        html.append("</table><h2>Decisions</h2>"
                    "<table><tr><th>marker</th><th>decision</th><th>reasons</th></tr>")
        for r in rows:
            html.append("<tr><td>%s</td><td class='%s'>%s</td><td>%s</td></tr>"
                        % (esc(r["marker"]), esc(r["decision"]),
                           esc(r["decision"]), esc(r["reasons"])))
        html.append("</table><h2>Clade &times; marker recovery</h2><table><tr><th>clade</th>")
        html += ["<th>%s</th>" % esc(m) for m in markers]
        html.append("</tr>")
        for clade in clades:
            html.append("<tr><td>%s</td>" % esc(clade))
            html += ["<td>%s</td>" % esc(clade_recovery[clade].get(m, "NA"))
                     for m in markers]
            html.append("</tr>")
        html.append("</table>")
        with open(output.html, "w") as fh:
            fh.write("".join(html))

        print("[report] %d KEEP / %d DROP written to %s"
              % (len(keep), len(drop), output.html))
