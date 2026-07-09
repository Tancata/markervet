# =============================================================================
# MarkerVet  --  config-driven screening of phylogenetic marker genes
# =============================================================================
#
# Purpose
# -------
# Screen candidate single-copy marker genes for *congruence* and *horizontal
# gene transfer (HGT)* so that only phylogenetically well-behaved markers are
# passed on to concatenation (supermatrix) and coalescence (gene-tree) analyses
# of deep relationships:
#
#   (a) prokaryote domain-level trees (Archaea + Bacteria), and
#   (b) the placement of eukaryotes relative to Archaea and Bacteria.
#
# Given a set of candidate markers sampled across a fixed taxon set, MarkerVet
# returns, per marker, a keep/drop decision plus cleaned and trimmed alignments
# and collapsed gene trees ready for downstream analysis.
#
# Design principles (see README.md for the full rationale)
# --------------------------------------------------------
#   * Collapse low-support branches BEFORE any tree comparison -- single-gene
#     trees at these depths are dominated by estimation error, so uncollapsed
#     comparisons measure noise, not real conflict.
#   * Run the artefact screen (composition, saturation, rogue tips) BEFORE HGT
#     calling, so that compositional attraction and saturation -- the classic
#     eukaryote-placement failure modes -- are not mislabelled as transfer.
#   * Take two independent HGT / outlier votes: PhylteR (broad, distance-based
#     outlier detection) and reconciliation (transfer localisation).
#   * Treat monophyly recovery as both a filter and a headline scientific
#     result.
#
# Drop policy is AGGRESSIVE by default (conservative-keep / aggressive-drop):
# a marker is dropped if it fails ANY enabled criterion.  Every threshold is
# exposed in the regime YAML under config/.
#
# Usage
# -----
#   snakemake --configfile config/prokaryote_deep.yaml --cores 8 --use-conda
#   snakemake --configfile config/three_domain.yaml   --cores 8 --use-conda
#
# The pipeline is a correctly-wired scaffold: the pure-Python / R steps contain
# real logic, and the external-tool rules (MAFFT, ClipKIT, IQ-TREE, ...) carry
# runnable shell templates.  It is not guaranteed to run end-to-end until those
# external tools and the marker/reference resources are installed.
# =============================================================================

import os
import glob
import sys

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
# A regime preset must be supplied on the command line, e.g.
#   --configfile config/prokaryote_deep.yaml
# If the user forgets, fall back to the prokaryote_deep preset so that a bare
# `snakemake -n` still parses.
if not config:
    configfile: "config/prokaryote_deep.yaml"

# Root output directory for this run (kept separate per regime so the two
# presets never clobber each other's results).
RESULTS = config.get("workdir", "results")

# Convenience accessors used across the included rule files.
SCRIPTS = os.path.join(workflow.basedir, "scripts")


# -----------------------------------------------------------------------------
# Marker discovery
# -----------------------------------------------------------------------------
# Two mutually exclusive input modes, selected by config["markers"]["extract"]:
#
#   extract: false  -- the user supplies one unaligned per-marker FASTA per
#                      marker in markers.fasta_dir  ({marker}{fasta_ext}).
#                      The extraction stage is skipped entirely.
#
#   extract: true   -- the user supplies proteomes plus an HMM per marker in
#                      markers.hmm_dir ({marker}.hmm).  The extraction stage
#                      runs hmmsearch + single-copy enforcement to build the
#                      per-marker FASTAs.
#
# In both modes MARKERS is the authoritative list of marker names that the rest
# of the pipeline fans out over.
def discover_markers():
    m = config["markers"]
    if m.get("hmm_list"):
        # Explicit list wins if provided (a plain text file, one name per line,
        # or an inline YAML list).
        if isinstance(m["hmm_list"], list):
            return list(m["hmm_list"])
        with open(m["hmm_list"]) as fh:
            return [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
    if m.get("extract", False):
        # Discover from the HMM directory.
        hmms = sorted(glob.glob(os.path.join(m["hmm_dir"], "*.hmm")))
        return [os.path.splitext(os.path.basename(h))[0] for h in hmms]
    # Discover from supplied per-marker FASTAs.
    ext = m.get("fasta_ext", ".faa")
    fastas = sorted(glob.glob(os.path.join(m["fasta_dir"], "*" + ext)))
    return [os.path.basename(f)[: -len(ext)] for f in fastas]


MARKERS = discover_markers()

# A degenerate marker set (e.g. resources not yet staged) should not crash a
# dry-run; warn instead so the DAG still builds.
if not MARKERS:
    sys.stderr.write("[MarkerVet] WARNING: no markers discovered -- check the "
                     "markers.* paths in your configfile.\n")


# -----------------------------------------------------------------------------
# Wildcard constraints
# -----------------------------------------------------------------------------
# Marker names may contain dots/underscores but never a path separator; pin the
# wildcard so Snakemake's rule matching stays unambiguous.
wildcard_constraints:
    marker=r"[^/]+",


# -----------------------------------------------------------------------------
# Stage includes -- one rules/*.smk per pipeline stage
# -----------------------------------------------------------------------------
include: "rules/extraction.smk"      # 1. hmmsearch + single-copy enforcement
include: "rules/align_trim.smk"      # 2. MAFFT then ClipKIT / BMGE
include: "rules/artefact_screen.smk" # 3. composition, saturation, TreeShrink
include: "rules/gene_trees.smk"      # 4. IQ-TREE + support-based collapse
include: "rules/congruence.smk"      # 5. gCF/sCF + RF + quartet vs reference
include: "rules/monophyly.smk"       # 6. MonoPhy-style clade recovery
include: "rules/outliers_hgt.smk"    # 7. PhylteR + reconciliation
include: "rules/aggregate.smk"       # 8. aggressive scorecard
include: "rules/outputs.smk"         # 9. cleaned outputs + supermatrix + report


# -----------------------------------------------------------------------------
# Default target
# -----------------------------------------------------------------------------
# `rule all` requests the end-of-pipeline deliverables.  outputs.smk expresses
# the passer-only alignments / trees and the concatenated supermatrix in terms
# of the scorecard, so requesting the report + supermatrix pulls the whole DAG.
rule all:
    input:
        scorecard=os.path.join(RESULTS, "aggregate", "scorecard.tsv"),
        report=os.path.join(RESULTS, "report", "markervet_report.html"),
        supermatrix=os.path.join(RESULTS, "outputs", "supermatrix.fasta"),
        partitions=os.path.join(RESULTS, "outputs", "supermatrix.partitions"),
        genetrees=os.path.join(RESULTS, "outputs", "gene_trees.done"),
