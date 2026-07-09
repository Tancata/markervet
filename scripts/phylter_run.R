#!/usr/bin/env Rscript
# =============================================================================
# phylter_run.R  --  PhylteR outlier detection across markers (MarkerVet stage 7)
# -----------------------------------------------------------------------------
# PhylteR (Comte et al. 2023) is the modern successor to Phylo-MCOA: it places
# every gene's patristic-distance matrix into a common space with DISTATIS and
# flags gene x taxon *cells* that are outliers -- i.e. taxa whose position in a
# given gene disagrees with their consensus position across all genes.  This is
# MarkerVet's broad, distance-based HGT/outlier vote, complementary to the
# reconciliation vote which localises individual transfers.
#
# It runs ONCE over the whole collapsed-gene-tree set (PhylteR needs the full
# gene x taxon block to define the consensus), then writes:
#
#   * <out_cells>   : long table of removed (gene, taxon) outlier cells;
#   * <out_summary> : one row per marker -- cells removed, fraction of that
#                     marker's taxa removed, and a gene-level-outlier flag --
#                     which the aggregator turns into keep/drop calls.
#
# A marker is a candidate drop if PhylteR flags it as a gene-level outlier OR it
# loses more than the configured fraction of its taxon cells (both thresholds
# are applied downstream in aggregate_scorecard.py; here we just report).
# =============================================================================

suppressWarnings(suppressMessages({
  library(optparse)
  library(ape)
}))

option_list <- list(
  make_option("--tree-dir", type = "character",
              help = "directory of per-marker collapsed gene trees"),
  make_option("--tree-ext", type = "character", default = ".nwk",
              help = "gene-tree filename extension [default %default]"),
  make_option("--out-cells", type = "character",
              help = "output TSV: removed (marker, taxon) outlier cells"),
  make_option("--out-summary", type = "character",
              help = "output TSV: per-marker outlier summary"),
  make_option("--gene-outlier-frac", type = "double", default = 0.5,
              help = paste("mark a gene a gene-level outlier when this fraction",
                           "of its taxa are removed [default %default]"))
)
opt <- parse_args(OptionParser(option_list = option_list))

# --- gather the gene trees ---------------------------------------------------
files <- list.files(opt$`tree-dir`, pattern = paste0("\\",
                    opt$`tree-ext`, "$"), full.names = TRUE)
if (length(files) == 0) {
  stop(sprintf("phylter_run.R: no '*%s' trees in %s", opt$`tree-ext`,
               opt$`tree-dir`))
}
marker_names <- sub(paste0("\\", opt$`tree-ext`, "$"), "", basename(files))

# Read each Newick tree; PhylteR accepts a list of ape::phylo objects.
trees <- lapply(files, function(f) {
  t <- tryCatch(read.tree(f), error = function(e) NULL)
  t
})
ok <- !vapply(trees, is.null, logical(1))
trees <- trees[ok]
marker_names <- marker_names[ok]

# --- run PhylteR -------------------------------------------------------------
# phylter() builds per-gene distance matrices, aligns them with DISTATIS, and
# iteratively removes outlier cells.  We keep its defaults except for handing in
# our marker names so the outputs are labelled.
suppressWarnings(suppressMessages(library(phylter)))
results <- phylter(trees, gene.names = marker_names)

# $Final$Outliers is a two-column matrix/data.frame: (gene, species) cells that
# PhylteR removed.  Guard against the no-outliers case (NULL / zero rows).
outliers <- results$Final$Outliers
if (is.null(outliers) || nrow(outliers) == 0) {
  outliers <- data.frame(gene = character(0), taxon = character(0),
                         stringsAsFactors = FALSE)
} else {
  outliers <- data.frame(gene = as.character(outliers[, 1]),
                         taxon = as.character(outliers[, 2]),
                         stringsAsFactors = FALSE)
}

# --- write the long cell table ----------------------------------------------
dir.create(dirname(opt$`out-cells`), showWarnings = FALSE, recursive = TRUE)
write.table(outliers, file = opt$`out-cells`, sep = "\t", quote = FALSE,
            row.names = FALSE)

# --- per-marker summary ------------------------------------------------------
# Number of taxa originally in each gene (denominator for the removed fraction).
taxa_count <- vapply(trees, function(t) length(t$tip.label), integer(1))
names(taxa_count) <- marker_names

removed_per_gene <- table(factor(outliers$gene, levels = marker_names))

# PhylteR also reports whole genes it considers globally discordant; if the
# installed version exposes them under $Final$CompleteOutliers use that,
# otherwise fall back to the removed-fraction heuristic below.
complete_genes <- tryCatch(as.character(results$Final$CompleteOutliers),
                           error = function(e) character(0))

summary_df <- data.frame(
  marker = marker_names,
  n_taxa = as.integer(taxa_count[marker_names]),
  cells_removed = as.integer(removed_per_gene[marker_names]),
  stringsAsFactors = FALSE
)
summary_df$cells_removed[is.na(summary_df$cells_removed)] <- 0L
summary_df$frac_removed <- ifelse(summary_df$n_taxa > 0,
                                  summary_df$cells_removed / summary_df$n_taxa, 0)
summary_df$gene_outlier <- (summary_df$marker %in% complete_genes) |
                           (summary_df$frac_removed >= opt$`gene-outlier-frac`)

dir.create(dirname(opt$`out-summary`), showWarnings = FALSE, recursive = TRUE)
write.table(summary_df, file = opt$`out-summary`, sep = "\t", quote = FALSE,
            row.names = FALSE)

cat(sprintf("[phylter_run] %d markers, %d outlier cells removed, %d gene-level outliers\n",
            length(marker_names), nrow(outliers), sum(summary_df$gene_outlier)),
    file = stderr())
