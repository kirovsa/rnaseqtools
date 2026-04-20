#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  library(optparse)
  library(data.table)
}))

option_list <- list(
  make_option(c("--counts"), type = "character", default = NULL,
              help = "Multi-sample counts file (wide format: genes as rows, samples as columns).",
              metavar = "FILE"),
  make_option(c("--fpkm"), type = "character", default = NULL,
              help = "Multi-sample FPKM file (wide format: genes as rows, samples as columns).",
              metavar = "FILE"),
  make_option(c("-m", "--metadata"), type = "character", default = NULL,
              help = paste0("Metadata TSV/CSV with at least: sample, group columns. ",
                            "If omitted, all samples are treated as one group."),
              metavar = "FILE"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output TSV file of genes that pass the group prevalence filter.", metavar = "FILE"),

  make_option(c("--meta-sep"), type = "character", default = "\t",
              help = "Metadata delimiter: \\t for TSV (default) or , for CSV. [default %default]"),
  make_option(c("--sample-col"), type = "character", default = "sample",
              help = "Metadata column name for sample ID. [default %default]"),
  make_option(c("--group-col"), type = "character", default = "group",
              help = "Metadata column name for group label. [default %default]"),

  make_option(c("--counts-gene-col"), type = "character", default = NULL,
              help = paste0("Column name in the counts file used as the gene identifier. ",
                            "Defaults to the first column.")),
  make_option(c("--counts-skip-cols"), type = "character", default = "geneID",
              help = paste0("Comma-separated list of extra (non-sample) columns in the counts file ",
                            "to ignore (e.g. gene-symbol column). [default %default]")),
  make_option(c("--fpkm-gene-col"), type = "character", default = "Ensembl_ID",
              help = "Column name in the FPKM file used as the gene identifier. [default %default]"),
  make_option(c("--fpkm-skip-cols"), type = "character", default = "GeneSymbol",
              help = paste0("Comma-separated list of extra (non-sample) columns in the FPKM file ",
                            "to ignore (e.g. gene-symbol column). [default %default]")),

  make_option(c("--min-count"), type = "double", default = 10,
              help = "Minimum count for a gene to be considered 'passing' in a sample. [default %default]"),
  make_option(c("--min-fpkm"), type = "double", default = 1,
              help = "Minimum FPKM for a gene to be considered 'passing' in a sample. [default %default]"),
  make_option(c("--logic"), type = "character", default = "AND",
              help = "Combine count and FPKM thresholds within a sample: AND or OR. [default %default]"),

  make_option(c("--min-samples"), type = "integer", default = 2,
              help = paste0("Keep a gene if it passes the per-sample filter in at least this many ",
                            "samples in ANY ONE group. [default %default]")),

  make_option(c("--keep-nonfinite"), action = "store_true", default = FALSE,
              help = "If set, do not drop genes with non-finite count/FPKM during per-sample evaluation."),
  make_option(c("--write-summary"), action = "store_true", default = FALSE,
              help = "If set, print a summary table to stdout.")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Normalise flag options to strict TRUE/FALSE scalars
opt$keep_nonfinite <- isTRUE(opt$keep_nonfinite)
opt$write_summary  <- isTRUE(opt$write_summary)

if (is.null(opt$counts) || is.null(opt$fpkm) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("--counts, --fpkm, and --output are required.", call. = FALSE)
}

logic <- toupper(opt$logic)
if (!logic %in% c("AND", "OR")) stop("--logic must be AND or OR", call. = FALSE)
if (opt$`min-samples` < 1)      stop("--min-samples must be >= 1", call. = FALSE)

# ---- Helper: parse a comma-separated list, trimming whitespace ----
parse_list <- function(x) {
  if (is.null(x) || !nzchar(trimws(x))) return(character(0))
  trimws(strsplit(x, ",")[[1]])
}

counts_skip <- parse_list(opt$`counts-skip-cols`)
fpkm_skip   <- parse_list(opt$`fpkm-skip-cols`)

# ---- Read counts ----
counts_dt <- fread(opt$counts, sep = "\t", header = TRUE, data.table = TRUE)

# Rename any blank column names (e.g. when the file has a leading delimiter in the header row)
blank_cols <- which(!nzchar(names(counts_dt)))
if (length(blank_cols) > 0) {
  names(counts_dt)[blank_cols] <- paste0("V", blank_cols)
}

counts_gene_col <- opt$`counts-gene-col`
if (is.null(counts_gene_col) || !nzchar(counts_gene_col)) {
  counts_gene_col <- names(counts_dt)[1]
}
if (!counts_gene_col %in% names(counts_dt)) {
  stop(sprintf("counts-gene-col '%s' not found in counts file.\nAvailable columns: %s",
               counts_gene_col, paste(names(counts_dt), collapse = ", ")), call. = FALSE)
}

counts_non_sample <- unique(c(counts_gene_col, counts_skip))
counts_samples    <- setdiff(names(counts_dt), counts_non_sample)
if (length(counts_samples) == 0) {
  stop("No sample columns found in the counts file after removing gene-ID and skip columns.", call. = FALSE)
}

# ---- Read FPKM ----
fpkm_dt <- fread(opt$fpkm, sep = "\t", header = TRUE, data.table = TRUE)

fpkm_gene_col <- opt$`fpkm-gene-col`
if (!fpkm_gene_col %in% names(fpkm_dt)) {
  stop(sprintf("fpkm-gene-col '%s' not found in FPKM file.\nAvailable columns: %s",
               fpkm_gene_col, paste(names(fpkm_dt), collapse = ", ")), call. = FALSE)
}

fpkm_non_sample <- unique(c(fpkm_gene_col, fpkm_skip))
fpkm_samples    <- setdiff(names(fpkm_dt), fpkm_non_sample)
if (length(fpkm_samples) == 0) {
  stop("No sample columns found in the FPKM file after removing gene-ID and skip columns.", call. = FALSE)
}

# ---- Resolve shared samples ----
common_samples <- intersect(counts_samples, fpkm_samples)
if (length(common_samples) == 0) {
  stop(sprintf(
    "No common sample columns found between counts and FPKM files.\nCounts samples: %s\nFPKM samples: %s",
    paste(counts_samples, collapse = ", "),
    paste(fpkm_samples,   collapse = ", ")
  ), call. = FALSE)
}

if (length(common_samples) < length(counts_samples) || length(common_samples) < length(fpkm_samples)) {
  warning(sprintf(
    "Some samples are present in only one file and will be ignored.\nUsing %d common samples: %s",
    length(common_samples), paste(common_samples, collapse = ", ")
  ))
}

# ---- Load or construct sample-to-group mapping ----
if (!is.null(opt$metadata)) {
  meta <- fread(opt$metadata, sep = opt$`meta-sep`, header = TRUE, data.table = TRUE)

  need_meta <- c(opt$`sample-col`, opt$`group-col`)
  missing_meta <- setdiff(need_meta, names(meta))
  if (length(missing_meta) > 0) {
    stop(sprintf(
      "Metadata is missing required column(s): %s\nAvailable columns: %s",
      paste(missing_meta, collapse = ", "),
      paste(names(meta), collapse = ", ")
    ), call. = FALSE)
  }

  meta <- unique(meta[, .(
    sample = as.character(get(opt$`sample-col`)),
    group  = as.character(get(opt$`group-col`))
  )])

  if (anyNA(meta$sample) || any(meta$sample == "")) stop("Metadata has empty/NA sample IDs.", call. = FALSE)
  if (anyNA(meta$group)  || any(meta$group == ""))  stop("Metadata has empty/NA group labels.", call. = FALSE)

  # Restrict to samples present in both files
  meta <- meta[sample %in% common_samples]
  if (nrow(meta) == 0) {
    stop("No metadata entries match the common samples found in counts and FPKM files.", call. = FALSE)
  }

  analysis_samples <- meta$sample
} else {
  meta <- data.table(sample = common_samples, group = "all")
  analysis_samples <- common_samples
}

# ---- Build per-sample pass/fail table ----
# Extract relevant columns from counts and fpkm
counts_sub <- counts_dt[, c(counts_gene_col, analysis_samples), with = FALSE]
setnames(counts_sub, counts_gene_col, "gene_id")

fpkm_sub <- fpkm_dt[, c(fpkm_gene_col, analysis_samples), with = FALSE]
setnames(fpkm_sub, fpkm_gene_col, "gene_id")

counts_sub[, gene_id := as.character(gene_id)]
fpkm_sub[,  gene_id := as.character(gene_id)]

# Keep only genes present in both files
all_genes <- intersect(counts_sub$gene_id, fpkm_sub$gene_id)
if (length(all_genes) == 0) {
  stop("No shared gene IDs found between counts and FPKM files.", call. = FALSE)
}
counts_sub <- counts_sub[gene_id %in% all_genes]
fpkm_sub   <- fpkm_sub[gene_id %in% all_genes]

# Melt to long format for vectorised evaluation
counts_long <- melt(counts_sub, id.vars = "gene_id",
                    variable.name = "sample", value.name = "count")
fpkm_long   <- melt(fpkm_sub,   id.vars = "gene_id",
                    variable.name = "sample", value.name = "fpkm")

counts_long[, count := as.numeric(count)]
fpkm_long[,  fpkm  := as.numeric(fpkm)]

# Merge counts and FPKM long tables
ps <- merge(counts_long, fpkm_long, by = c("gene_id", "sample"), all = FALSE)

# Attach group labels
ps <- merge(ps, meta, by = "sample", all.x = TRUE)
ps <- ps[!is.na(group)]   # drop samples without a group assignment

if (!opt$keep_nonfinite) {
  ps <- ps[is.finite(count) & is.finite(fpkm)]
}

# Evaluate per-sample pass/fail
if (logic == "AND") {
  ps[, pass := (count >= opt$`min-count`) & (fpkm >= opt$`min-fpkm`)]
} else {
  ps[, pass := (count >= opt$`min-count`) | (fpkm >= opt$`min-fpkm`)]
}

# Count passing samples per gene per group
group_counts <- ps[pass == TRUE, .(n_pass = uniqueN(sample)), by = .(group, gene_id)]

# For each gene, find the max passing-sample count across all groups
max_pass_by_gene <- group_counts[, .(max_n_pass_in_any_group = max(n_pass, na.rm = TRUE)), by = gene_id]

# Genes that never pass in any sample won't appear in group_counts; fill with 0
max_pass_by_gene <- merge(
  data.table(gene_id = all_genes),
  max_pass_by_gene,
  by = "gene_id",
  all.x = TRUE
)
max_pass_by_gene[is.na(max_n_pass_in_any_group), max_n_pass_in_any_group := 0L]

kept <- max_pass_by_gene[max_n_pass_in_any_group >= opt$`min-samples`]

setorder(kept, -max_n_pass_in_any_group, gene_id)
fwrite(kept, file = opt$output, sep = "\t", quote = FALSE, na = "NA")

if (opt$write_summary) {
  cat(sprintf("Samples analysed: %d\n", length(analysis_samples)))
  cat(sprintf("Groups: %s\n", paste(sort(unique(meta$group)), collapse = ", ")))
  cat(sprintf("Unique genes in both files: %d\n", length(all_genes)))
  cat(sprintf("Kept genes: %d\n", nrow(kept)))
  cat(sprintf("Filtered genes: %d\n", length(all_genes) - nrow(kept)))
  cat(sprintf("Per-sample thresholds: %s  count >= %g, FPKM >= %g\n",
              logic, opt$`min-count`, opt$`min-fpkm`))
  cat(sprintf("Group prevalence rule: keep if passes in >= %d samples in ANY ONE group\n",
              opt$`min-samples`))
}

invisible(TRUE)
