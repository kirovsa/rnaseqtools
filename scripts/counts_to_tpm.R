#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  library(optparse)
}))

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Input tab-delimited file with gene, count, and length columns.",
              metavar = "FILE"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output tab-delimited file with TPM column in place of count.",
              metavar = "FILE"),

  make_option(c("--gene-col"), type = "character", default = "gene",
              help = "Column name for gene identifiers. [default %default]"),
  make_option(c("--count-col"), type = "character", default = "count",
              help = "Column name for raw counts. [default %default]"),
  make_option(c("--length-col"), type = "character", default = "length",
              help = "Column name for gene/transcript lengths (in bp). [default %default]"),
  make_option(c("--tpm-col"), type = "character", default = "TPM",
              help = "Column name to use for TPM values in the output. [default %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# ---- Validate required arguments ----
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("--input is required.", call. = FALSE)
}
if (is.null(opt$output)) {
  print_help(opt_parser)
  stop("--output is required.", call. = FALSE)
}

gene_col   <- opt[["gene-col"]]
count_col  <- opt[["count-col"]]
length_col <- opt[["length-col"]]
tpm_col    <- opt[["tpm-col"]]

# ---- Read input ----
dt <- read.table(opt$input, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

for (col in c(gene_col, count_col, length_col)) {
  if (!col %in% names(dt)) {
    stop(sprintf("Column '%s' not found in input file.", col), call. = FALSE)
  }
}

# ---- Compute TPM ----
# Validate lengths
if (any(is.na(dt[[length_col]]) | dt[[length_col]] <= 0, na.rm = FALSE)) {
  stop("All length values must be positive and non-NA.", call. = FALSE)
}

# Step 1: reads per kilobase (RPK)
rpk <- dt[[count_col]] / (dt[[length_col]] / 1000)

# Step 2: per-million scaling factor
scaling_factor <- sum(rpk, na.rm = TRUE) / 1e6
if (!is.finite(scaling_factor) || scaling_factor == 0) {
  stop("Cannot compute TPM: scaling factor is zero or non-finite (all counts may be zero or NA).",
       call. = FALSE)
}

# Step 3: TPM
tpm <- rpk / scaling_factor

# ---- Build output: replace count column with TPM ----
col_order <- names(dt)
col_order[col_order == count_col] <- tpm_col
out <- dt
out[[tpm_col]]   <- tpm
out[[count_col]] <- NULL
out <- out[, col_order, drop = FALSE]

# ---- Write output ----
write.table(out, opt$output, sep = "\t", quote = FALSE, row.names = FALSE)
