#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  library(optparse)
  library(data.table)
}))

option_list <- list(
  make_option(c("-m", "--metadata"), type = "character", default = NULL,
              help = "Metadata TSV/CSV with at least: sample, group, rsem_path (path to each sample's RSEM genes.results).",
              metavar = "FILE"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output TSV file of genes that pass group prevalence filter.", metavar = "FILE"),

  make_option(c("--meta-sep"), type = "character", default = "\t",
              help = "Metadata delimiter: \\t for TSV (default) or , for CSV. [default %default]"),
  make_option(c("--sample-col"), type = "character", default = "sample",
              help = "Metadata column name for sample ID. [default %default]"),
  make_option(c("--group-col"), type = "character", default = "group",
              help = "Metadata column name for group label. [default %default]"),
  make_option(c("--path-col"), type = "character", default = "rsem_path",
              help = "Metadata column name for path to RSEM file. [default %default]"),

  make_option(c("--id-col"), type = "character", default = "gene_id",
              help = "Gene ID column name in RSEM files. [default %default]"),
  make_option(c("--count-col"), type = "character", default = "expected_count",
              help = "Count column name in RSEM files. [default %default]"),
  make_option(c("--tpm-col"), type = "character", default = "TPM",
              help = "TPM column name in RSEM files. [default %default]"),

  make_option(c("--min-count"), type = "double", default = 10,
              help = "Minimum count threshold for a gene to be considered 'passing' in a sample. [default %default]"),
  make_option(c("--min-tpm"), type = "double", default = 1,
              help = "Minimum TPM threshold for a gene to be considered 'passing' in a sample. [default %default]"),
  make_option(c("--logic"), type = "character", default = "AND",
              help = "Combine thresholds within a sample: AND or OR. [default %default]"),

  make_option(c("--min-samples"), type = "integer", default = 2,
              help = "A gene is kept if it passes the per-sample filter in at least this many samples in ANY ONE group. [default %default]"),

  make_option(c("--keep-nonfinite"), action = "store_true", default = FALSE,
              help = "If set, do not drop rows with non-finite count/TPM in per-sample evaluation."),
  make_option(c("--write-summary"), action = "store_true", default = FALSE,
              help = "If set, print a summary to stdout.")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# ---- FIX: normalize flag options to strict TRUE/FALSE scalars ----
# This avoids errors like: Error in !opt$keep_nonfinite : invalid argument type
opt$keep_nonfinite <- isTRUE(opt$keep_nonfinite)
opt$write_summary  <- isTRUE(opt$write_summary)
# -----------------------------------------------------------------

if (is.null(opt$metadata) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("Both --metadata and --output are required.", call. = FALSE)
}

logic <- toupper(opt$logic)
if (!logic %in% c("AND", "OR")) stop("--logic must be AND or OR", call. = FALSE)
if (opt$`min-samples` < 1) stop("--min-samples must be >= 1", call. = FALSE)

meta <- fread(opt$metadata, sep = opt$`meta-sep`, header = TRUE, data.table = TRUE)

need_meta <- c(opt$`sample-col`, opt$`group-col`, opt$`path-col`)
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
  group  = as.character(get(opt$`group-col`)),
  path   = as.character(get(opt$`path-col`))
)])

if (anyNA(meta$sample) || any(meta$sample == "")) stop("Metadata has empty/NA sample IDs.", call. = FALSE)
if (anyNA(meta$group)  || any(meta$group == ""))  stop("Metadata has empty/NA group labels.", call. = FALSE)
if (anyNA(meta$path)   || any(meta$path == ""))   stop("Metadata has empty/NA rsem_path values.", call. = FALSE)

# Read each sample and compute pass/fail per gene per sample
per_sample_list <- vector("list", nrow(meta))

for (i in seq_len(nrow(meta))) {
  s <- meta$sample[i]
  g <- meta$group[i]
  p <- meta$path[i]

  dt <- fread(p, sep = "\t", header = TRUE, data.table = TRUE)

  required_cols <- c(opt$`id-col`, opt$`count-col`, opt$`tpm-col`)
  missing <- setdiff(required_cols, names(dt))
  if (length(missing) > 0) {
    stop(sprintf(
      "RSEM file for sample '%s' (%s) is missing required column(s): %s\nAvailable columns: %s",
      s, p, paste(missing, collapse = ", "), paste(names(dt), collapse = ", ")
    ), call. = FALSE)
  }

  dt <- dt[, .(
    gene_id = as.character(get(opt$`id-col`)),
    count   = as.numeric(get(opt$`count-col`)),
    tpm     = as.numeric(get(opt$`tpm-col`))
  )]

  # Use isTRUE() pattern (robust to NULL/non-logical values)
  if (!isTRUE(opt$keep_nonfinite)) {
    dt <- dt[is.finite(count) & is.finite(tpm)]
  }

  if (logic == "AND") {
    dt[, pass := (count >= opt$`min-count`) & (tpm >= opt$`min-tpm`)]
  } else {
    dt[, pass := (count >= opt$`min-count`) | (tpm >= opt$`min-tpm`)]
  }

  per_sample_list[[i]] <- dt[, .(gene_id, pass, sample = s, group = g)]
}

ps <- rbindlist(per_sample_list, use.names = TRUE, fill = TRUE)

# Count passing samples per gene within each group
group_counts <- ps[pass == TRUE, .(n_pass = uniqueN(sample)), by = .(group, gene_id)]

# For each gene, find the max number of passing samples across groups
max_pass_by_gene <- group_counts[, .(max_n_pass_in_any_group = max(n_pass, na.rm = TRUE)), by = gene_id]

# Genes that never pass in any sample won't appear in group_counts; add them with 0
all_genes <- unique(ps$gene_id)
max_pass_by_gene <- merge(
  data.table(gene_id = all_genes),
  max_pass_by_gene,
  by = "gene_id",
  all.x = TRUE
)
max_pass_by_gene[is.na(max_n_pass_in_any_group), max_n_pass_in_any_group := 0L]

kept <- max_pass_by_gene[max_n_pass_in_any_group >= opt$`min-samples`]

# Write kept genes + supporting statistic
setorder(kept, -max_n_pass_in_any_group, gene_id)
fwrite(kept, file = opt$output, sep = "\t", quote = FALSE, na = "NA")

if (isTRUE(opt$write_summary)) {
  cat(sprintf("Samples in metadata: %d\n", nrow(meta)))
  cat(sprintf("Groups: %s\n", paste(sort(unique(meta$group)), collapse = ", ")))
  cat(sprintf("Unique genes observed across all samples: %d\n", length(all_genes)))
  cat(sprintf("Kept genes: %d\n", nrow(kept)))
  cat(sprintf("Filtered genes: %d\n", length(all_genes) - nrow(kept)))
  cat(sprintf("Per-sample thresholds: %s count >= %g, TPM >= %g\n", logic, opt$`min-count`, opt$`min-tpm`))
  cat(sprintf("Group prevalence rule: keep if passes in >= %d samples in ANY ONE group\n", opt$`min-samples`))
}

invisible(TRUE)