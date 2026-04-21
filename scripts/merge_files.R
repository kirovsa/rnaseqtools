#!/usr/bin/env Rscript
#
# merge_files.R
#
# Description:
#   Merge a list of tab-delimited (or CSV) files by a specified key column
#   (identified by 1-based index) and retain only a user-specified set of
#   columns (also identified by 1-based index).  The merge is performed as
#   a sequential full outer join (or inner join) across all input files.
#
# Usage example:
#   Rscript scripts/merge_files.R \
#     --files  sample1.txt,sample2.txt,sample3.txt \
#     --key    1 \
#     --columns 1,2,3 \
#     --output  merged.txt
#
# Author: <author>
#

suppressWarnings(suppressMessages({
  library(optparse)
}))

option_list <- list(
  make_option(c("-f", "--files"), type = "character", default = NULL,
              help = "Comma-separated list of input file paths to merge.",
              metavar = "FILE1,FILE2,..."),
  make_option(c("-k", "--key"), type = "integer", default = 1L,
              help = "1-based column index to use as the merge key (must be present in all files). [default %default]",
              metavar = "INT"),
  make_option(c("-c", "--columns"), type = "character", default = NULL,
              help = paste("Comma-separated list of 1-based column indices to retain from each file",
                           "(the key column is always included). If omitted, all columns are kept."),
              metavar = "INT,INT,..."),
  make_option(c("-o", "--output"), type = "character", default = "merged_output.txt",
              help = "Output file path. [default %default]",
              metavar = "FILE"),
  make_option(c("-s", "--sep"), type = "character", default = "\t",
              help = "Field separator for input and output files. [default: tab]",
              metavar = "SEP"),
  make_option(c("--header"), action = "store_true", default = TRUE,
              help = "Input files have a header row. [default %default]"),
  make_option(c("--join"), type = "character", default = "outer",
              help = "Join type: 'outer' (full outer join, fills missing values with NA) or 'inner' (keep only rows present in all files). [default %default]",
              metavar = "outer|inner"
)
)

opt_parser <- OptionParser(
  option_list = option_list,
  description  = "Merge tab-delimited files by a key column, retaining selected columns."
)
opt <- parse_args(opt_parser)

# ---- Validate required arguments ----
if (is.null(opt$files)) {
  print_help(opt_parser)
  stop("--files is required.", call. = FALSE)
}

join_type <- tolower(opt$join)
if (!join_type %in% c("outer", "inner")) {
  stop("--join must be 'outer' or 'inner'.", call. = FALSE)
}

# ---- Parse file list ----
file_paths <- trimws(strsplit(opt$files, ",")[[1]])
file_paths <- file_paths[nzchar(file_paths)]
if (length(file_paths) < 2L) {
  stop("At least two input files must be supplied via --files.", call. = FALSE)
}

# Check all files exist and are readable
for (fp in file_paths) {
  if (!file.exists(fp)) {
    stop(sprintf("Input file not found: %s", fp), call. = FALSE)
  }
  if (file.access(fp, mode = 4L) != 0L) {
    stop(sprintf("Input file is not readable: %s", fp), call. = FALSE)
  }
}

# ---- Parse column indices ----
key_idx <- opt$key
if (is.na(key_idx) || key_idx < 1L) {
  stop("--key must be a positive integer.", call. = FALSE)
}

keep_idx <- NULL
if (!is.null(opt$columns)) {
  keep_idx <- suppressWarnings(as.integer(trimws(strsplit(opt$columns, ",")[[1]])))
  if (any(is.na(keep_idx)) || any(keep_idx < 1L)) {
    stop("--columns must be a comma-separated list of positive integers.", call. = FALSE)
  }
}

# ---- Read and validate each file ----
tables <- vector("list", length(file_paths))

for (i in seq_along(file_paths)) {
  fp <- file_paths[[i]]
  dt <- tryCatch(
    read.table(fp, sep = opt$sep, header = opt$header, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) stop(sprintf("Error reading file '%s': %s", fp, conditionMessage(e)), call. = FALSE)
  )

  ncols <- ncol(dt)

  # Validate key index
  if (key_idx > ncols) {
    stop(sprintf(
      "File '%s' has only %d column(s); --key %d is out of bounds.",
      fp, ncols, key_idx
    ), call. = FALSE)
  }

  # Determine which column indices to keep for this file
  if (is.null(keep_idx)) {
    idx <- seq_len(ncols)
  } else {
    idx <- keep_idx
  }

  # Always include the key column
  idx <- unique(c(key_idx, idx))

  # Validate all requested indices are within bounds
  oob <- idx[idx > ncols]
  if (length(oob) > 0L) {
    stop(sprintf(
      "File '%s' has only %d column(s); column index/indices out of bounds: %s",
      fp, ncols, paste(oob, collapse = ", ")
    ), call. = FALSE)
  }

  # Subset columns
  dt <- dt[, idx, drop = FALSE]

  # Identify the key column name in the subsetted table
  # (idx[1] is always key_idx because we prepended it above)
  key_pos_in_subset <- which(idx == key_idx)[1L]
  key_colname <- names(dt)[key_pos_in_subset]

  # For non-first files, suffix duplicate non-key column names to avoid collisions
  if (i > 1L) {
    non_key_cols <- setdiff(names(dt), key_colname)
    tables_so_far_names <- unique(unlist(lapply(tables[seq_len(i - 1L)], names)))
    dups <- intersect(non_key_cols, tables_so_far_names)
    for (d in dups) {
      names(dt)[names(dt) == d] <- paste0(d, ".", i)
    }
  }

  # Store the key column name (may have been renamed above, but key column isn't renamed)
  attr(dt, "key_colname") <- key_colname

  tables[[i]] <- dt
}

# ---- Merge all tables ----
# Use the key column name from the first table as the common join key.
# All files must share the same key column name after reading (we keep the original name).
key_name <- attr(tables[[1L]], "key_colname")

# Verify key column is named consistently or rename to a common name
for (i in seq_along(tables)) {
  kn <- attr(tables[[i]], "key_colname")
  if (kn != key_name) {
    names(tables[[i]])[names(tables[[i]]) == kn] <- key_name
  }
}

merged <- tables[[1L]]

for (i in seq(2L, length(tables))) {
  dt_right <- tables[[i]];
  if (join_type == "outer") {
    merged <- merge(merged, dt_right, by = key_name, all = TRUE)
  } else {
    merged <- merge(merged, dt_right, by = key_name, all = FALSE)
  }
}

# ---- Write output ----
write.table(merged, opt$output, sep = opt$sep, quote = FALSE, row.names = FALSE, na = "NA")

# ---- Summary to stderr ----
message(sprintf(
  "merge_files.R: merged %d file(s) -> %d row(s) x %d column(s) written to '%s'",
  length(file_paths), nrow(merged), ncol(merged), opt$output
))
