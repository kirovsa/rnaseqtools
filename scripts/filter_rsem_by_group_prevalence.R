# Updated R Script: filter_rsem_by_group_prevalence.R

# Setting optparse flags
opt$keep_nonfinite <- isTRUE(opt$keep_nonfinite)
op$write_summary <- isTRUE(opt$write_summary)

# Assuming the rest of your existing logic here...

# Replace conditions accordingly
if (!isTRUE(opt$keep_nonfinite)) {
    # Your logic here
}

# Guard summary output
if (isTRUE(opt$write_summary)) {
    # Summary logic here
}