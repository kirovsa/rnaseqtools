# rnaseqtools

Tools for preprocessing bulk RNASeq data.
Starting with a filtering script that is flexible and can combine both TPM and counts from RSEM output.

---

## filter_rsem_by_group_prevalence.R

Filter genes from RSEM output based on group-level prevalence: a gene is **kept** if it passes
per-sample count/TPM thresholds in at least a minimum number of samples within **at least one**
biological group.

### Requirements

- R (≥ 4.0)
- R packages: `optparse`, `data.table`

```r
install.packages(c("optparse", "data.table"))
```

### Inputs

| File | Description |
|------|-------------|
| Metadata TSV/CSV | One row per sample; must contain sample ID, group label, and path to RSEM output |
| RSEM `genes.results` files | One per sample; standard RSEM tab-delimited output |

#### Metadata format (default column names)

| sample | group | rsem_path |
|--------|-------|-----------|
| ctrl_1 | control | path/to/ctrl_1_genes.results |
| ctrl_2 | control | path/to/ctrl_2_genes.results |
| treat_1 | treatment | path/to/treat_1_genes.results |
| treat_2 | treatment | path/to/treat_2_genes.results |

#### RSEM `genes.results` format (relevant columns)

| gene_id | … | expected_count | TPM | … |
|---------|---|----------------|-----|---|

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `-m`, `--metadata` | *(required)* | Metadata file (TSV or CSV) |
| `-o`, `--output` | *(required)* | Output TSV of genes passing the filter |
| `--meta-sep` | `\t` | Metadata delimiter (`\t` for TSV, `,` for CSV) |
| `--sample-col` | `sample` | Metadata column for sample ID |
| `--group-col` | `group` | Metadata column for group label |
| `--path-col` | `rsem_path` | Metadata column for path to RSEM file |
| `--id-col` | `gene_id` | Gene ID column name in RSEM files |
| `--count-col` | `expected_count` | Count column name in RSEM files |
| `--tpm-col` | `TPM` | TPM column name in RSEM files |
| `--min-count` | `10` | Minimum count for a sample to be considered "passing" |
| `--min-tpm` | `1` | Minimum TPM for a sample to be considered "passing" |
| `--logic` | `AND` | Combine thresholds within a sample: `AND` or `OR` |
| `--min-samples` | `2` | Minimum number of passing samples required in any one group |
| `--keep-nonfinite` | `FALSE` | Do not drop rows with non-finite count/TPM values |
| `--write-summary` | `FALSE` | Print a brief summary to stdout |

### Filtering logic

For each gene and each sample:

1. A sample is **passing** if:
   - `--logic AND` (default): `expected_count >= --min-count` **and** `TPM >= --min-tpm`
   - `--logic OR`: `expected_count >= --min-count` **or** `TPM >= --min-tpm`

2. Within each group, count the number of passing samples for each gene.

3. A gene is **kept** if, in **at least one** group, it is passing in ≥ `--min-samples` samples.

### Output

A tab-delimited file with two columns:

| gene_id | max_n_pass_in_any_group |
|---------|------------------------|
| GENE_A | 4 |
| GENE_B | 2 |
| … | … |

Rows are ordered by `max_n_pass_in_any_group` (descending), then `gene_id`.

### Example

The `examples/` directory contains ready-to-use input files:

```
examples/
├── metadata.tsv          # 4 samples across 2 groups
├── ctrl_1_genes.results  # RSEM output for control sample 1
├── ctrl_2_genes.results  # RSEM output for control sample 2
├── treat_1_genes.results # RSEM output for treatment sample 1
└── treat_2_genes.results # RSEM output for treatment sample 2
```

**Run with default settings** (AND logic, min-count 10, min-TPM 1, min-samples 2):

```bash
Rscript scripts/filter_rsem_by_group_prevalence.R \
  --metadata examples/metadata.tsv \
  --output   filtered_genes.tsv \
  --write-summary
```

Expected output (`filtered_genes.tsv`):

```
gene_id	max_n_pass_in_any_group
GENE_A	2
GENE_B	2
GENE_E	2
```

- `GENE_A` passes in both control samples and both treatment samples (2 passing in each group; `max_n_pass_in_any_group` = 2) → kept.
- `GENE_B` passes in both control samples (2/2, count ≥ 10 **and** TPM ≥ 1) but 0 treatment samples (count < 10 and TPM < 1) → kept because the control group has ≥ 2 passing samples.
- `GENE_E` passes in both treatment samples, 0 control samples → kept.
- `GENE_C` passes in only 1 sample per group → **filtered**.
- `GENE_D` never reaches the thresholds → **filtered**.

**Use OR logic** (keep a gene if it passes count **or** TPM threshold in a sample):

```bash
Rscript scripts/filter_rsem_by_group_prevalence.R \
  --metadata  examples/metadata.tsv \
  --output    filtered_genes_or.tsv \
  --logic     OR \
  --min-count 10 \
  --min-tpm   1
```

**Require passing in at least 3 samples within a group:**

```bash
Rscript scripts/filter_rsem_by_group_prevalence.R \
  --metadata    examples/metadata.tsv \
  --output      filtered_genes_strict.tsv \
  --min-samples 3
```

**Use a CSV metadata file with custom column names:**

```bash
Rscript scripts/filter_rsem_by_group_prevalence.R \
  --metadata   my_samples.csv \
  --meta-sep   , \
  --sample-col SampleID \
  --group-col  Condition \
  --path-col   FilePath \
  --output     filtered_genes.tsv
```
