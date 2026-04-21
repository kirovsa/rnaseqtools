# rnaseqtools

Tools for preprocessing bulk RNASeq data.
Starting with a filtering script that is flexible and can combine both TPM and counts from RSEM output.

---

## counts_to_tpm.R

Convert raw counts to TPM (Transcripts Per Million) given per-gene lengths.
Reads a single tab-delimited file with `gene`, `count`, and `length` columns and writes an
equivalent file with the `count` column replaced by `TPM`.

### Requirements

- R (≥ 4.0)
- R packages: `optparse`, `data.table`

```r
install.packages(c("optparse", "data.table"))
```

### Inputs

| File | Description |
|------|-------------|
| Input TSV | One row per gene; must contain a gene ID column, a raw-count column, and a gene-length column (in bp) |

#### Input format (default column names)

| gene | count | length |
|------|-------|--------|
| GENE_A | 250 | 1500 |
| GENE_B | 180 | 2000 |

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `-i`, `--input` | *(required)* | Input tab-delimited file |
| `-o`, `--output` | *(required)* | Output tab-delimited file |
| `--gene-col` | `gene` | Column name for gene identifiers |
| `--count-col` | `count` | Column name for raw counts |
| `--length-col` | `length` | Column name for gene/transcript lengths (bp) |
| `--tpm-col` | `TPM` | Column name to use for TPM values in the output |

### TPM calculation

For each gene *i*:

1. **RPK** = count_i / (length_i / 1,000)
2. **Scaling factor** = Σ RPK / 1,000,000
3. **TPM_i** = RPK_i / scaling factor

By construction, TPM values sum to 1,000,000 across all genes.

### Output

A tab-delimited file with the same columns as the input except the count column is replaced by TPM:

| gene | TPM | length |
|------|-----|--------|
| GENE_A | 481347.77 | 1500 |
| … | … | … |

### Example

```bash
Rscript scripts/counts_to_tpm.R \
  --input  examples/example_counts_length.tsv \
  --output tpm_output.tsv
```

**Custom column names:**

```bash
Rscript scripts/counts_to_tpm.R \
  --input      my_data.tsv \
  --output     my_data_tpm.tsv \
  --gene-col   gene_id \
  --count-col  expected_count \
  --length-col transcript_length \
  --tpm-col    tpm
```

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

---

## filter_multi_sample_by_group_prevalence.R

Filter genes from pre-merged wide-format count and FPKM matrices based on group-level prevalence:
a gene is **kept** if it passes per-sample count/FPKM thresholds in at least a minimum number of
samples within **at least one** biological group.

Use this script when you already have multi-sample count and FPKM matrices (genes as rows,
samples as columns) rather than individual per-sample RSEM files.

### Requirements

- R (≥ 4.0)
- R packages: `optparse`, `data.table`

```r
install.packages(c("optparse", "data.table"))
```

### Inputs

| File | Description |
|------|-------------|
| Counts matrix | Wide-format TSV; genes as rows, samples as columns |
| FPKM matrix | Wide-format TSV; genes as rows, samples as columns |
| Metadata TSV/CSV *(optional)* | One row per sample; must contain a sample ID column and a group label column. If omitted, all samples are treated as a single group. |

#### Counts matrix format

The first column is used as the gene identifier by default. Any additional non-sample columns
(e.g. a gene-symbol column) are listed in `--counts-skip-cols`.

| *(gene_id)* | geneID | sample_1 | sample_2 | … |
|-------------|--------|----------|----------|---|
| ENSG00000160072 | ATAD3B | 2230 | 2035 | … |

#### FPKM matrix format (default column names)

| GeneSymbol | Ensembl_ID | sample_1 | sample_2 | … |
|------------|------------|----------|----------|---|
| ATAD3B | ENSG00000160072 | 17.42 | 16.84 | … |

#### Metadata format (default column names)

| sample | group |
|--------|-------|
| sample_1 | control |
| sample_2 | control |
| sample_3 | treatment |
| sample_4 | treatment |

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--counts` | *(required)* | Multi-sample counts file (wide format) |
| `--fpkm` | *(required)* | Multi-sample FPKM file (wide format) |
| `-o`, `--output` | *(required)* | Output TSV of genes passing the filter |
| `-m`, `--metadata` | *(optional)* | Metadata file (TSV or CSV); if omitted all samples form one group |
| `--meta-sep` | `\t` | Metadata delimiter (`\t` for TSV, `,` for CSV) |
| `--sample-col` | `sample` | Metadata column for sample ID |
| `--group-col` | `group` | Metadata column for group label |
| `--counts-gene-col` | *(first column)* | Column name in the counts file used as gene identifier |
| `--counts-skip-cols` | `geneID` | Comma-separated non-sample columns in the counts file to ignore |
| `--fpkm-gene-col` | `Ensembl_ID` | Column name in the FPKM file used as gene identifier |
| `--fpkm-skip-cols` | `GeneSymbol` | Comma-separated non-sample columns in the FPKM file to ignore |
| `--min-count` | `10` | Minimum count for a sample to be considered "passing" |
| `--min-fpkm` | `1` | Minimum FPKM for a sample to be considered "passing" |
| `--logic` | `AND` | Combine thresholds within a sample: `AND` or `OR` |
| `--min-samples` | `2` | Minimum number of passing samples required in any one group |
| `--keep-nonfinite` | `FALSE` | Do not drop rows with non-finite count/FPKM values |
| `--write-summary` | `FALSE` | Print a brief summary to stdout |

### Filtering logic

For each gene and each sample:

1. A sample is **passing** if:
   - `--logic AND` (default): `count >= --min-count` **and** `FPKM >= --min-fpkm`
   - `--logic OR`: `count >= --min-count` **or** `FPKM >= --min-fpkm`

2. Within each group, count the number of passing samples for each gene.

3. A gene is **kept** if, in **at least one** group, it is passing in ≥ `--min-samples` samples.

Only genes whose identifier appears in **both** the counts and FPKM files are evaluated. Only
samples whose name appears in **both** files are used.

### Output

A tab-delimited file with two columns:

| gene_id | max_n_pass_in_any_group |
|---------|------------------------|
| ENSG00000160072 | 5 |
| … | … |

Rows are ordered by `max_n_pass_in_any_group` (descending), then `gene_id`.

### Example

The `examples/` directory contains ready-to-use wide-format input files:

```
examples/
├── multi_sample_metadata.tsv         # 20 samples across 4 groups (PINT, PINU, SCRT, SCRU)
├── multi_sample_example_counts.txt   # counts matrix (genes × samples)
└── multi_sample_example_fpkm.txt    # FPKM matrix (genes × samples)
```

**Run with default settings** (AND logic, min-count 10, min-FPKM 1, min-samples 2):

```bash
Rscript scripts/filter_multi_sample_by_group_prevalence.R \
  --counts   examples/multi_sample_example_counts.txt \
  --fpkm     examples/multi_sample_example_fpkm.txt \
  --metadata examples/multi_sample_metadata.tsv \
  --output   filtered_genes.tsv \
  --write-summary
```

Expected output (`filtered_genes.tsv`):

```
gene_id	max_n_pass_in_any_group
ENSG00000160072	5
```

`ENSG00000160072` passes count ≥ 10 **and** FPKM ≥ 1 in all 5 samples of every group → kept
with `max_n_pass_in_any_group` = 5.

**Use OR logic** (keep a gene if it passes count **or** FPKM threshold in a sample):

```bash
Rscript scripts/filter_multi_sample_by_group_prevalence.R \
  --counts   examples/multi_sample_example_counts.txt \
  --fpkm     examples/multi_sample_example_fpkm.txt \
  --metadata examples/multi_sample_metadata.tsv \
  --output   filtered_genes_or.tsv \
  --logic    OR \
  --min-count 10 \
  --min-fpkm  1
```

**Require passing in at least 4 samples within a group:**

```bash
Rscript scripts/filter_multi_sample_by_group_prevalence.R \
  --counts      examples/multi_sample_example_counts.txt \
  --fpkm        examples/multi_sample_example_fpkm.txt \
  --metadata    examples/multi_sample_metadata.tsv \
  --output      filtered_genes_strict.tsv \
  --min-samples 4
```

**Run without a metadata file** (all samples treated as one group):

```bash
Rscript scripts/filter_multi_sample_by_group_prevalence.R \
  --counts  examples/multi_sample_example_counts.txt \
  --fpkm    examples/multi_sample_example_fpkm.txt \
  --output  filtered_genes_nogroup.tsv
```

**Use custom gene-ID and skip-column names:**

```bash
Rscript scripts/filter_multi_sample_by_group_prevalence.R \
  --counts           my_counts.txt \
  --fpkm             my_fpkm.txt \
  --metadata         my_samples.tsv \
  --counts-gene-col  gene_id \
  --counts-skip-cols symbol,biotype \
  --fpkm-gene-col    ensembl_id \
  --fpkm-skip-cols   gene_name \
  --output           filtered_genes.tsv
```
