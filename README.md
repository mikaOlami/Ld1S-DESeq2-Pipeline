# Ld1S DESeq2 Pipeline

DESeq2 analysis workflow for *Leishmania donovani* (Ld1S), starting from raw paired-end FASTQ files and producing:
- genome-aligned, sorted, indexed BAM files
- gene-level count matrix
- sample metadata (colData)
- PCA plot (`PCA_plot.png`)
- differential expression results CSV

The workflow uses SMALT for genome mapping, `bedtools multicov` for gene-level counting, and DESeq2 for differential expression.

---

## Requirements

### Command-line tools
- bash (Linux environment)
- smalt **v0.7.6**
- samtools **v1.9**
- bedtools (requires `bedtools multicov`)

### R
- Rscript **v4.5.2** (or compatible)
- The DESeq2 script installs missing R packages automatically (CRAN + Bioconductor):
  - CRAN: `dplyr`, `ggplot2`, `ggrepel`
  - Bioconductor: `DESeq2`, plus installer `BiocManager`

---

## Inputs

### FASTQ files
Place raw paired-end reads in:
- `./FASTQ/`

Supported naming patterns include:
- `*_R1.fastq.gz` / `*_R2.fastq.gz`
- `*_R1_001.fastq.gz` / `*_R2_001.fastq.gz`

---

## Usage

### Step 0 - Clone the repository and move into it:
```bash
git clone https://github.com/mikaOlami/Ld1S-DESeq2-Pipeline.git
cd Ld1S-DESeq2-Pipeline
```

### Step 1 - Map FASTQs to the Ld genome (SMALT) and generate sorted/indexed BAMs
From your working directory (the directory that contains `FASTQ/`), run:
```bash
bash scripts/LD_PolyA_Analysis.sh
```

Alternatively, use:
```bash
chmod +x scripts/LD_PolyA_Analysis.sh # only once
./scripts/LD_PolyA_Analysis.sh
```

Outputs (created in the working directory):
- Bams/*.sorted.bam
- Bams/*.sorted.bam.bai
- Logs/*.log

### Step 2 - Generate gene-level multicov table from BAMs
```bash
bash scripts/multiBamCov_Ld.sh
```

Output:
- `multicov_counts.tsv`

### Step 3 - Build DESeq2 inputs (counts matrix + colData template)
```bash
bash scripts/prepare_deseq_inputs.sh multicov_counts.tsv
```

Outputs:
- `counts.tsv` (gene x samples count matrix)
- `colData.tsv` (template to fill)
Edit `colData.tsv`:
  - set Condition (and optionally Batch) for each sample
  - to exclude samples, comment them out by prefixing the line with #
    - (commented lines are ignored when `colData.tsv` is read)

### Step 4 - Run DESeq2 comparison (produces PCA and results table)
```bash
Rscript scripts/DESEQ2_Analysis_Ld.R \
  counts.tsv colData.tsv \
  --cond CONDITION_NAME_AS_IN_COLDATA --ctrl CTRL_NAME_AS_IN_COLDATA \
  --cond_title CONDITION_TITLE --ctrl_title CTRL_TITLE
```

Outputs:
- `PCA_plot.png`
- `<cond_title>_vs_<ctrl_title>-results-with-normalized.csv` (DESeq2 results file)

### Notes
- The DESeq2 script uses only samples listed in `colData.tsv` (commented lines starting with `#` are ignored).
- Counts columns are automatically reordered to match the sample order in `colData.tsv`.
- Gene descriptions are added from `DB/LD_Id_Description.txt`.



