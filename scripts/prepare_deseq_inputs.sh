#!/usr/bin/env bash
set -euo pipefail

IN="${1:-multicov_counts.tsv}"
COUNTS_OUT="${2:-counts.tsv}"
COLDATA_OUT="${3:-colData.tsv}"

[[ -f "$IN" ]] || { echo "Error: input not found: $IN" >&2; exit 1; }

# multicov_counts.tsv columns:
# chrom start end gene score strand sample1 sample2 ...

awk -F'\t' 'BEGIN{OFS="\t"}
NR==1{
  # header
  printf "gene";
  for (i=7; i<=NF; i++) printf OFS $i;
  printf "\n";
  next
}
{
  printf $4;
  for (i=7; i<=NF; i++) printf OFS $i;
  printf "\n";
}' "$IN" > "$COUNTS_OUT"

# colData template from the header sample names
# Default Condition=TODO, Batch=1
{
  echo -e "SampleID\tCondition\tBatch"
  awk -F'\t' 'NR==1{for(i=7;i<=NF;i++) print $i"\tTODO\t1"}' "$IN"
} > "$COLDATA_OUT"

echo "Wrote: $COUNTS_OUT"
echo "Wrote: $COLDATA_OUT (edit Condition/Batch; comment out samples with #)"
