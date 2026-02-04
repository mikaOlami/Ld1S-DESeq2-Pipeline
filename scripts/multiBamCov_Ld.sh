#!/usr/bin/env bash
set -euo pipefail

# Repo-root discovery (this script lives in scripts/)
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

DEFAULT_BED="${REPO_DIR}/DB/LD_mRNAs_w_merged_w_UTRs.bed"
OUT="multicov_counts.tsv"

BED=""

usage() {
  echo "Usage: $0 [BED_FILE] [-o output.tsv]" >&2
  exit 1
}

while (( "$#" )); do
  case "$1" in
    -o) shift; OUT="${1:-}"; [[ -n "$OUT" ]] || usage; shift ;;
    -h|--help) usage ;;
    -*)
      echo "Unknown option: $1" >&2
      usage
      ;;
    *)
      if [[ -z "$BED" ]]; then
        BED="$1"; shift
      else
        echo "Unexpected extra argument: $1" >&2
        usage
      fi
      ;;
  esac
done

BED="${BED:-$DEFAULT_BED}"

command -v samtools >/dev/null 2>&1 || { echo "Error: samtools not found on PATH" >&2; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo "Error: bedtools not found on PATH" >&2; exit 1; }

[[ -f "$BED" ]] || { echo "Error: BED file not found: $BED" >&2; exit 1; }

# Find BAMs
BAM_DIR="."
[[ -d "Bams" ]] && BAM_DIR="Bams"

mapfile -t BAMS < <(find "$BAM_DIR" -maxdepth 1 -type f -name "*.bam" | sort)
((${#BAMS[@]} > 0)) || { echo "No BAM files found in $BAM_DIR" >&2; exit 1; }

echo -e "Using BED:\n  $BED"
echo -e "Writing:\n  $OUT"

# Ensure indexes exist
for b in "${BAMS[@]}"; do
  [[ -e "${b}.bai" ]] || samtools index "$b"
done

# Write header + multicov output
{
  printf "chrom\tstart\tend\tgene\tscore\tstrand"
  for b in "${BAMS[@]}"; do
    printf "\t%s" "$(basename "${b%.bam}")"
  done
  printf "\n"
  bedtools multicov -bed "$BED" -bams "${BAMS[@]}"
} > "$OUT"

echo "Done: $OUT"