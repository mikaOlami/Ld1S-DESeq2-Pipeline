#!/usr/bin/env bash
set -euo pipefail

# L. donovani (Ld1S) PolyA RNA-seq mapping pipeline
# FASTQ/ (paired-end .fastq.gz)  ->  Bams/*.sorted.bam (+ .bai)
#
# Repo layout expected:
#   <repo>/scripts/LD_PolyA_Analysis.sh   (this script)
#   <repo>/DB/LD_smalt_index.smi/.sma     (SMALT genome index; prefix: DB/LD_smalt_index)
#
# Run from a working directory that contains FASTQ/:
#   bash /path/to/repo/scripts/LD_PolyA_Analysis.sh
#
# Optionally override tools/threads:
#   threads=16 max_jobs=10 SMALT=smalt SAMTOOLS=samtools bash scripts/LD_PolyA_Analysis.sh
#
# Optionally override repo root if auto-detection fails:
#   REPO_DIR=/path/to/repo bash scripts/LD_PolyA_Analysis.sh

# ---- Config ----
threads="${threads:-8}"
max_jobs="${max_jobs:-25}"   # how many samples to process concurrently

# ---- Resolve repo root (so the script works from any working directory) ----
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${REPO_DIR:-"$(cd -- "${SCRIPT_DIR}/.." && pwd)"}"

# ---- Reference paths inside repo ----
LD_DB="${REPO_DIR}/DB/LD_smalt_index"   # SMALT index prefix (expects .smi/.sma)

# ---- Tools (assumed to be on PATH unless overridden) ----
SMALT="${SMALT:-smalt}"
SAMTOOLS="${SAMTOOLS:-samtools}"

# ---- Sanity ----
command -v "$SMALT" >/dev/null 2>&1 || { echo "Error: smalt not found on PATH (or SMALT not set)" >&2; exit 1; }
command -v "$SAMTOOLS" >/dev/null 2>&1 || { echo "Error: samtools not found on PATH (or SAMTOOLS not set)" >&2; exit 1; }

[[ -f "${LD_DB}.smi" && -f "${LD_DB}.sma" ]] || {
  echo "Error: SMALT index not found at:" >&2
  echo "  ${LD_DB}.smi / ${LD_DB}.sma" >&2
  exit 1
}

[[ -d FASTQ ]] || {
  echo -e "Error: FASTQ directory does not exist in the current working directory.\n" \
          "Please create FASTQ/ and place paired-end *.fastq.gz files inside it." >&2
  exit 1
}

mkdir -p Bams Logs

# ---- simple job limiter ----
joblim_wait() {
  while (( $(jobs -rp | wc -l) >= max_jobs )); do
    wait -n || true
  done
}

# ---- per-sample pipeline ----
process_sample() {
  local r1="$1"
  local threads="$2"

  # derive R2 (supports *_R1.fastq.gz and *_R1_001.fastq.gz)
  local r1_base r2_base r2
  r1_base="$(basename "$r1")"
  r2_base="$(
    echo "$r1_base" | sed -E 's/(.*)_R1(_[0-9]+)?\.fastq\.gz$/\1_R2\2.fastq.gz/'
  )"
  r2="FASTQ/$r2_base"

  if [[ ! -e "$r2" ]]; then
    echo "[WARN] Skipping: paired file not found for $r1 -> expected $r2" >&2
    return 0
  fi

  local base
  base="$(
    echo "$r1_base" | sed -E 's/_R1(_[0-9]+)?\.fastq\.gz$//'
  )"

  local ubam="Bams/${base}.bam"
  local ubam_tmp="${ubam}.tmp"
  local sbam="Bams/${base}.sorted.bam"
  local sbam_tmp="${sbam}.tmp"
  local sbai="${sbam}.bai"

  local log_map="Logs/${base}.smalt.log"
  local log_sort="Logs/${base}.sort.log"

  echo "[START] $base"

  # --- MAP ---
  if [[ ! -s "$ubam" || "$ubam" -ot "$r1" || "$ubam" -ot "$r2" ]]; then
    echo "[MAP: Ld PolyA vs genome] $base"
    "$SMALT" map -n "${threads}" "${LD_DB}" "$r1" "$r2" 2> "$log_map" \
      | "$SAMTOOLS" view -@ "${threads}" -b -f 0x02 -F 4 -o "$ubam_tmp" -
    "$SAMTOOLS" quickcheck -v "$ubam_tmp"
    mv -f "$ubam_tmp" "$ubam"
    [[ -s "$log_map" ]] || rm -f "$log_map"
  else
    echo "[MAP] SKIP (up-to-date) $base"
  fi

  [[ -s "$ubam" ]] || { echo "[ERROR] BAM missing for $base" >&2; return 1; }

  # --- SORT ---
  if [[ ! -s "$sbam" || "$sbam" -ot "$ubam" ]]; then
    echo "[SORT] Sorting $base"
    "$SAMTOOLS" sort -@ "${threads}" -o "$sbam_tmp" "$ubam" 2> "$log_sort"
    "$SAMTOOLS" quickcheck -v "$sbam_tmp"
    mv -f "$sbam_tmp" "$sbam"
    rm -f "$ubam"  # delete unsorted BAM once sorted is verified
    [[ -s "$log_sort" ]] || rm -f "$log_sort"
  else
    echo "[SORT] SKIP (up-to-date) $base"
  fi

  # --- INDEX ---
  if [[ ! -s "$sbai" || "$sbai" -ot "$sbam" ]]; then
    echo "[INDEX] Indexing $base"
    # Note: -@ may be ignored by older samtools; harmless if unsupported
    "$SAMTOOLS" index -@ "${threads}" "$sbam" 2>> "$log_sort" || "$SAMTOOLS" index "$sbam"
  else
    echo "[INDEX] SKIP (up-to-date) $base"
  fi

  echo "[DONE] $base"
}

# ---- discover & process ----
shopt -s nullglob
found_any=false

for r1 in FASTQ/*_R1_001.fastq.gz FASTQ/*_R1.fastq.gz; do
  [[ -e "$r1" ]] || continue
  found_any=true
  joblim_wait
  (
    set -euo pipefail
    process_sample "$r1" "$threads"
  ) &
done

if ! $found_any; then
  echo "No FASTQ/*_R1*.fastq.gz files found in FASTQ/" >&2
  exit 0
fi

wait
echo "[STAGE] All per-sample chains finished."

# clean up empty logs
find Logs -type f -name "*.log" -size 0 -delete 2>/dev/null || true
echo "SCRIPT FINISHED."
