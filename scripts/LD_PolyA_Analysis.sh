#!/usr/bin/env bash
set -euo pipefail

# L. donovani (Ld1S) PolyA RNA-seq mapping pipeline
# FASTQ/ (paired-end .fastq.gz)  ->  Bams/*.sorted.bam (+ .bai)
#
#
# Run from a working directory that contains FASTQ/:
#   bash /path/to/repo/scripts/LD_PolyA_Analysis.sh
#
# Optional overrides:
#   threads=16 max_jobs=10 SMALT=smalt SAMTOOLS=samtools bash scripts/LD_PolyA_Analysis.sh
#   REPO_DIR=/path/to/repo bash scripts/LD_PolyA_Analysis.sh

# ---- Config ----
threads="${threads:-8}"
max_jobs="${max_jobs:-25}"

# ---- Resolve repo root (works from any working directory) ----
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${REPO_DIR:-"$(cd -- "${SCRIPT_DIR}/.." && pwd)"}"
DB_DIR="${REPO_DIR}/DB"

# ---- Tools (assumed on PATH unless overridden) ----
SMALT="${SMALT:-smalt}"
SAMTOOLS="${SAMTOOLS:-samtools}"

command -v "$SMALT" >/dev/null 2>&1 || { echo "Error: smalt not found on PATH (or SMALT not set)" >&2; exit 1; }
command -v "$SAMTOOLS" >/dev/null 2>&1 || { echo "Error: samtools not found on PATH (or SAMTOOLS not set)" >&2; exit 1; }

# ---- Reference FASTA + SMALT index prefix (inside repo DB/) ----
LD_FASTA="${DB_DIR}/Leishmania_donovani_sudanese.fa"
LD_DB="${DB_DIR}/LD_smalt_index"   # prefix; expects ${LD_DB}.smi and ${LD_DB}.sma

[[ -d "$DB_DIR" ]] || { echo "Error: DB directory not found: $DB_DIR" >&2; exit 1; }
[[ -f "$LD_FASTA" ]] || { echo "Error: genome FASTA not found: $LD_FASTA" >&2; exit 1; }

# ---- Build SMALT index if missing ----
if [[ ! -f "${LD_DB}.smi" || ! -f "${LD_DB}.sma" ]]; then
  echo "[DB] SMALT index not found. Building it now..."

  # Build inside DB so outputs land in DB/
  (
    cd "$DB_DIR"
    "$SMALT" index -k 11 -s 1 "$(basename "$LD_DB")" "$(basename "$LD_FASTA")"
  )

  [[ -f "${LD_DB}.smi" && -f "${LD_DB}.sma" ]] || {
    echo "Error: SMALT index build failed (missing ${LD_DB}.smi/.sma)" >&2
    exit 1
  }
  echo "[DB] SMALT index built successfully."
else
  echo "[DB] SMALT index found: ${LD_DB}.smi/.sma"
fi

# ---- Sanity: working directory must contain FASTQ/ ----
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

  local r1_base r2_base r2
  r1_base="$(basename "$r1")"

  # derive R2 (supports *_R1.fastq.gz and *_R1_001.fastq.gz)
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
    echo "[MAP: LD PolyA vs genome] $base"
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
