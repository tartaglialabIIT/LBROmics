#!/usr/bin/env bash
# Paths: override via environment variables (defaults are relative placeholders).
# Requires: bedtools, bigWigToBedGraph (UCSC), python3 + stats_peaks_vs_random.py
set -euo pipefail

# ============================================================
# Compare SAMMY signal inside peaks vs random matched intervals
# Random = bedtools shuffle on peaks (same #intervals, same lengths)
#
# Special-case for H3K9me3: GC-matched random sampling
#   - build a large candidate pool by shuffling peaks
#   - compute GC% of peaks and candidates with bedtools nuc
#   - bin peaks by GC deciles and sample matching counts from candidates
#
# Outputs per mark in OUTROOT/<MARK>/:
#   peaks.filtered.bed
#   sammy_on_peaks.tsv
#   random.filtered.seedXXXX.bed
#   sammy_on_random.seedXXXX.tsv
#   stats.seedXXXX.txt
#   plot.seedXXXX.pdf
# ============================================================

# ----------------------------
# INPUT PATHS (from your convo)
# ----------------------------
SAMMY_DIR="${SAMMY_DIR:-./spp_mle}"
BW1="${SAMMY_DIR}/NPCwt_n1_S2SvsS3_mle.bigWig"
BW2="${SAMMY_DIR}/NPCwt_n2_S2SvsS3_mle.bigWig"
BW3="${SAMMY_DIR}/NPCwt_IRE_R3_S2SvsS3_mle.bigWig"

H3K27AC_PEAKS="${H3K27AC_PEAKS:-./peaks/H3K27ac.consensus_peaks.bed}"
H3K4ME3_PEAKS="${H3K4ME3_PEAKS:-./peaks/H3K4me3.consensus_peaks.bed}"
H3K9ME3_PEAKS="${H3K9ME3_PEAKS:-./peaks/H3K9me3.consensus_peaks.bed}"

GENOME_SIZES="${GENOME_SIZES:-./mm10.chrom.sizes.from_bam}"

# mm10 FASTA for GC matching (bedtools nuc)
MM10_FASTA="${MM10_FASTA:-./mm10.fa}"

# optional blacklist (leave empty "" if none)
BLACKLIST_BED=""

# Output root
OUTROOT="${OUTROOT:-./sammy_vs_peaks_stats}"
mkdir -p "$OUTROOT"

# Randomization
N_SEEDS=20
SEED0=1001

# deepTools threads
NPROC=24

# Python stats script
STATS_PY="stats_peaks_vs_random.py"

# GC-matching candidate pool multiplier for H3K9me3
GC_CAND_MULT=10

# Which marks to run (set to "all" or "H3K9me3" etc.)
RUN_WHICH="all"

# ----------------------------
# Checks
# ----------------------------
for f in "$BW1" "$BW2" "$BW3" "$GENOME_SIZES" "$STATS_PY"; do
  [[ -s "$f" ]] || { echo "ERROR: missing file: $f" >&2; exit 1; }
done
command -v bedtools >/dev/null || { echo "ERROR: bedtools not found in PATH" >&2; exit 1; }
command -v multiBigwigSummary >/dev/null || { echo "ERROR: multiBigwigSummary not found in PATH" >&2; exit 1; }
command -v python3 >/dev/null || { echo "ERROR: python3 not found in PATH" >&2; exit 1; }

if [[ -n "$BLACKLIST_BED" ]]; then
  [[ -s "$BLACKLIST_BED" ]] || { echo "ERROR: blacklist set but missing: $BLACKLIST_BED" >&2; exit 1; }
fi

# GC matching requires FASTA
if [[ ! -s "$MM10_FASTA" ]]; then
  echo "WARNING: MM10_FASTA missing ($MM10_FASTA). H3K9me3 will fall back to plain shuffle." >&2
fi

# ----------------------------
# Helpers
# ----------------------------
filter_bed() {
  # Keep only standard chromosomes, drop weird contigs; output ONLY first 3 columns
  local inbed="$1"
  local outbed="$2"
  awk 'BEGIN{OFS="\t"}
       $1 ~ /^chr/ &&
       $1 != "chrM" &&
       $1 !~ /_random$/ &&
       $1 !~ /^chrUn/ &&
       $1 !~ /^Un_/ &&
       $1 !~ /alt$/ &&
       $1 !~ /fix$/ {print $1,$2,$3}' "$inbed" \
    | sort -k1,1 -k2,2n > "$outbed"
}

extract_sammy() {
  local bed="$1"
  local outprefix="$2"
  multiBigwigSummary BED-file \
    -b "$BW1" "$BW2" "$BW3" \
    --BED "$bed" \
    --labels WT1 WT2 WT3 \
    --outFileName "${outprefix}.npz" \
    --outRawCounts "${outprefix}.tsv" \
    -p "$NPROC"
}

make_excl_bed() {
  local peaks="$1"
  local excl_out="$2"
  if [[ -n "$BLACKLIST_BED" ]]; then
    cat "$peaks" "$BLACKLIST_BED" \
      | awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' \
      | sort -k1,1 -k2,2n > "$excl_out"
  else
    cp "$peaks" "$excl_out"
  fi
}

make_random_shuffle() {
  local peaks="$1"
  local excl="$2"
  local seed="$3"
  local outbed="$4"
  local noOverlap="$5"   # "true" or "false"

  local extra=()
  if [[ "$noOverlap" == "true" ]]; then
    extra+=("-noOverlapping")
  fi

  bedtools shuffle \
    -i "$peaks" \
    -g "$GENOME_SIZES" \
    -excl "$excl" \
    -chrom \
    -seed "$seed" \
    -maxTries 20000 \
    "${extra[@]}" \
  | awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' \
  | sort -k1,1 -k2,2n > "$outbed"
}

run_stats() {
  local peaks_tsv="$1"
  local rand_tsv="$2"
  local mark="$3"
  local out_txt="$4"
  local out_pdf="$5"
  python3 "$STATS_PY" "$peaks_tsv" "$rand_tsv" "$mark" "$out_txt" "$out_pdf"
}

ensure_fai() {
  # Create FASTA index once if missing
  local fa="$1"
  local fai="${fa}.fai"
  if [[ -s "$fa" && ! -s "$fai" ]]; then
    if command -v samtools >/dev/null; then
      echo "[INFO] FASTA index missing; creating: $fai"
      samtools faidx "$fa"
    else
      echo "WARNING: samtools not found; bedtools may create .fai itself (slower / sometimes buggy)." >&2
    fi
  fi
}

# ---- H3K9me3-specific: GC-matched random sampling ----
# Notes:
# - Works with bedtools v2.27.1 header "8_pct_gc"
# - We parse the header robustly (do NOT rely on pandas comment='#')
make_random_gc_matched_H3K9me3() {
  local peaks="$1"
  local excl="$2"
  local seed="$3"
  local outbed="$4"

  if [[ ! -s "$MM10_FASTA" ]]; then
    echo "[WARN] No FASTA; fallback to plain shuffle for H3K9me3" >&2
    make_random_shuffle "$peaks" "$excl" "$seed" "$outbed" "false"
    return 0
  fi

  ensure_fai "$MM10_FASTA"

  local tmpdir
  tmpdir="$(mktemp -d)"
  trap 'rm -rf "'"$tmpdir"'"' RETURN

  local n_peaks
  n_peaks=$(wc -l < "$peaks")
  local cand_n=$((n_peaks * GC_CAND_MULT))

  local cand_raw="$tmpdir/candidates.raw.bed"
  local cand_bed="$tmpdir/candidates.bed"

  : > "$cand_raw"
  local tries=0
  while [[ $(wc -l < "$cand_raw") -lt "$cand_n" ]]; do
    tries=$((tries+1))
    local s=$((seed + tries*37))
    bedtools shuffle \
      -i "$peaks" \
      -g "$GENOME_SIZES" \
      -excl "$excl" \
      -chrom \
      -seed "$s" \
      -maxTries 20000 \
    >> "$cand_raw"

    if [[ "$tries" -gt 80 ]]; then
      echo "ERROR: cannot generate enough candidates for GC matching (tries=$tries)." >&2
      exit 1
    fi
  done

  head -n "$cand_n" "$cand_raw" \
    | awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' \
    | sort -k1,1 -k2,2n > "$cand_bed"

  # Compute nuc tables
  bedtools nuc -fi "$MM10_FASTA" -bed "$peaks"    > "$tmpdir/peaks.nuc.tsv"
  bedtools nuc -fi "$MM10_FASTA" -bed "$cand_bed" > "$tmpdir/cand.nuc.tsv"

  # Build GC-matched sample
  python3 - <<'PY' "$tmpdir/peaks.nuc.tsv" "$tmpdir/cand.nuc.tsv" "$seed" "$outbed"
import sys, re
import numpy as np
import pandas as pd

peaks_nuc, cand_nuc, seed_s, outbed = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
seed = int(seed_s)

def read_bedtools_nuc(path: str) -> pd.DataFrame:
    # Read first line as header (starts with "#1_usercol ... 8_pct_gc ...")
    with open(path, "r") as f:
        header = f.readline().rstrip("\n")
    cols = header.split("\t")
    cols = [c.lstrip("#").strip() for c in cols]
    df = pd.read_csv(path, sep="\t", header=None, skiprows=1)
    if df.shape[1] != len(cols):
        # sometimes bedtools collapses columns visually; handle by re-splitting on whitespace
        df = pd.read_csv(path, sep=r"\t+|\s+", engine="python", header=None, skiprows=1)
    df.columns = cols[:df.shape[1]]
    return df

peaks = read_bedtools_nuc(peaks_nuc)
cand  = read_bedtools_nuc(cand_nuc)

# Find GC column (bedtools v2.27.1: "8_pct_gc")
gc_candidates = [c for c in peaks.columns if "pct_gc" in c.lower()]
if not gc_candidates:
    raise SystemExit(f"ERROR: cannot find GC column. Columns: {list(peaks.columns)}")
gc_col = gc_candidates[0]

# Identify coord columns (first 3 user columns)
# In v2.27.1 they are "1_usercol","2_usercol","3_usercol"
c1, c2, c3 = peaks.columns[0], peaks.columns[1], peaks.columns[2]

peaks_gc = pd.to_numeric(peaks[gc_col], errors="coerce").to_numpy()
cand_gc  = pd.to_numeric(cand[gc_col],  errors="coerce").to_numpy()

# Drop rows with bad GC
peaks = peaks[np.isfinite(peaks_gc)]
cand  = cand[np.isfinite(cand_gc)]

peaks_gc = pd.to_numeric(peaks[gc_col], errors="coerce").to_numpy()
cand_gc  = pd.to_numeric(cand[gc_col],  errors="coerce").to_numpy()

n_peaks = len(peaks)
if n_peaks == 0 or len(cand) == 0:
    raise SystemExit("ERROR: empty peaks or candidate GC table after filtering.")

# GC deciles based on peaks distribution
qs = np.quantile(peaks_gc, np.linspace(0, 1, 11))
qs[0]  = -1.0
qs[-1] =  2.0

peaks["gc_bin"] = pd.cut(peaks_gc, bins=qs, include_lowest=True, labels=False)
cand["gc_bin"]  = pd.cut(cand_gc,  bins=qs, include_lowest=True, labels=False)

rng = np.random.default_rng(seed)
out_parts = []

# For each GC bin, sample matching count
for b in range(10):
    need = int((peaks["gc_bin"] == b).sum())
    if need == 0:
        continue
    pool = cand[cand["gc_bin"] == b]
    if len(pool) < need:
        # if pool too small in that bin, fall back to whole candidate set
        pool = cand
    take = pool.sample(n=need, replace=False, random_state=seed + b*101)
    out_parts.append(take[[c1, c2, c3]])

out = pd.concat(out_parts, axis=0)

# If for any reason we got slightly off (shouldn't), trim/pad
if len(out) > n_peaks:
    out = out.sample(n=n_peaks, random_state=seed)
elif len(out) < n_peaks:
    # last resort: sample remaining from candidate set
    rem = n_peaks - len(out)
    extra = cand[[c1,c2,c3]].sample(n=rem, replace=False, random_state=seed+999)
    out = pd.concat([out, extra], axis=0)

out = out.rename(columns={c1:"chr", c2:"start", c3:"end"}).copy()
out["start"] = out["start"].astype(int)
out["end"]   = out["end"].astype(int)
out = out.sort_values(["chr","start","end"])

out.to_csv(outbed, sep="\t", header=False, index=False)
print(f"[OK] GC-matched random written: {outbed} (n={len(out)})", file=sys.stderr)
print(f"[OK] GC column used: {gc_col}", file=sys.stderr)
PY

  # sanity count
  local n_out
  n_out=$(wc -l < "$outbed")
  if [[ "$n_out" -ne "$n_peaks" ]]; then
    echo "ERROR: GC-matched random count mismatch ($n_out != $n_peaks)" >&2
    exit 1
  fi
}

run_mark() {
  local mark="$1"
  local peaks="$2"

  [[ -s "$peaks" ]] || { echo "SKIP: missing peaks for $mark: $peaks"; return 0; }

  local outdir="${OUTROOT}/${mark}"
  mkdir -p "$outdir"

  echo "=============================="
  echo "MARK: $mark"
  echo "Peaks: $peaks"
  echo "Out:   $outdir"
  echo "=============================="

  local peaks_filt="${outdir}/peaks.filtered.bed"
  filter_bed "$peaks" "$peaks_filt"
  local n_peaks
  n_peaks=$(wc -l < "$peaks_filt")
  echo "[OK] Filtered peaks: $peaks_filt (n=$n_peaks)"

  [[ "$n_peaks" -ge 50 ]] || { echo "ERROR: too few peaks after filtering for $mark (n=$n_peaks)" >&2; return 1; }

  local peaks_prefix="${outdir}/sammy_on_peaks"
  extract_sammy "$peaks_filt" "$peaks_prefix"
  echo "[OK] Extracted SAMMY on peaks -> ${peaks_prefix}.tsv"

  local excl="${outdir}/exclusion.tmp.bed"
  make_excl_bed "$peaks_filt" "$excl"

  for i in $(seq 0 $((N_SEEDS-1))); do
    local seed=$((SEED0 + i))
    local rand_bed="${outdir}/random.filtered.seed${seed}.bed"

    if [[ "$mark" == "H3K9me3" ]]; then
      echo "[INFO] H3K9me3: using GC-matched random (seed=$seed)"
      make_random_gc_matched_H3K9me3 "$peaks_filt" "$excl" "$seed" "$rand_bed"
    else
      make_random_shuffle "$peaks_filt" "$excl" "$seed" "$rand_bed" "true"
      local n_rand
      n_rand=$(wc -l < "$rand_bed")
      if [[ "$n_rand" -lt "$n_peaks" ]]; then
        echo "[WARN] shuffle(noOverlapping) produced fewer ($n_rand < $n_peaks). Retrying allowing overlap."
        make_random_shuffle "$peaks_filt" "$excl" "$seed" "$rand_bed" "false"
      fi
    fi

    local n_rand
    n_rand=$(wc -l < "$rand_bed")
    if [[ "$n_rand" -lt "$n_peaks" ]]; then
      echo "ERROR: fewer random intervals than peaks ($n_rand < $n_peaks) for seed=$seed" >&2
      exit 1
    fi
    echo "[OK] Random matched bed (seed=$seed): $rand_bed (n=$n_rand)"

    local rand_prefix="${outdir}/sammy_on_random.seed${seed}"
    extract_sammy "$rand_bed" "$rand_prefix"
    echo "[OK] Extracted SAMMY on random -> ${rand_prefix}.tsv"

    local stats_txt="${outdir}/stats.seed${seed}.txt"
    local plot_pdf="${outdir}/plot.seed${seed}.pdf"
    run_stats "${peaks_prefix}.tsv" "${rand_prefix}.tsv" "$mark" "$stats_txt" "$plot_pdf"
    echo "[OK] Stats: $stats_txt"
    echo "[OK] Plot:  $plot_pdf"
  done

  rm -f "$excl"
  echo "[DONE] $mark"
}

# ----------------------------
# Run
# ----------------------------
case "$RUN_WHICH" in
  all)
    run_mark "H3K27ac" "$H3K27AC_PEAKS"
    run_mark "H3K4me3" "$H3K4ME3_PEAKS"
    run_mark "H3K9me3" "$H3K9ME3_PEAKS"
    ;;
  H3K27ac) run_mark "H3K27ac" "$H3K27AC_PEAKS" ;;
  H3K4me3) run_mark "H3K4me3" "$H3K4ME3_PEAKS" ;;
  H3K9me3) run_mark "H3K9me3" "$H3K9ME3_PEAKS" ;;
  *) echo "ERROR: RUN_WHICH must be one of: all, H3K27ac, H3K4me3, H3K9me3" >&2; exit 1 ;;
esac

echo "All done. Results in: $OUTROOT"
