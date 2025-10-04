#!/bin/bash
set -euo pipefail

ASSEMBLY_DIR="assemblies"
REF_DB="ref_mash_db.msh"
THRESHOLD=0.001 # should be > 100/<est_genome_len>
OUT_CSV="passed_dist_filter.csv"
TMP_DIR=$(mktemp -d)
PARALLEL=4   # number of parallel jobs
LOCKFILE="$TMP_DIR/prog.lock"
COUNTER_FILE="$TMP_DIR/progress.count"

rm -rf $OUT_CSV

cleanup() {
  rm -rf "$TMP_DIR"
}
trap cleanup EXIT

# Prepare output and progress
echo "assembly_file,mash_distance" > "$OUT_CSV"
echo 0 > "$COUNTER_FILE"

# Build array of assembly files (handles spaces/newlines) and count reliably
mapfile -d '' -t ASMS < <(find "$ASSEMBLY_DIR" -maxdepth 1 -type f -print0)
TOTAL=${#ASMS[@]}

if [ "$TOTAL" -eq 0 ]; then
  echo "No files found in $ASSEMBLY_DIR"
  exit 0
fi

# Progress display (runs in background)
show_progress() {
  while :; do
    if [ ! -f "$COUNTER_FILE" ]; then break; fi
    count=$(cat "$COUNTER_FILE")
    # ensure count is an integer
    count=${count:-0}
    if [ "$TOTAL" -gt 0 ]; then
      percent=$(( 100 * count / TOTAL ))
    else
      percent=0
    fi
    # build bar (50 chars)
    filled=$(( percent * 50 / 100 ))
    empty=$(( 50 - filled ))
    # produce filled and empty segments without seq when zero
    if [ "$filled" -gt 0 ]; then
      filled_seg=$(printf '#%.0s' $(seq 1 $filled))
    else
      filled_seg=""
    fi
    if [ "$empty" -gt 0 ]; then
      empty_seg=$(printf ' %.0s' $(seq 1 $empty))
    else
      empty_seg=""
    fi
    printf '\rProgress: [%s%s] %d/%d (%d%%)' "$filled_seg" "$empty_seg" "$count" "$TOTAL" "$percent"
    if [ "$count" -ge "$TOTAL" ]; then
      printf '\n'
      break
    fi
    sleep 0.5
  done
}
show_progress & progress_pid=$!

# Worker function for xargs
process_one() {
  asm="$1"
  asm_base=$(basename "$asm")
  sketch_prefix="$TMP_DIR/${asm_base%.*}"
  sketch_file="${sketch_prefix}.msh"
  dist_file="$TMP_DIR/${asm_base}.dist"

  mash sketch -o "$sketch_prefix" "$asm" >/dev/null 2>&1
  mash dist "$REF_DB" "$sketch_file" > "$dist_file"

  awk -v asm="$asm_base" -v t="$THRESHOLD" '
    $3 <= t { printf "%s,%.6f\n", asm, $3; exit }
  ' "$dist_file" >> "$OUT_CSV"

  # increment progress counter atomically
  (
    flock 200
    if [ ! -f "$COUNTER_FILE" ]; then
      echo 0 > "$COUNTER_FILE"
    fi
    cnt=$(cat "$COUNTER_FILE")
    cnt=$((cnt + 1))
    echo "$cnt" > "$COUNTER_FILE"
  ) 200>"$LOCKFILE"
}
export -f process_one

export REF_DB THRESHOLD TMP_DIR OUT_CSV LOCKFILE COUNTER_FILE

# Run workers via xargs using the array of files to avoid re-running find
printf '%s\0' "${ASMS[@]}" | xargs -0 -n1 -P "$PARALLEL" -I{} bash -c 'process_one "$@"' _ {}

# Wait for progress background job to finish
wait "$progress_pid"

echo "Filtered results saved to $OUT_CSV"
