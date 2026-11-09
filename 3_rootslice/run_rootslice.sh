#!/usr/bin/env bash
# Run RootSlice (Jagdeep's Apple Silicon build) on the 216 configs in 3_rootslice/configs/, then Jagdeep's R volume analysis.
# Raw run folders: $ROOTSLICE_HOME/results_costcond_v5/ (archive these on Box).
# Collected volumes: <outputs>/rootslice/volume_analysis_configs.csv ; failures: <outputs>/rootslice/failed_runs.txt
# Usage:  bash 3_rootslice/run_rootslice.sh        (about 1.5 h; safe to re-run, finished runs are skipped)
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"; RS="${ROOTSLICE_HOME:-$HOME/root_models/RootSlice}"
OUTD="${COSTCOND_OUT:-$HERE/../outputs}/rootslice"; mkdir -p "$OUTD"; LOG="$OUTD/rootslice_log.txt"
exec > >(tee -a "$LOG") 2>&1
echo "######## RootSlice batch started $(date) ########"
[ -x "$RS/app/RootSlice" ] || { echo "RootSlice not found at $RS - see env/setup_rootslice_mac.sh and data/README.md"; exit 1; }
SETS="configs"
for set in $SETS; do
  RES="$RS/results_costcond_v5${RS_TAG:-}"; mkdir -p "$RES"   # RS_TAG: optional suffix for a fresh results folder
  n=0; total=$(ls "$HERE/$set"/*.xml | wc -l | tr -d ' ')
  for cfg in "$HERE/$set"/*.xml; do
    n=$((n+1)); name="$(basename "$cfg" .xml)"
    if find "$RES/$name" -name AllOutputName.csv 2>/dev/null | grep -q .; then echo "[$n/$total] $name already done"; continue; fi
    echo "[$n/$total] $name ..."
    out="$RES/$name"; mkdir -p "$out"; cp "$cfg" "$out/InputData.xml"
    start=$(date +%s)
    ok=0
    for attempt in 1 2; do   # time limit: RootSlice can hang on some geometries; kill after LIMIT s, retry once
      (cd "$out" && exec "$RS/app/RootSlice" > run_stdout.txt 2>&1) & pid=$!
      t=0; while kill -0 $pid 2>/dev/null && [ $t -lt ${LIMIT:-180} ]; do sleep 2; t=$((t+2)); done
      if kill -0 $pid 2>/dev/null; then kill $pid 2>/dev/null; sleep 1; kill -9 $pid 2>/dev/null; echo "   attempt $attempt timed out after ${LIMIT:-180} s"; continue; fi
      wait $pid && find "$out" -name AllOutputName.csv | grep -q . && { ok=1; break; }
    done
    if [ $ok = 1 ]; then echo "   ok ($(( $(date +%s)-start )) s)"
    else echo "   FAILED/TIMED OUT - skipped"; echo "$name" >> "$OUTD/failed_runs.txt"; tail -3 "$out/run_stdout.txt"; fi
  done
  (cd "$RS" && Rscript analysis/RootSlice_volume_analysis.R "$RES" >/dev/null 2>&1) && echo "R analysis done for $set" || echo "R analysis FAILED for $set"
  cp "$RES/volume_analysis.csv" "$OUTD/volume_analysis_configs.csv" && echo "copied volume_analysis_configs.csv"
done
echo "######## RootSlice batch finished $(date) ########"
