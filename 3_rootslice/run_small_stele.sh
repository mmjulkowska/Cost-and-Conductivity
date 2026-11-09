#!/usr/bin/env bash
# Run the small-stele fallback configs (3_rootslice/make_configs_smallstele.py) for the 19 replicates that hang or crash in RootSlice.
# For each replicate, tries variant A, then B, then C (120 s limit each) and keeps the first that finishes, as R5_<id>.
# Also runs 4 controls (variant A on replicates that ran with the full template) as CTRL_<id>.
# Raw run folders: $ROOTSLICE_HOME/results_costcond_v5_smallstele/. Collected volumes: <outputs>/rootslice/volume_analysis_smallstele.csv
# Usage:  bash 3_rootslice/run_small_stele.sh        (about 10-40 min)
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"; RS="${ROOTSLICE_HOME:-$HOME/root_models/RootSlice}"; CFG="$HERE/configs_smallstele"
OUTD="${COSTCOND_OUT:-$HERE/../outputs}/rootslice"; mkdir -p "$OUTD"
RES="$RS/results_costcond_v5_smallstele"; TRY="$RS/tries_costcond_v5_smallstele"; LOG="$OUTD/smallstele_log.txt"; LIMIT=${LIMIT:-120}
exec > >(tee -a "$LOG") 2>&1
echo "######## small-stele batch started $(date) ########"
[ -x "$RS/app/RootSlice" ] || { echo "RootSlice not found at $RS"; exit 1; }
mkdir -p "$RES" "$TRY"; : > "$OUTD/smallstele_variant_used.txt"
runone() {  # $1 config, $2 try dir -> 0 if finished
  local cfg="$1" out="$2"; rm -rf "$out"; mkdir -p "$out"; cp "$cfg" "$out/InputData.xml"
  (cd "$out" && exec "$RS/app/RootSlice" > run_stdout.txt 2>&1) & local pid=$! t=0
  while kill -0 $pid 2>/dev/null && [ $t -lt $LIMIT ]; do sleep 2; t=$((t+2)); done
  if kill -0 $pid 2>/dev/null; then kill $pid 2>/dev/null; sleep 1; kill -9 $pid 2>/dev/null; echo "   timed out (${LIMIT} s)"; return 1; fi
  wait $pid 2>/dev/null; find "$out" -name AllOutputName.csv | grep -q .
}
ids=$(ls "$CFG" | grep -v '^CTRL_' | sed 's/_[ABC]\.xml$//' | sort -u)
for id in $ids; do
  name="R5_$id"
  if find "$RES/$name" -name AllOutputName.csv 2>/dev/null | grep -q .; then echo "$id already done"; continue; fi
  ok=0
  for v in A B C; do
    echo "$id variant $v ..."; s=$(date +%s)
    if runone "$CFG/${id}_$v.xml" "$TRY/${id}_$v"; then
      rm -rf "$RES/$name"; cp -R "$TRY/${id}_$v" "$RES/$name"; sed -i '' "s/<Main name=\"[^\"]*\"/<Main name=\"$name\"/" "$RES/$name/InputData.xml" 2>/dev/null
      echo "$id $v" >> "$OUTD/smallstele_variant_used.txt"; echo "   ok ($(( $(date +%s)-s )) s)"; ok=1; break
    else echo "   variant $v failed"; fi
  done
  [ $ok = 1 ] || { echo "$id none" >> "$OUTD/smallstele_variant_used.txt"; echo "   ALL VARIANTS FAILED"; }
done
for c in "$CFG"/CTRL_*.xml; do
  id=$(basename "$c" _A.xml); echo "$id (control) ..."
  if runone "$c" "$TRY/$id"; then rm -rf "$RES/$id"; cp -R "$TRY/$id" "$RES/$id"; echo "   ok"; else echo "   control failed"; fi
done
(cd "$RS" && Rscript analysis/RootSlice_volume_analysis.R "$RES" >/dev/null 2>&1) && echo "R analysis done" || echo "R analysis FAILED"
cp "$RES/volume_analysis.csv" "$OUTD/volume_analysis_smallstele.csv" && echo "copied volume_analysis_smallstele.csv"
echo "variants used:"; cat "$OUTD/smallstele_variant_used.txt"
echo "######## small-stele batch finished $(date) ########"
