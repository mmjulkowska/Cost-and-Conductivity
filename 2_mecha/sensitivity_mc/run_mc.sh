#!/usr/bin/env bash
# Monte Carlo sensitivity of kr (V. Couvreur's suggestion). Runs the three designs in turn; resumable (re-run the same command).
#   bash 2_mecha/sensitivity_mc/run_mc.sh     (several hours; results in <outputs>/mecha/sensitivity_mc/)
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"; OUTD="${COSTCOND_OUT:-$HERE/../../outputs}/mecha/sensitivity_mc"; mkdir -p "$OUTD"
exec > >(tee -a "$OUTD/mc_log.txt") 2>&1
for c in /opt/anaconda3/bin/conda /opt/homebrew/Caskroom/miniforge/base/bin/conda "$HOME/miniforge3/bin/conda"; do
  [ -x "$c" ] && eval "$("$c" shell.bash hook)" && break; done
echo "######## MC started $(date), $(sysctl -n hw.ncpu 2>/dev/null || nproc) CPUs ########"
for d in "rep 2" "rep 10" "all 2"; do
  conda run --no-capture-output -n mecha_dev python -W ignore "$HERE/mc_kr.py" $d 2>&1 | grep -v -i warning
done
echo "######## MC finished $(date) ########"
conda run --no-capture-output -n mecha_dev python "$HERE/analyse_mc.py" && conda run --no-capture-output -n mecha_dev python "$HERE/fig_mc.py"
