#!/usr/bin/env bash
# Rerun the whole pipeline from Emily Miller's measurement sheet and compare with the frozen reference.
# Before:  conda activate mecha_dev  (env/install_mecha_mac.sh), RootSlice installed (env/setup_rootslice_mac.sh),
#          data/AllDatammAndpixelsNov2026_2.xlsx downloaded (and data/Cross Sections/ for stage 5), see data/README.md.
# Usage:   bash reproduce/run_all.sh                 default stages: 1 2 3 4 5 6 8   (about 2 h, mostly RootSlice)
#          bash reproduce/run_all.sh 1 2             only GRANAR and MECHA (~10 min)
#          bash reproduce/run_all.sh 7               Monte Carlo of kr (several hours; optional)
# Stages:  1 GRANAR inputs + anatomies   2 MECHA + statistics + figure   3 RootSlice configs + runs
#          4 costs, pairing, Arjun table  5 GRANAR drawings + figure icons (icons need photos)  6 checks (lignin cap, lateral endodermis)
#          7 Monte Carlo of kr (run + analyse)                          8 compare with reference ($COSTCOND_REF)
set -u
HERE="$(cd "$(dirname "$0")/.." && pwd)"; PY="${PYTHON:-python}"; export PYTHONWARNINGS=ignore
OUTD="${COSTCOND_OUT:-$HERE/outputs}"; mkdir -p "$OUTD"; LOG="$OUTD/run_all_log.txt"
exec > >(tee -a "$LOG") 2>&1
STAGES="${*:-1 2 3 4 5 6 8}"
echo "######## run_all started $(date); stages: $STAGES ########"; git -C "$HERE" log -1 --format='commit %h %cd' 2>/dev/null
run() { echo "---- $*"; "$@" || { echo "FAILED: $*"; exit 1; }; }
for s in $STAGES; do case $s in
  1) run $PY "$HERE/1_granar/make_granar_inputs.py"; run Rscript "$HERE/1_granar/build_anatomies.R" ;;
  2) run $PY "$HERE/2_mecha/batch.py"; run $PY "$HERE/2_mecha/analyse.py" ;;
  3) run $PY "$HERE/3_rootslice/make_configs.py"; run bash "$HERE/3_rootslice/run_rootslice.sh"
     run $PY "$HERE/3_rootslice/make_configs_smallstele.py"; run bash "$HERE/3_rootslice/run_small_stele.sh" ;;
  4) run $PY "$HERE/4_costs/rep_costs.py"; run $PY "$HERE/4_costs/rep_stats.py"; run $PY "$HERE/4_costs/fig_rep.py"
     if ls "$OUTD/mecha/sensitivity_mc/mc_rep_x2.csv" >/dev/null 2>&1; then run $PY "$HERE/2_mecha/sensitivity_mc/analyse_mc.py"; run $PY "$HERE/2_mecha/sensitivity_mc/fig_mc.py"
     else echo "   (no Monte Carlo results: the Arjun table is built without the kr uncertainty sheet; run stage 7 or copy them from Box)"; fi
     run $PY "$HERE/4_costs/build_xlsx.py" ;;
  5) run $PY "$HERE/figures/plot_granar.py"; run $PY "$HERE/figures/make_icons.py"; run $PY "$HERE/figures/fig_icons.py" ;;
  6) run $PY "$HERE/2_mecha/checks/lignin_cap.py"; run $PY "$HERE/2_mecha/checks/bracket_laterals.py" ;;
  7) run bash "$HERE/2_mecha/sensitivity_mc/run_mc.sh" ;;
  8) run $PY "$HERE/reproduce/compare_to_reference.py" ;;
  *) echo "unknown stage $s" ;;
esac; done
echo "######## run_all finished $(date) ########"
