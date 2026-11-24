#!/bin/bash
# run_all.sh -- step 6, CPlantBox root architectures of the three topologies, their Krs and their water uptake in drying soil.
# Needs: the step-2 outputs (outputs/mecha/), the step-5 inputs in data/krs/ (made by 5_krs_panel/prep/) and the PimpiBig2
# day-5 RSML ($PIMPIBIG2_RSML). Outputs go to outputs/cplantbox/.
#   bash 6_cplantbox/run_all.sh              stages 1 3 4 5 with the fitted parameters stored in params/ (about 15 min)
#   FIT=1 bash 6_cplantbox/run_all.sh        refit every parameter set first (stage 2, about 1 h)
#   bash 6_cplantbox/run_all.sh 6            drought runs in pots (162 runs, about 1 h on 2 cores)
#   Stages: 1 targets from the tracings  2 fits (v1, v2, v3, v4)  3 simulations, galleries, checks  4 growth animations
#           5 Krs with the v5 conductivities  6 soil water and drought
set -euo pipefail
cd "$(dirname "$0")"; HERE="$(pwd)"
PY="${PYTHON:-python3}"; export PYTHONWARNINGS=ignore
OUTD="$($PY -c 'import cpb_paths as C; print(C.CPB_OUT)')"
STAGES="${*:-1 2 3 4 5}"
q() { "$@" 2>&1 | grep -v -e 'organTypeNumber' -e '^reference$' -e '^Organism::write' || true; }
for s in $STAGES; do case $s in
  1) echo '---- 1 targets'
     q $PY parse_rsml.py; q $PY kinetics.py; q $PY track_laterals.py; q $PY targets_v3.py ;;
  2) if [ "${FIT:-0}" = 1 ]; then echo '---- 2 fits'
       q $PY fit.py; q $PY fit_angle.py; q $PY fit_angle_plate.py; q $PY make_params_exaggerated.py
       for t in B C dT; do q $PY fit_v2.py $t 150; done; q $PY fit_angle_v2.py
       for t in B C dT; do q $PY fit_v3.py $t 150; done; q $PY fit_angle_v3.py
       q $PY make_params_v4.py; q $PY refit_v4.py
     else echo '---- 2 fitted parameters copied from params/ (FIT=1 to refit)'
       cp "$HERE"/params/*.json "$OUTD"/
     fi ;;
  3) echo '---- 3 simulations at 30 cm total root length, galleries, checks'
     q $PY simulate_equal_trl.py
     q env OUT=sims FIG=figs $PY plot_gallery.py; q env OUT=sims FIG=figs $PY plot_validation.py
     q env PARAMS=params_exaggerated.json OUT=sims_exaggerated $PY simulate_equal_trl.py
     q env OUT=sims_exaggerated FIG=figs_exaggerated TAG=' (exaggerated topologies)' $PY plot_gallery.py; q $PY plot_profiles.py
     q env PARAMS=params_v2.json OUT=sims_v2 $PY simulate_equal_trl.py
     q env OUT=sims_v2 FIG=figs_v2 TAG=' (refit with day 1–5 growth kinetics)' $PY plot_gallery.py
     q env OUT=sims_v2 FIG=figs_v2 $PY plot_validation.py; q $PY plot_kinetics.py ;;
  4) echo '---- 4 growth animations (v3, v4)'
     q env PARAMS=params_v3.json ANIM=anim $PY animate.py; q env PARAMS=params_v4.json ANIM=anim_v4 $PY animate.py ;;
  5) echo '---- 5 Krs'
     q $PY krs_analysis.py; q $PY plot_krs.py; q $PY animate_krs.py ;;
  6) echo '---- 6 soil water and drought'
     q $PY soil_batch.py; q $PY plot_soil.py; q $PY animate_soil.py loam drydown ;;
  *) echo "unknown stage $s" ;;
esac; done
echo "done: outputs in $OUTD"
