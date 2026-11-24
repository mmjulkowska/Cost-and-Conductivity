#!/bin/bash
# run_krs_gwas.sh -- EMMAX GWAS of the Krs traits on VASH, reusing the Oct 2026 pipeline unchanged.
#
# Made by 5_krs_panel/gwas/krs_gwas_prep.py in outputs/krs/krs_gwas/. From that outputs/krs folder:
#   scp -r krs_gwas vash:~/pareto_gwas_krs_v2
# Then on VASH:
#   tmux new -s krs
#   bash ~/pareto_gwas_krs_v2/run_krs_gwas.sh
#
# Own WORK folder, so nothing in ~/pareto_gwas (the gravitropism traits) is touched.
# v2 is re-run because the MAF filter is computed on the phenotyped subset, and the Krs
# analysis set (232 accessions) is not the same as the Pareto-trait set.
# v3 times one job and asks before running the rest.
set -euo pipefail
export PATH=$HOME/anaconda3/bin:$PATH
export WORK=$(cd "$(dirname "$0")" && pwd)   # the folder this script sits in
export BUNDLE=${BUNDLE:-$HOME/pimp_gwas_for_ruth}
SCRIPTS=${SCRIPTS:-$HOME/pareto_gwas/scripts}
mkdir -p "$WORK/bin"
[[ -x "$WORK/bin/emmax-intel64" || ! -x "$HOME/pareto_gwas/bin/emmax-intel64" ]] || ln -sf "$HOME/pareto_gwas/bin/emmax-intel64" "$WORK/bin/emmax-intel64"   # else v3 looks for emmax on PATH
cat "$WORK/pheno/alignment_report.txt"
bash    "$SCRIPTS/v2_filter_markers.sh"
bash    "$SCRIPTS/v3_run_emmax.sh"
bash    "$SCRIPTS/v4_annotate.sh"
Rscript "$SCRIPTS/v5_check_calibration.R"
echo "Results: $WORK/tables, $WORK/plots (check every trait reads 'ok' before interpreting)"
