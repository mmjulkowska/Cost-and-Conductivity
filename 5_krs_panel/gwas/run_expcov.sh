#!/bin/bash
# run_expcov.sh -- EMMAX with experiment covariates, for the inflated salt Krs traits (+ control ln Krs as reference).
# Reuses the MAF-filtered markers of ~/pareto_gwas_krs_v2 (same 232-accession analysis set).
#   From outputs/krs (made by krs_gwas_expcov.py):  scp -r krs_gwas_expcov vash:~/pareto_gwas_krs_expcov
#   On VASH, in tmux:                         bash ~/pareto_gwas_krs_expcov/run_expcov.sh
set -euo pipefail
export PATH=$HOME/anaconda3/bin:$PATH
export WORK=$(cd "$(dirname "$0")" && pwd)
V2=$HOME/pareto_gwas_krs_v2
G=$HOME/pimp_gwas_for_ruth/genotypes
SCRIPTS=${SCRIPTS:-$HOME/pareto_gwas/scripts}
EMMAX=$V2/bin/emmax-intel64
KIN=$G/ALLCHR.filt.prune.aBN.kinf
ln -sfn "$V2/tped" "$WORK/tped"
# covariates = intercept + 5 PCs (as before) + 5 experiment-share columns
paste "$G/ALLCHR.filt.prune.PCA.txt" "$WORK/exp_covariates.txt" > "$WORK/cov_pca_exp.txt"
awk '{print NF}' "$WORK/cov_pca_exp.txt" | sort -u | sed 's/^/covariate columns: /'
wc -l < "$WORK/cov_pca_exp.txt" | sed 's/^/covariate rows (must be 265): /'
mkdir -p "$WORK/out" "$WORK/logs"
for p in "$WORK"/pheno/emmax/trait.*.INT.pheno; do
  for C in 01 02 03 04 05 06 07 08 09 10 11 12; do printf "%s\t%s\n" "$p" "$C"; done
done | xargs -P 6 -n 2 bash -c '
  t=$(basename "$0" .INT.pheno); t=${t#trait.}
  '"$EMMAX"' -v -d 10 -t '"$WORK"'/tped/CHR$1 -p "$0" -k '"$KIN"' -c '"$WORK"'/cov_pca_exp.txt \
     -o '"$WORK"'/out/$t.CHR$1 > '"$WORK"'/logs/$t.CHR$1.log 2>&1'
echo ".ps: $(ls "$WORK"/out/*.ps | wc -l) of 36   .reml: $(ls "$WORK"/out/*.reml | wc -l) of 36"
bash    "$SCRIPTS/v4_annotate.sh"
Rscript "$SCRIPTS/v5_check_calibration.R"
cat "$WORK/calibration_summary.tsv"
