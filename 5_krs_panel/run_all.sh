#!/bin/bash
# run_all.sh -- step 5, Krs of the traced panel, from inputs in data/krs/ to outputs/krs/.
# Needs the step-2 outputs (outputs/mecha/) and the inputs listed in README.md.
# NBOOT=20 bash run_all.sh for a quick check (the bootstrap of the heritability CIs is the slow part, ~20 min at 300).
set -euo pipefail
cd "$(dirname "$0")"
python3 krs_panel.py            # Krs for every plant and tracing day
python3 krs_heritability.py     # H2 (day 5), G x E, accession values
python3 plot_krs_panel.py
python3 krs_project.py          # projection to day 10, both lateral rules
python3 hindcast.py             # fit days <= 3, predict days 4-5
python3 plot_projection.py
python3 gwas/krs_gwas_prep.py   # EMMAX phenotype files -> outputs/krs/krs_gwas/
python3 gwas/krs_gwas_expcov.py # experiment-covariate check -> outputs/krs/krs_gwas_expcov/
python3 krs_geo.py              # geography (needs pcs_with_krs_traits.tsv from the server for the groups)
python3 plot_krs_map.py
echo "done: outputs in $(cd ../outputs/krs 2>/dev/null && pwd || echo outputs/krs)"
