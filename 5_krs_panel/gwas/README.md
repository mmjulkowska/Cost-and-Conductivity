# Krs GWAS — S. pimpinellifolium panel (negative result), Nov 2026

This README describes the archived results (`krs_gwas_bundle_20261124`, 902 MB; Box / lab Zenodo deposit). The scripts in
this folder regenerate the inputs (`krs_gwas_prep.py`, `krs_gwas_expcov.py`) and run the analysis on the server
(`run_krs_gwas.sh`, `run_expcov.sh`, `krs_vs_structure.R`); see `../README.md`.

GWAS of simulated root-system hydraulic conductance (Krs) traits, computed from SmartRoot tracings of PimpiBig1–6
(day 5, control and 100 mM NaCl) with the v5 MECHA segment conductivities. See `../README.md` for how Krs was computed.

## Result
No association reached Bonferroni (0.05 / 4.77 M markers ≈ 1.0e-8) or FDR 5% in any trait or model.

| Run | Trait | n | λ | Best P |
|---|---|---|---|---|
| main_232acc | Krs beyond size, control (KrsResid.C) | 231 | 0.97 | 1.2e-7 |
| main_232acc | Krs beyond size, salt (KrsResid.S) | 228 | 0.99 | 2.4e-7 |
| main_232acc | ln Krs, control (logKrs.C) | 231 | 1.03 | 3.7e-6 |
| main_232acc | ln Krs, salt (logKrs.S) | 228 | 1.28 (inflated) | 1.2e-6 |
| main_232acc | ln Krs salt − control (logKrs.SminusC) | 227 | 1.28 (inflated) | 4.9e-6 |
| expcov_check | ln Krs, control, raw means + experiment covariates | 231 | 1.05 | 2.4e-6 |
| expcov_check | ln Krs, salt, raw means + experiment covariates | 228 | 1.29 (inflated) | 2.2e-6 |
| expcov_check | ln Krs salt − control, raw means + experiment covariates | 227 | 1.26 (inflated) | 1.2e-6 |

The inflation of the salt traits is not caused by experiment batches (unchanged with experiment covariates). It follows
genetic structure: two minority genetic groups, which contain 17 of the 19 Ecuadorian accessions, lose more Krs under salt
(`main_232acc/structure/`; `outputs/krs/krs_accession_geo.csv`). Salt-trait P-values should not be interpreted.

## Methods
- Markers: 265 accessions, 17,962,579 SNPs (LA2093 reference), MAF ≥ 0.05 and missingness ≤ 0.10 computed on the
  232 phenotyped accessions → 4,769,914 markers.
- EMMAX (`emmax-intel64`), IBS kinship `ALLCHR.filt.prune.aBN.kinf`, covariates = intercept + 5 PCs
  (expcov_check: + share of each accession's plants in PimpiBig2–6).
- Phenotypes: experiment-adjusted accession means (main) or raw accession means (expcov_check) of day-5 plants,
  rank-based inverse normal transformed (Blom). Accession link: tracing code → TGRC code (panel geometry table) → `.tfam` IID.

## Contents
Each run folder:
- `summary_stats/<trait>.tsv.gz`: all tested markers, columns CHR, BP, SNP, STAT (β), SE, P
- `tables/`: Bonferroni and FDR hit tables (empty) and the 20 best markers per trait
- `plots/`: Manhattan and QQ plots
- `calibration_summary.tsv`: λ, verdict, tail excess, hit counts, best P
- `reml/`: EMMAX variance components; h² = line 4 / (line 4 + line 5)
- `pheno/`: the exact phenotype files analyzed; `logs.tgz`: per-job EMMAX logs
- `main_232acc/structure/`: PCs joined to the traits, `krs_vs_structure.R`
- `expcov_check/cov_pca_exp.txt`: the covariate file used

`scripts/`: the marker filter, EMMAX run, annotation and calibration scripts (v2–v5) and the two run scripts.
