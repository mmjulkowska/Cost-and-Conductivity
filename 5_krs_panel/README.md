# Step 5: root-system conductance (Krs) of the traced *S. pimpinellifolium* panel

Uses the step-2 segment conductivities to compute whole-root-system hydraulic conductance (Krs) for every SmartRoot-traced
seedling of PimpiBig1–6 (control and 100 mM NaCl, tracing days 1–5), then:
- summarizes Krs over time and under salt, and estimates its heritability and genotype × salt interaction;
- projects each plant's architecture to day 10 from its own growth rates and computes Krs each day;
- prepares and runs an EMMAX GWAS on the server (no significant loci; kept as a negative result);
- relates the salt response to genetic group and collection site.

## How Krs is computed (`krs_tree.py`)
- **Network from the tracing:** main root of its traced length; each first-order lateral inserted at its traced position with
  its traced length; 0.1 cm segments. With uniform soil water potential and no gravity, Krs depends only on topology and
  segment lengths, so the traced plant is used directly; no CPlantBox growth simulation is needed.
- **Segment classes** by distance from each root's own tip, boundaries at the midpoints between the sampled positions
  (tip, 3, 6, 9 cm): main root MT0 (0–1.5 cm, no barrier), MT (1.5–4.5, B1), MM (4.5–7.5, B6), MP (> 7.5, B6);
  laterals LT (0–1.5, no barrier), LM (1.5–4.5, B10), LP (> 4.5, B10).
- **Properties:** per-class medians of kr and Kx from `outputs/mecha/conductivities_per_anatomy.csv`; radius from Emily Miller's
  measured cross-section area (`outputs/mecha/anatomy_by_segment.csv`, MT0 = MT). Variant `lat_endo`: LM/LP kr × the B6/B10
  ratio from `outputs/mecha/checks/bracket_laterals_summary.csv`.
- **Solver:** Doussan-type network; radial conductance 2πa·l·kr split between a segment's two nodes, axial Kx/l; collar at the
  main-root base. Krs (cm³ hPa⁻¹ d⁻¹) = uptake / (soil − collar potential). Matches the analytic single-root solution
  and CPlantBox-based Krs on the same plants to within 0.15%.

## Inputs (not in git; put them in `data/krs/` or point `KRS_DATA` at their folder)

| File | What | Source |
|---|---|---|
| `pooled_plants_by_day.csv`, `pooled_laterals.csv` | one row per plant and day / per lateral (≥ 0.1 cm) and day, PimpiBig1–6 | `prep/pool_exports.py` from the SmartRoot exports (`Tomato_Pareto/ZENODO_deposit/data/06_smartroot_root_exports`) |
| `timeseries_roots.csv` | PimpiBig2 days 1–5 from the RSML | `prep/parse_timeseries.py` from Dropbox `Tomato_Pimpinellifolium/Pimp_RSA/SmartRooted/pimpiBIG2` |
| `01_panel_root_geometry.csv` | panel table: accession code, TGRC code, condition, architectural type | Pareto-front panel (CPlantBox_RSA `input/`) |
| `sample_key.csv`, `ALLCHR.filt.prune.tfam` | genotype sample order (265) | GWAS bundle `Pimp_Pareto_GWAS_for_Ruth/samples/` |
| `S.pimpinellifolium GPS data.xlsx` | latitude/longitude of 227 accessions | Dropbox `Tomato_Pimpinellifolium/Pimp_SNPs/` |
| `pcs_with_krs_traits.tsv` | PCs and k-means genetic group per accession (optional, for `krs_geo.py`) | written on the server by `gwas/krs_vs_structure.R` |

Also needed: the step-2 outputs in `outputs/mecha/` (from Box, see `data/README.md`).
Python: numpy, pandas, scipy, statsmodels, matplotlib, openpyxl; `basemap` for the map only.

## Run
```bash
bash 5_krs_panel/run_all.sh          # ~25 min; NBOOT=20 bash 5_krs_panel/run_all.sh for a quick check
```
Outputs go to `outputs/krs/`.

| Script | Does | Main outputs |
|---|---|---|
| `krs_panel.py` | assembles the tracings (PimpiBig5 days remapped by image date: 20 Oct = d2 … 23 Oct = d5; 19 Oct dropped; PimpiBig2 days 1–3 from the RSML), Krs for every plant-day, plus tip-zone (0.5 / 3 cm) and lateral-endodermis variants | `krs_per_plant_day.csv` (5,037 plant-days) |
| `krs_heritability.py` | day 5: H² from `y ~ experiment + (1|accession)` per condition, 300 accession bootstraps; G × E model; experiment-adjusted accession values | `krs_heritability.csv`, `krs_gxe.csv`, `krs_accession_values.csv` |
| `krs_project.py` | per-plant growth parameters and projection to day 10 with two lateral rules (accelerating, plant-scaled; constant per-plant rate); Krs and H² per day | `growth_params_fit5.csv`, `krs_projected_fit5.csv`, `krs_projected_H2.csv` |
| `hindcast.py` | fit on days ≤ 3, predict days 4–5 (760 PimpiBig2/5/6 plants) | `hindcast_kinetic.csv` |
| `gwas/krs_gwas_prep.py`, `gwas/krs_gwas_expcov.py` | EMMAX phenotype files (265 rows, `.tfam` order, Blom INT) and the run scripts | `krs_gwas/`, `krs_gwas_expcov/` |
| `krs_geo.py` | salt response vs latitude, longitude, Ecuador/Peru and genetic group | `krs_accession_geo.csv` |
| `plot_krs_panel.py`, `plot_projection.py`, `plot_krs_map.py` | figures | `krs_panel_summary`, `krs_projection_summary`, `krs_salt_map` (.png/.pdf) |

Accession link for the GWAS: tracing code → TGRC code (panel table) → `.tfam` IID. Do not use `sample_key.genotype_code`: it
covers only the 199 accessions of the Pareto-trait file and drops 44 genotyped accessions (all PimpiBig5).

## GWAS on the server (`gwas/`)
Reuses the October 2026 EMMAX pipeline (`v2_filter_markers.sh` … `v5_check_calibration.R`, in `~/pareto_gwas/scripts`)
with the September genotype bundle (`~/pimp_gwas_for_ruth`).
```bash
scp -r outputs/krs/krs_gwas vash:~/pareto_gwas_krs_v2                   # then on VASH, in tmux:
bash ~/pareto_gwas_krs_v2/run_krs_gwas.sh                                # MAF filter, 60 EMMAX jobs, annotation, calibration
scp -r outputs/krs/krs_gwas_expcov vash:~/pareto_gwas_krs_expcov        # experiment-covariate check (needs the run above)
bash ~/pareto_gwas_krs_expcov/run_expcov.sh
WORK=~/pareto_gwas_krs_v2 Rscript krs_vs_structure.R                     # PCs, genetic groups, pseudo-h2
```
v3 asks once after a test job: check that the sample count in the log matches `pheno/alignment_report.txt`, then type `y`.
The full results of the November 2026 run are archived as `krs_gwas_bundle_20261124` (see `gwas/README.md`).

## Results (November 2026)
- **Krs over time and salt.** Day-5 median Krs: control 8.5, salt 4.1 (×10⁻⁴ cm³ hPa⁻¹ d⁻¹); salt keeps ~0.48 of control
  from day 3 on. Mostly root size (ln Krs ~ ln total root length, R² = 0.90); salted plants also have fewer laterals,
  so fewer barrier-free tips. Control Krs at 30 cm total root length matches the CPlantBox architectures (7.8–10.2).
- **Heritability, day 5, accession means:** ln Krs 0.72 (control) / 0.76 (salt); Krs beyond size (residual of ln Krs on
  ln total root length) 0.73 / 0.72. Genotype × salt is large for Krs beyond size (accession means correlate r = 0.18
  between conditions).
- **Projection to day 10** (hindcast two days ahead: median predicted/traced Krs 0.91, log r 0.76–0.80 within treatment):
  median Krs 21.3 (control) and 10.5 (salt), accelerating rule; the two lateral rules rank plants identically
  (Spearman 0.998, values within ~6%) although total root length differs by ~25%. Krs grows ×2.5 while total root length
  grows ×5 (day 5 → 10): new length is mostly low-kr length behind the tips, and main-root xylem starts to limit.
  H² of projected Krs falls from 0.70 to 0.49 (control) and from 0.78 to 0.72 (salt) between days 6 and 10.
- **GWAS:** no Bonferroni or FDR 5% hits for any trait (best P 1.2e-7, n = 228–231). Salt ln Krs and the salt/control
  ratio are inflated (λ ≈ 1.28), unchanged with experiment covariates, so not batch.
- **Structure and geography:** the inflation follows genetic structure. Two minority genetic groups, containing 17 of the 19
  Ecuadorian accessions, lose the most Krs under salt; Ecuadorian accessions keep 37% of control Krs vs 47% for
  Peruvian ones (Mann–Whitney P = 0.002).

## Caveats
- Salt plants use the control segment conductivities: the salt effect on Krs is architectural only.
- H² may include shared-environment effects: most accessions were grown in one experiment only (experiment fitted as fixed);
  plates are not modelled.
- Krs depends on the barrier-free tip length (0.5 cm tip −31%, 3 cm tip +11–21%) and tip Kx assumes mature xylem.
- Projection: no second-order laterals, no maximal lateral length, no plate edges, constant main-root rate, unchanged
  conductivities with age; lateral emergence per plant from two to five traced days (5 cm spacing when a plant made no
  new laterals; group median when traced on one day).
- Ecuador/Peru split at 3.5°S is approximate; 182 of 240 accessions have coordinates.
