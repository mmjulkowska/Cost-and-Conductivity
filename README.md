# Cost and Conductivity: tomato root segment anatomy, hydraulics and construction cost

Code and model inputs used to compute the radial conductivity, axial conductance and carbon and nitrogen cost of root segments of *Solanum pimpinellifolium*, for a Pareto-front analysis of root system architecture.
Julkowska Lab, Boyce Thompson Institute.

Each measured cell in each root cross-section photographed by Emily Miller is one **replicate**. Each replicate goes through three models, all built from the same measured sizes:

| Step | Model | What it gives | Folder |
|---|---|---|---|
| 1 | **GRANAR** (R) | 2D cross-section anatomy | [`1_granar/`](1_granar/README.md) |
| 2 | **MECHA** (Python) | radial conductivity kr and axial conductance Kx, with the barrier set for each segment | [`2_mecha/`](2_mecha/README.md) |
| 3 | **RootSlice** | cortical cell-wall, cytoplasm and vacuole volumes | [`3_rootslice/`](3_rootslice/README.md) |
| 4 | Python | carbon and nitrogen costs, pairing with conductivity, statistics, table for the Pareto model | [`4_costs/`](4_costs/README.md) |
| – | Python | drawings of every GRANAR anatomy, segment galleries, and the figure with half-photo / half-GRANAR cross-sections | [`figures/`](figures/README.md) |

```
Box: Emily Miller's measurement sheet (pixels)
        │  1_granar/make_granar_inputs.py   (1.03214286e-7 mm² px⁻¹; d = 2√(A/π))
        ▼
1_granar/v5_replicates.csv (268 rows, 216 used) ─────────────┐
        │                                                    │ 3_rootslice/make_configs.py
1_granar/granar_params/ (216)                                ▼
        │ build_anatomies.R                    3_rootslice/configs/ (216) + configs_smallstele/ (61)
        ▼                                                    │ RootSlice
outputs/granar/anatomies/ (208)                              ▼
        │ 2_mecha/batch.py                     outputs/rootslice/volume_analysis_*.csv (197 + 13)
        ▼                                                    │
outputs/mecha/conductivities_per_anatomy.csv (228) ──────────┤ 4_costs/rep_costs.py
                                                             ▼
                                   outputs/costs/replicate_level_costs_conductivity.csv (207 paired)
                                                             │ 4_costs/build_xlsx.py
                                                             ▼
                                   outputs/deliverables/Root_costs_conductivity_for_Arjun_v5_Nov2026.xlsx
```

## Numbers at each step

| Step | n |
|---|---|
| Section rows in the sheet × 4 measured cells | 268 (67 sections) |
| Replicates used (row 4 torn, row 60 no stele, and unmeasured cells excluded) | 216 (55 sections) |
| GRANAR anatomies (8 failed to tessellate with 11 seeds) | 208 |
| MECHA rows (208 anatomies + 20 main-tip anatomies rerun as the true tip, MT0) | 228 |
| RootSlice replicates (197 full template + 13 with a simplified stele interior; 6 failed) | 210 |
| Replicates with both cost and conductivity | 207 (53 sections) |

Segments (sheet row → segment): MP 1–9, MM 10–20, MT 21–26, LP 27–40, LM 41–50, LT 51–67. Positions are 3 cm apart along the main root (MP, MM, MT, tip) and a lateral root (LP, LM, LT).

## What is where

- **In this repository:** all code; the model inputs derived from the sheet (GRANAR parameter files, RootSlice configs, `v5_replicates.csv`); READMEs.
- **On Box:** raw data, model outputs, raw RootSlice runs, the frozen reference and the deliverable table. Links are in [`data/README.md`](data/README.md).
- **Paths:** set once in [`paths.py`](paths.py). Downloaded data goes in `data/`, outputs are written to `outputs/`. Both can be moved with the environment variables `COSTCOND_DATA`, `COSTCOND_OUT` and `ROOTSLICE_HOME`.

## Software

| Tool | Version |
|---|---|
| GRANAR | R package `granar` 0.0.1.2 |
| MECHA | `openalea.mecha` from [MECHA-dev](https://github.com/openalea-incubator/MECHA-dev) `develop` @ 945b2be (26 Oct 2026), with a one-line rename fix applied by the installer |
| GRANAP | GRANAP-dev @ 248adca |
| Python | 3.12, conda env `mecha_dev` (numpy, scipy, pandas, statsmodels, openpyxl, matplotlib, scikit-image, pillow) |
| RootSlice | Jagdeep Sidhu's Apple Silicon build (`RootSlice_for_Ruth.zip`, on Box) with its `RootSlice_volume_analysis.R` |

Install: `bash env/install_mecha_mac.sh` and `bash env/setup_rootslice_mac.sh /path/to/RootSlice_for_Ruth.zip`.

## Rerun and check

```
conda activate mecha_dev
bash reproduce/run_all.sh                  # stages 1–6 and 8, about 2 h (RootSlice is most of it)
```
Stage 8 compares every output table with the frozen reference from Box (`COSTCOND_REF=/path/to/reference`). A rerun from this repository reproduced all 41 reference tables to within 2 × 10⁻¹⁰ relative difference (9 Nov 2026).

## Known limitations

- **Barriers come from toluidine-blue brightfield images,** which cannot resolve a Casparian strip. In the laterals (LP, LM) the endodermal barrier is "not detected", not "absent". `2_mecha/checks/` shows that adding one halves lateral kr.
- **Xylem is GRANAR's single central metaxylem,** but tomato is diarch. Kx depends on this choice.
- **The true-tip conductivity (MT0) uses the MT anatomy with mature xylem,** so its Kx is an upper bound.

## Citation and licence

Please cite this repository (see [`CITATION.cff`](CITATION.cff)) and the models it uses: GRANAR (Heymans et al. 2020, *Plant Physiology*), MECHA (Couvreur et al. 2018, *Plant Physiology*) and RootSlice (Sidhu et al. 2023).
- **Code:** MIT licence ([`LICENSE`](LICENSE)).
- **Data, model inputs, outputs, figures and documentation:** CC BY 4.0 ([`LICENSE-data.md`](LICENSE-data.md)). Reuse requires attribution.
