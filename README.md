# Cost and Conductivity: tomato root segment anatomy, hydraulics and construction cost

Code and model inputs used to compute the radial conductivity, axial conductance and carbon and nitrogen cost of root segments of *Solanum pimpinellifolium*, for a Pareto-front analysis of root system architecture.
Julkowska Lab, Boyce Thompson Institute.

Each measured cell in each root cross-section photographed by Emily Miller is one **replicate**. Each replicate goes through three models, all built from the same measured sizes (steps 1–4). Step 5 then uses the segment conductivities to compute whole-root-system conductance (Krs) for every traced seedling of the panel, and step 6 builds the three root-system topologies of the panel in CPlantBox at equal size to ask how much topology itself matters for Krs and for water uptake in drying soil:

| Step | Model | What it gives | Folder |
|---|---|---|---|
| 1 | **GRANAR** (R) | 2D cross-section anatomy | [`1_granar/`](1_granar/README.md) |
| 2 | **MECHA** (Python) | radial conductivity kr and axial conductance Kx, with the barrier set for each segment | [`2_mecha/`](2_mecha/README.md) |
| 3 | **RootSlice** | cortical cell-wall, cytoplasm and vacuole volumes | [`3_rootslice/`](3_rootslice/README.md) |
| 4 | Python | carbon and nitrogen costs, pairing with conductivity, statistics, table for the Pareto model | [`4_costs/`](4_costs/README.md) |
| 5 | Python | root-system conductance (Krs) of every traced seedling and day, projection to day 10, heritability, GWAS inputs and server scripts, geography of the salt response | [`5_krs_panel/`](5_krs_panel/README.md) |
| 6 | **CPlantBox** (Python) | Christmas tree, Broomstick and Droopy telephone pole architectures fitted to the traced panel (shape and day 1–5 growth), compared at equal total root length: galleries and growth animations, Krs with the step-2 conductivities, water uptake from pots of three soils with and without drought | [`6_cplantbox/`](6_cplantbox/README.md) |
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
                                   outputs/deliverables/20261109_Root_costs_conductivity_for_Arjun_v5.xlsx
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

## Step 5 in numbers

| What | n |
|---|---|
| Traced plant-days (PimpiBig1–6, control and 100 mM NaCl, tracing days 1–5) | 5,037 |
| Plants on day 5 / accessions | 1,562 / 240 |
| Accessions traced and genotyped (GWAS) | 232 |
| Plants projected to day 10 (two lateral-growth rules) / hindcast check | 1,573 / 760 |

Day-5 median Krs is 8.5 (control) and 4.1 (salt) × 10⁻⁴ cm³ hPa⁻¹ d⁻¹, with accession-mean H² 0.72–0.76. The GWAS found no significant locus; the salt-response traits follow genetic structure and geography (Ecuadorian accessions lose more Krs under salt). Details in [`5_krs_panel/README.md`](5_krs_panel/README.md).

## Step 6 in numbers

| What | n |
|---|---|
| Topologies × replicates per parameter set, grown to 30 cm total root length (2D plate and 3D) | 3 × 12 |
| Control plants behind the fitted shape / growth time courses (PimpiBig1–6 / Big2, 5, 6) | 121–197 / 71–129 per type |
| Pot simulations (3 topologies × 3 soils × 3 water scenarios × fixed or growing roots × 3 replicates, 10 days each) | 162 |

At 30 cm total root length, Krs is 10.2 (Droopy telephone pole), 8.2 (Christmas tree) and 7.8 (Broomstick) × 10⁻⁴ cm³ hPa⁻¹ d⁻¹; 78–88% of the water enters through barrier-free lateral tips. In every soil and drought scenario the order of 10-day transpiration is the same (Droopy telephone pole > Christmas tree > Broomstick, 5–10 points of demand apart); continued root growth matters more (+15–25 points). Details in [`6_cplantbox/README.md`](6_cplantbox/README.md).

## What is where

- **In this repository:** all code; the model inputs derived from the sheet (GRANAR parameter files, RootSlice configs, `v5_replicates.csv`); READMEs.
- **Step 5 inputs** (tracing tables, genotype sample key, collection-site coordinates) go in `data/krs/` (git-ignored); sources are listed in [`5_krs_panel/README.md`](5_krs_panel/README.md). The archived GWAS results are described in [`5_krs_panel/gwas/README.md`](5_krs_panel/gwas/README.md).
- **Step 6** reads the step-5 inputs in `data/krs/` plus the PimpiBig2 day-5 SmartRoot RSML (`PIMPIBIG2_RSML`); its fitted parameters are kept in `6_cplantbox/params/`, its outputs go to `outputs/cplantbox/`.
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
| Step 5 | Python with numpy, scipy, pandas, statsmodels, matplotlib, openpyxl; `basemap` for the map. GWAS on the server: EMMAX (`emmax-intel64`), PLINK 1.90b6.21, R |
| Step 6 | CPlantBox 2.2.0 (`pip install CPlantBox`, Python 3.13), numpy, scipy, pandas, matplotlib; ffmpeg for the animations |

Install: `bash env/install_mecha_mac.sh` and `bash env/setup_rootslice_mac.sh /path/to/RootSlice_for_Ruth.zip`.

## Rerun and check

```
conda activate mecha_dev
bash reproduce/run_all.sh                  # stages 1–6 and 8, about 2 h (RootSlice is most of it)
bash 5_krs_panel/run_all.sh                # step 5, about 25 min (NBOOT=20 for a quick check); needs outputs/mecha/ and data/krs/
bash 6_cplantbox/run_all.sh               # step 6, about 15 min with the stored parameters (FIT=1 to refit, stage 6 for the drought runs, ~1 h)
```
Stage 8 compares every output table with the frozen reference from Box (`COSTCOND_REF=/path/to/reference`). A rerun from this repository reproduced all 41 reference tables to within 2 × 10⁻¹⁰ relative difference (9 Nov 2026).

## Known limitations

- **Barriers come from toluidine-blue brightfield images,** which cannot resolve a Casparian strip. In the laterals (LP, LM) the endodermal barrier is "not detected", not "absent". `2_mecha/checks/` shows that adding one halves lateral kr.
- **Xylem is GRANAR's single central metaxylem,** but tomato is diarch. Kx depends on this choice.
- **The true-tip conductivity (MT0) uses the MT anatomy with mature xylem,** so its Kx is an upper bound.
- **Step 5 uses the control anatomy for salt-grown plants,** so its salt effects on Krs are architectural only; see [`5_krs_panel/README.md`](5_krs_panel/README.md) for the projection and GWAS caveats.
- **Step 6 extrapolates beyond tracing day 5** (day-12 architectures, growing pot runs) without second-order laterals, and its soil model is 1D; see [`6_cplantbox/README.md`](6_cplantbox/README.md).

## Citation and licence

Please cite this repository (see [`CITATION.cff`](CITATION.cff)) and the models it uses: GRANAR (Heymans et al. 2020, *Plant Physiology*), MECHA (Couvreur et al. 2018, *Plant Physiology*) and RootSlice (Sidhu et al. 2023); for step 5, the root hydraulic network of Doussan et al. (1998, *Annals of Botany*) and EMMAX (Kang et al. 2010, *Nature Genetics*); for step 6, CPlantBox (Schnepf et al. 2018, *Annals of Botany*; Zhou et al. 2020, *in silico Plants*).
- **Code:** MIT licence ([`LICENSE`](LICENSE)).
- **Data, model inputs, outputs, figures and documentation:** CC BY 4.0 ([`LICENSE-data.md`](LICENSE-data.md)). Reuse requires attribution.
