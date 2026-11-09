# 3. RootSlice: cortical volumes

RootSlice builds a 3D slice of root cortex from cell sizes and layer counts. It reports cortical cell-wall, cytoplasm and vacuole volumes, which `4_costs/` converts into carbon and nitrogen costs.

## Software
- RootSlice, Jagdeep Sidhu's Apple Silicon build (`RootSlice_for_Ruth.zip`, Box link in `data/README.md`). It is third-party software and is not distributed here.
- Volumes are collected by the build's own `analysis/RootSlice_volume_analysis.R`.
- Install: `bash env/setup_rootslice_mac.sh /path/to/RootSlice_for_Ruth.zip` (installs to `~/root_models/RootSlice`, or `$ROOTSLICE_HOME`).

## `make_configs.py`
Reads `1_granar/v5_replicates.csv`, so every RootSlice replicate uses the same measured sizes as its GRANAR anatomy. It writes `configs/R5_MP<n>_cell<k>.xml` (216) and `configs_values.csv`.

Template `template_rootslice.xml`: no aerenchyma (RCA 0), cell length 100 µm, and Jagdeep's settings for xylem, exodermis and sclerenchyma. Per replicate:

| Field | Value |
|---|---|
| Cortical `BaseRadius` | stele radius (vasculature diameter / 2) |
| Cortical `CellDiameter` | cortex cell #k |
| Cortical `NumFiles` | layer count (exodermis counted as a layer) |
| Cortical `Thickness` | cell diameter × layers |
| Epidermis, Endodermis `CellDiameter` | cell #k |
| Stele `InnestCellDiameer` | vasculature cell #k (segment median if missing) |

## Runs
- `run_rootslice.sh` runs all 216 configs, with a 180 s limit and one retry per run.
  - Raw runs: `$ROOTSLICE_HOME/results_costcond_v5/` (archive on Box)
  - Volumes: `outputs/rootslice/volume_analysis_configs.csv` (197 finished)
- 19 replicates with small steles (radius 19–31 µm) hang or crash. `make_configs_smallstele.py` writes simplified stele interiors to `configs_smallstele/`, and `run_small_stele.sh` tries them in order (120 s each):
  - **A:** fewer xylem files.
  - **B:** A plus the inner stele cell × 0.7.
  - **C:** one metaxylem.
- The stele interior does not enter the cortical volumes: four controls (`CTRL_`) that ran with both templates gave identical cortical volumes.
- 13 replicates were recovered this way and 6 failed every variant (MP34_cell2, MP57_cell1–3, MP59_cell1, MP66_cell2). Records: `outputs/rootslice/smallstele_variant_used.txt` and `volume_analysis_smallstele.csv`.

```
python 3_rootslice/make_configs.py && bash 3_rootslice/run_rootslice.sh          # ~1.5 h
python 3_rootslice/make_configs_smallstele.py && bash 3_rootslice/run_small_stele.sh
```
