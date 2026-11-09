# 1. GRANAR: cross-section anatomies

GRANAR builds a 2D root cross-section (cells, walls, tissue labels) from a few measured sizes. MECHA uses these anatomies as its geometry.

## Input
`data/AllDatammAndpixelsNov2026_2.xlsx` holds Emily Miller's measurements in pixels (Box link in `data/README.md`): root, cortex and vasculature areas, the cortical layer count, and up to 4 measured cells per tissue per section.

Scale and size conversion:
- 311.26 px = 100 µm, i.e. 0.32127 µm px⁻¹ and **1.03214286 × 10⁻⁷ mm² px⁻¹**
- diameter d = 2√(A/π)

## `make_granar_inputs.py`
Writes one parameter file per section × measured cell k (k = 1–4):

| GRANAR parameter | Value |
|---|---|
| stele `layer_diameter` | 2√(A/π), A = vasculature total area |
| stele `cell_diameter` | vasculature cell #k (segment median of measured vasculature cells if #k is missing) |
| cortex `cell_diameter` | cortex cell #k |
| cortex `n_layers` | layer count − 1, plus 1 exodermis layer |
| endodermis, epidermis `cell_diameter` | endodermis / epidermis cell #k |
| pericycle | 0.021 mm |
| exodermis | 0.02 mm |
| xylem `max_size` | 0.035 mm (GRANAR's single central metaxylem) |
| phloem, root hair, pith | 0.01 mm; 0.07 mm; 0.03 / 0.01 mm |

**Excluded:** row 4 (DSC_00425PR, torn section), row 60 (DSC_0030RT, no stele visible), and any replicate missing the vasculature area, the layer count or cell #k.

**Outputs (tracked in git):**
- `v5_replicates.csv`: 268 rows; `used`, `reason`, and every size in mm. It is also the input for RootSlice.
- `granar_params/MP<n>_cell<k>.xml`: 216 parameter files.

## `build_anatomies.R`
- Seed = 10 × section + cell, plus the offset in `seed_overrides.txt`. MP8_cell3, MP10_cell2 and MP16_cell3 use +100, because their first anatomy made the MECHA solve singular.
- If GRANAR fails, the script retries with +100 … +1000 (11 seeds in total) and records the offset used.
- 8 replicates fail with every seed: MP19_cell4, MP37_cell2, MP57_cell1–3, MP59_cell1, MP66_cell1–2. They still get costs (RootSlice) but no conductivity.

**Outputs:**
- `outputs/granar/anatomies/MP<n>_cell<k>_root.xml` (208 anatomies)
- `outputs/granar/granar_anatomy_summary.csv` (seed offset, tissue areas and cell counts)

```
python 1_granar/make_granar_inputs.py
Rscript 1_granar/build_anatomies.R              # or: Rscript 1_granar/build_anatomies.R MP8_cell3
```
GRANAR's runtime field inside each anatomy file changes from run to run; the geometry does not.

## Look at them
`figures/plot_granar.py` draws each anatomy (`outputs/figures/granar/`) and one gallery per segment at a common scale (`outputs/figures/granar_gallery_<segment>.pdf`).

## Check against the photographs
`2_mecha/analyse.py` writes `outputs/mecha/anatomy_vs_emily.csv`. Median simulated / measured root area is 1.31; GRANAR packs layers more loosely than the real root.
