# Data on Box

Raw data and all outputs live on Cornell Box. Nothing in this table is stored in git.

## Inputs to download

Download these and place them in `data/`, or point `COSTCOND_DATA` at a folder that holds them.

| What | Put it at | Box location | Public link |
|---|---|---|---|
| Emily Miller's measurement sheet (areas in pixels; one row per section image) | `data/AllDatammAndpixelsNov2026_2.xlsx` | Plant.Architecture.Lab › Emily Miller › Cross Sectioning - Main and Lateral Root › [Data](https://cornell.app.box.com/folder/337849184727) › November2026 x 100 › AllDatammAndpixelsNov2026_2.xlsx | https://cornell.box.com/s/92t6m9hzx2lnzyq2wgohxmd3frrr28wu |
| Emily Miller's cross-section photographs (toluidine blue, brightfield). Only needed for the figure icons | `data/Cross Sections/` (keep the subfolders) | Plant.Architecture.Lab › Emily Miller › Cross Sectioning - Main and Lateral Root › [Cross Sections](https://cornell.app.box.com/folder/341762726334) | https://cornell.box.com/s/c3v406j5xdznxreyrprmz7zgxxh2fvfp |
| RootSlice build `RootSlice_for_Ruth.zip` (Jagdeep Sidhu; third-party, not public) | anywhere; pass to `env/setup_rootslice_mac.sh` | Plant.Architecture.Lab › Magda › NSF Math Bio - ParetoFront › 20261109_GRANAR_MECHA_RootSlice › inputs (zip as sent by Jagdeep Sidhu) | lab access only |

Scale used for the sheet: 311.26 px = 100 µm (stage micrometer), i.e. 1.03214286 × 10⁻⁷ mm² px⁻¹.

## Outputs and archives

Box project folder: Plant.Architecture.Lab › Magda › NSF Math Bio - ParetoFront › 20261109_GRANAR_MECHA_RootSlice, with `inputs/` (third-party software) and `outputs/` (everything the pipeline produces).

| What | Box location |
|---|---|
| Model outputs (contents of `outputs/`: anatomies, conductivities, statistics, costs, figures, Monte Carlo) | Plant.Architecture.Lab › Magda › NSF Math Bio - ParetoFront › 20261109_GRANAR_MECHA_RootSlice › outputs › `Cost-and-Conductivity_outputs_20261109.zip` |
| Raw RootSlice run folders of the reference run: folder `RootSlice_raw_runs_v5_20261109/` with one .tgz per run (`results_configs_v5rep_v5/`, 216; `results_v5rep_smallstele/`, 17) plus the collected volume_analysis.csv; a rerun from this repository writes `results_costcond_v5*` | Plant.Architecture.Lab › Magda › NSF Math Bio - ParetoFront › 20261109_GRANAR_MECHA_RootSlice › outputs › `RootSlice_raw_runs_v5_20261109/` |
| Deliverable sent to the Pareto model (`Root_costs_conductivity_for_Arjun_v5_Nov2026.xlsx`) | Plant.Architecture.Lab › Magda › NSF Math Bio - ParetoFront › 20261109_GRANAR_MECHA_RootSlice › outputs › `Cost-and-Conductivity_outputs_20261109.zip` › `outputs/deliverables/` |
| Frozen reference results, for `reproduce/compare_to_reference.py` (set `COSTCOND_REF` to the downloaded folder) | Plant.Architecture.Lab › Magda › NSF Math Bio - ParetoFront › 20261109_GRANAR_MECHA_RootSlice › outputs › `Cost-and-Conductivity_outputs_20261109.zip` (unzip; the `outputs/` folder inside is the reference) |

The data on Box are shared under CC BY 4.0 (see `../LICENSE-data.md`). Please cite this repository when you reuse them.
