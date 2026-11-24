# Step 6: CPlantBox root architectures, their Krs, and water uptake in drying soil

Builds three root-system topologies of *S. pimpinellifolium* seedlings in CPlantBox, from the measured panel:
**Christmas tree** (C), **Broomstick** (B) and **Droopy telephone pole** (dT), following the panel's `architectural_type` codes.
Every comparison is made at **equal total root length** (30 cm, about the control median at day 5), so the types differ
in shape, not size.

The step then:
- gives every root piece the step-2 (v5) segment conductivities and computes whole-root-system conductance (Krs);
- grows the root systems in pots of three soils, with and without drought.

Step 5 computes Krs for each traced plant of the panel. This step asks the complementary question: how much does the
topology itself matter, with size held equal?

## Inputs (not in git)

| What | Where |
|---|---|
| `01_panel_root_geometry.csv`, `pooled_plants_by_day.csv`, `pooled_laterals.csv`, `timeseries_roots.csv` | `data/krs/` (`$KRS_DATA`): the step-5 inputs, made by `5_krs_panel/prep/` |
| PimpiBig2 day-5 SmartRoot RSML (300 dpi) | `$PIMPIBIG2_RSML` (default `~/mnt/Pimp_RSA/SmartRooted/pimpiBIG2`; Dropbox `Tomato_Pimpinellifolium/Pimp_RSA/SmartRooted/pimpiBIG2`) |
| Segment conductivities, radii, lateral-endodermis ratio | step-2 outputs in `outputs/mecha/` (from Box, see `data/README.md`) |

The repository also holds `inputs/rsml_selected.csv`. It lists the 72 PimpiBig2 day-5 tracings used (B 14, C 28, dT 30; control and salt).

**Software:** CPlantBox 2.2.0 (`pip install CPlantBox`; Python 3.13 wheel), numpy, scipy, pandas, matplotlib, ffmpeg (animations).

## Run

```bash
bash 6_cplantbox/run_all.sh              # stages 1 3 4 5 with the stored fitted parameters (params/), about 15 min
FIT=1 bash 6_cplantbox/run_all.sh        # refit every parameter set first (about 1 h)
bash 6_cplantbox/run_all.sh 6            # drought runs in pots: 162 runs, about 1 h on 2 cores
```

Outputs go to `outputs/cplantbox/`.

**Reproduction check (24 November 2026):** a rerun from this folder with the stored parameters gave the original outputs exactly. That covers:
- the targets (`rsml_*.csv`, `kinetics_targets*.json`);
- the simulated replicate metrics of v1, v2 and the exaggerated set;
- the Krs tables;
- a pot run.

## Scripts

| Stage | Script | Does | Main outputs |
|---|---|---|---|
| 1 | `parse_rsml.py` | PimpiBig2 day-5 RSML: lateral insertion, length, angles | `rsml_plants.csv`, `rsml_laterals.csv` |
| 1 | `kinetics.py` | PimpiBig2 days 1–5, roots followed by SmartRoot ID: main-root and lateral growth by type | `timeseries_*.csv`, `kinetics_targets.json` |
| 1 | `track_laterals.py` | PimpiBig5/6: laterals followed day to day by insertion position (SmartRoot IDs not kept there) | `pooled_lateral_rates.csv` |
| 1 | `targets_v3.py` | pooled targets for the v3 fit | `kinetics_targets_v3.json` |
| 2 | `fit.py`, `fit_angle.py`, `fit_angle_plate.py` | v1: day-5 shape and angles | `fit_shape.json`, `params_fitted.json` |
| 2 | `make_params_exaggerated.py` (`try_exag.py` to explore) | exaggerated topologies | `params_exaggerated.json` |
| 2 | `fit_v2.py`, `fit_angle_v2.py` | v2: PimpiBig2 day 1–5 kinetics, accelerating lateral growth | `fit_v2_<type>.json`, `params_v2.json` |
| 2 | `fit_v3.py`, `fit_angle_v3.py` | v3: pooled PimpiBig1–6 | `fit_v3_<type>.json`, `params_v3.json` |
| 2 | `make_params_v4.py`, `refit_v4.py` | v4: constant-rate main root, lateral maximum by type | `params_v4.json` |
| 3 | `simulate_equal_trl.py`, `plot_gallery.py`, `plot_validation.py`, `plot_profiles.py`, `plot_kinetics.py` | 12 replicates per type grown to 30 cm, 2D plate and 3D; galleries and checks against the data | `sims*/`, `figs*/` |
| 4 | `animate.py` | growth from germination to tracing day 12 (v3, v4) | `anim/`, `anim_v4/` |
| 5 | `conductance.py`, `krs_analysis.py`, `plot_krs.py`, `animate_krs.py` | Krs with the v5 conductivities | `krs/` |
| 6 | `soil_sim.py`, `soil_batch.py`, `plot_soil.py`, `animate_soil.py` | water uptake from pots, three soils, three water scenarios | `soil_runs/`, `soil/` |

**Shared modules:**
- `cpb_paths.py`: paths; every script changes into `outputs/cplantbox/`.
- `rsa_model.py`: builds the CPlantBox plant.
- `growth.py`: the accelerating lateral growth function.

`params/` holds the fitted parameters of the November 2026 runs, which stage 2 copies when `FIT` is not set.

## Model

### Architecture

| Item | Setting |
|---|---|
| Root orders | Main root plus first-order laterals; the day-5 plants have no second-order laterals |
| Main root, free parameters | growth rate, apical unbranched zone (`la`), inter-lateral distance with a depth slope (`ln`, `lnk`) |
| Laterals, free parameters | growth, insertion angle, gravitropism |
| Geometry: 3D | the measured insertion angle is used |
| Geometry: 2D plate | a 0.5 mm slab; it deflects laterals, so on the plate the angles are calibrated to the plate image measurements |

**v1, day-5 shape** (`params_fitted.json`): fitted to the panel day-5 medians (main-root length, lateral number, total lateral length; n = 122–203 control plants per type) and to the PimpiBig2 RSML lateral-position profiles.
- **Weakness:** it can only reach the measured lengths with constant lateral growth rates of 0.6–1.5 cm d⁻¹.

**v2, growth kinetics** (`params_v2.json`): adds the PimpiBig2 day 1–5 tracings of the same plants.
- **Laterals speed up as they age:** about 0.5 cm d⁻¹ on their first day, rising to 0.9–2.0 by the fourth.
- **New growth law:** laterals grow at rate rL + sL·age (`growth.py`). The fitted rL of 0.3–0.4 cm d⁻¹ matches the ISA3 plates (about 0.36).
- **Root age at tracing day 0:** fitted at 0.6–0.9 d.

**v3, pooled PimpiBig1–6** (`params_v3.json`): lateral-position profiles from all six experiments (121–197 control plants per type); time courses from Big2, Big5 and Big6.
- **Big5 day numbering:** days follow the image dates (20 Oct = day 2, 22 Oct = day 4, 23 Oct = day 5). Its 19 Oct tracings are dropped because their main roots are longer than on the next day.
- **Lateral tracking:** matching by insertion position agrees with the SmartRoot ID for 99.8% of 3,801 laterals (Big1/2, days 4→5).

| v3 | Christmas tree | Broomstick | Droopy telephone pole |
|---|---|---|---|
| Root age at tracing day 0 (d) | 0.51 | 0.32 | 0.23 |
| Main-root initial rate (cm d⁻¹) | 1.50 | 1.58 | 2.67 |
| Main-root lmax (cm) | 31.7 | 9.9 | 16.8 |
| Apical unbranched zone (cm) | 3.12 | 2.38 | 4.12 |
| Mean inter-lateral distance (cm) | 0.77 | 0.45 | 0.41 |
| Lateral rate at emergence rL (cm d⁻¹) | 0.32 | 0.34 | 0.40 |
| Lateral acceleration sL (cm d⁻²) | 0.22 | 0.33 | 0.12 |

**v4, current set** (`params_v4.json`), as requested by M. Julkowska. Two changes from v3, each followed by a refit of the apical zone, rL and sL only (`refit_v4.py`), so that days 1–5 still match:
1. **Main root grows at a constant rate** (fit to days 1–5): C 1.33, B 1.11, dT 1.96 cm d⁻¹. It keeps branching beyond the old lmax, and the spacing keeps widening as before.
2. **Maximal lateral length by type:** B 8 cm > C 6 cm > dT 4 cm. These values were set, not fitted. They lie near the 99th percentile of measured day-5 lateral length (7.6, 5.3 and 5.1 cm).
- **Cost of the constant rate:** the dT main root overshoots its slowing between days 4 and 5 (10.3 cm vs 8.8 cm measured).

**Exaggerated set** (`params_exaggerated.json`): hand-set from v1 to sharpen the contrasts, so it is not fitted to data:
- **Broomstick:** 6 laterals, mean 3.7 cm.
- **Christmas tree:** lateral length grades from about 4 cm at the top to 0.1 cm at the bottom.
- **Droopy telephone pole:** all laterals about 1.2 cm.

### Root conductance (Krs)

**Which conductivity each root piece gets** (`conductance.py`): the v5 class, by distance from that root's own tip. Samples were taken at the tip and at 3, 6 and 9 cm from it, so the boundaries are the midpoints.
- **Main root:** MT0 (0–1.5 cm, no barrier), MT (1.5–4.5 cm, B1), MM (4.5–7.5 cm, B6), MP (> 7.5 cm, B6).
- **Laterals:** LT (0–1.5 cm, no barrier), LM (1.5–4.5 cm, B10), LP (> 4.5 cm, B10).

**Values:**
- **kr and Kx:** the medians per segment from `outputs/mecha/conductivities_per_anatomy.csv`.
- **Radius:** Emily Miller's measured cross-section area.
- **Lateral endodermis variant (`lat_endo`):** LP/LM kr is multiplied by the B6/B10 ratio from `outputs/mecha/checks/`.

**Solver:** a Doussan-type network. Radial conductance is 2πa·l·kr, split between a segment's two nodes; axial conductance is Kx/l. It matches the analytic single-root solution exactly. Step 5's `krs_tree.py` agrees with it to within 0.15% on the same plants.

### Soil

**Plant and pot** (`soil_sim.py`):
- **Plants:** v4 architectures, grown to 30 cm total root length inside a pot 2.5 cm wide × 20 cm deep, and confined to it.
- **Two modes:** the root system is either fixed from then on, or keeps growing.

**Soil water:** a 1D vertical Richards equation in 0.5 cm layers, with a closed bottom.
- **Soils:** sandy loam, loam and silty clay loam, with the van Genuchten–Mualem parameters of Carsel & Parrish (1988).

**Root–soil coupling:** each root piece connects to its layer's soil through a perirhizal conductance, 2π l K(h)/ln(r_b/a).
- **r_b:** half the mean root spacing in that layer.
- **K:** taken at the layer's bulk pressure head.

**Demand:**
- **Transpiration:** 3 cm³ d⁻¹ for every plant, on a day/night sinusoid.
- **Stress:** when the collar would fall below −15000 hPa, it is held there and transpiration drops.

**Water scenarios:**

| Name | What it means |
|---|---|
| Well watered | Water added at the surface as fast as it is taken up |
| Dry-down | Starts from a drained pot (container capacity), no water added |
| Dry top | Top 6 cm starts at −15000 hPa, wet below, no water added |

**Check:** the water balance closes to < 0.1 cm³ in all 162 runs.

## Results (November 2026)

**Architecture at equal size**
- **The three types keep their measured differences at 30 cm total root length** (v2 parameters, 3D, median of 12 replicates):

  | Type | Main root (cm) | Laterals | Mean lateral (cm) |
  |---|---|---|---|
  | Christmas tree | 7.4 | 16 | 1.4 |
  | Broomstick | 5.3 | 12 | 2.1 |
  | Droopy telephone pole | 8.8 | 20 | 1.1 |

- **Equal age as well as equal size:** with the v2/v3 kinetics, all three reach 30 cm at nearly the same age.

**Krs at 30 cm total root length** (×10⁻⁴ cm³ hPa⁻¹ d⁻¹): droopy telephone pole 10.2 > Christmas tree 8.2 > broomstick 7.8.
- **The step-5 panel agrees:** the traced control panel at day 5 has a median of 8.5.
- **Lateral tips take up most of the water:** 78–88% enters through barrier-free lateral tips and 4–6% through the main-root tip.
- **A lateral endodermal barrier changes Krs little:** −1 to −4% at 30 cm, rising to −8 to −16% by day 12 (extrapolated).
- **The largest uncertainty is the length of the barrier-free tip zone:** about −25 to −30% for 0.5 cm and +30% for 3 cm, at day 12.

**Drought:** the ranking is the same in every soil and drought scenario: droopy telephone pole > Christmas tree > broomstick.
- **Size of the topology effect:** 5–10 percentage points of 10-day transpiration (relative to demand) at equal root length.
- **Continued root growth matters more:** it adds 15–25 points.
- **Well watered:** every plant meets demand, but the collar reaches about −10400 hPa at peak, so the roots' own conductance is close to limiting.
- **Dry top in small pots:** capillary rise rewets the dry top within about a day in loam and sandy loam.

| 10-day transpiration (% of demand), root system fixed at 30 cm | Christmas tree | Broomstick | Droopy telephone pole |
|---|---|---|---|
| Loam, dry-down | 67 | 64 | 70 |
| Sandy loam, dry-down | 52 | 49 | 56 |
| Silty clay loam, dry-down | 45 | 40 | 49 |
| Loam, dry top | 43 | 39 | 48 |
| Silty clay loam, dry top | 24 | 20 | 29 |

## Limitations

**Architecture**
- **No second-order laterals.** Tomato forms them in the second week.
- **Beyond tracing day 5 everything is extrapolation:** day-12 architectures, Krs and the growing pot runs.
- **Fitted to the median plant of each type.** Replicates vary less than the panel does.

**Hydraulics**
- **Tip Kx (MT0, LT) assumes mature xylem,** so it is an upper bound.
- **kr depends only on distance from the tip,** not on age, soil or salt.

**Soil**
- **The soil model is 1D.** Depletion around single roots enters only through the perirhizal term, with K at bulk pressure.
- **Pot size and demand are seedling-scale choices.**
- **No root growth response to drought, and no osmotic adjustment.**

## References
- CPlantBox: Schnepf et al. 2018, *Annals of Botany* 121:1033–1053; Zhou et al. 2020, *in silico Plants* 2:diaa001.
- Root hydraulic network: Doussan et al. 1998, *Annals of Botany* 81:213–223.
- Soil hydraulic parameters: Carsel & Parrish 1988, *Water Resources Research* 24:755–769.
