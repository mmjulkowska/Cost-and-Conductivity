# 2. MECHA: radial conductivity and axial conductance

MECHA solves water flow through the cell network of each GRANAR anatomy. It returns:
- **kr**: radial conductivity, cm hPa⁻¹ d⁻¹, per unit root surface
- **Kx**: axial conductance, cm⁴ hPa⁻¹ d⁻¹, per unit pressure gradient (Poiseuille flow in the xylem)

## Software
- `openalea.mecha` from MECHA-dev `develop` @ 945b2be.
- Install with `env/install_mecha_mac.sh`. The installer applies one fix: in `utils/hydraulic_solver.py` the helper `_contact_nodes()` was merged under the wrong name.

## Barrier per segment
Assigned from toluidine-blue staining of Emily Miller's sections.

| Segment | Endodermis | Exodermis | MECHA code |
|---|---|---|---|
| MP, MM | lignin (Casparian strip) | lignin + suberin | B6 |
| MT | lignin | none | B1 |
| MT0, main true tip (MT anatomies) | none | none | no barrier* |
| LP, LM | not detected | lignin + suberin | B10 |
| LT | none | none | no barrier* |

\*No barrier = a mature-xylem stage with every barrier wall conductivity set to the primary-wall kw (recipe from V. Couvreur). Because the xylem is mature, MT0 Kx is an upper bound.

## Hydraulic parameters (MECHA-dev defaults)
| Parameter | Value |
|---|---|
| Primary wall kw | 2.4 × 10⁻⁴ cm² hPa⁻¹ d⁻¹ |
| Barrier walls | 1 × 10⁻¹⁶ |
| Septa | 1.2 × 10⁻⁴ |
| Membrane kmb | 3 × 10⁻⁵ cm hPa⁻¹ d⁻¹ |
| Aquaporin kAQP | 4.3 × 10⁻⁴ cm hPa⁻¹ d⁻¹ |
| Plasmodesma Kpl | 5.3 × 10⁻¹² cm³ hPa⁻¹ d⁻¹ |
| Cell length | 200 µm |
| `im_scale` | 1000 (anatomies are in mm) |

## Scripts
| Script | Writes |
|---|---|
| `run_one.py` | one MECHA run (anatomy, barrier code, no-barrier option) |
| `batch.py` | `outputs/mecha/mecha_raw.csv`, `conductivities_per_anatomy.csv` (228 rows) |
| `analyse.py` | `outputs/mecha/stats/` (ANOVA on log10, Tukey HSD, Levene, Welch), `anatomy_vs_emily.csv`, `anatomy_by_segment.csv`, `outputs/figures/Fig_conductivity_barriers_v5*` |
| `sensitivity_mc/` | Monte Carlo of kr against kw, kmb, kAQP and Kpl (Latin hypercube, ×½–2 and ×0.1–10). Run `run_mc.sh` (several hours); results in `outputs/mecha/sensitivity_mc/` |
| `checks/lignin_cap.py` | exodermis as a lignin cap (B9) instead of suberin + lignin |
| `checks/bracket_laterals.py` | LP, LM with an endodermal strip (B6) instead of B10 |

```
python 2_mecha/batch.py && python 2_mecha/analyse.py
```

## Results (default parameters)
- **Statistics:** kr F(6,221) = 1419 and Kx F(5,202) = 53, both on log10 values. MT0 is left out of the Kx test because it shares the MT anatomies.
- **Monte Carlo:** segment kr spans ×2.3–2.5 (95 % range). MP and MM are the lowest segments in every draw. Kx does not depend on the cell-scale parameters.
- **Lignin-cap check:** kr ×1.09–1.10 for MP and MM, ×1.22 for LP and LM.
- **Lateral endodermis bracket:** lateral kr falls to ×0.53–0.56, the same level as MM and MP.
