# 4. Costs, pairing with conductivity, and the table for the Pareto model

All outputs are written to `outputs/costs/`, and the workbook to `outputs/deliverables/`.

## `rep_costs.py`
Scales every replicate to 1 mm of root: k = 1000 / TotalHeight (µm) / cortex2SectionRatio.

| Output | Formula |
|---|---|
| Building C (ng C mm⁻¹) | cortical cell-wall volume × 1.51 g cm⁻³ × 0.45 g C g⁻¹ |
| Nitrogen (ng N mm⁻¹) | 1.8348777 × 10⁻¹⁵ g µm⁻³ × cytoplasm volume + 4.20201 × 10⁻¹⁵ g µm⁻³ × vacuole volume |
| Maintenance (ng C mm⁻¹ d⁻¹) | cytoplasm volume × 50 nmol CO₂ s⁻¹ cm⁻³ × 12.011 g mol⁻¹ × 86 400 s d⁻¹ |

Writes:
- `replicate_costs.csv`: all 210 RootSlice replicates.
- `replicate_level_costs_conductivity.csv`: 207 replicates paired with kr and Kx. MT rows also carry the MT0 kr.
- `section_costs.csv`: replicate means per section (53 sections).
- `costs_<variable>_summary_letters.csv`, `_tukey.csv` and `costs_anova.json`: one-way ANOVA on log10 section means, with n = sections.

Building C by section, ng C mm⁻¹ (Tukey group): MP 5540 (a), MM 2723 (b), MT 2268 (bc), LP 2026 (bc), LM 1990 (bc), LT 1760 (c).

## `rep_stats.py`, `fig_rep.py`
`rep_stats.py` computes cost vs conductivity correlations three ways:
- Spearman's rank correlation;
- a mixed model with cross-section as a random effect;
- a repeated-measures correlation within cross-sections.

Within sections, building C vs kr gives r_rm = −0.58; building C vs Kx is about 0. `fig_rep.py` plots the figure, written to `Fig_cost_vs_conductivity_replicates*`.

## `build_xlsx.py`
Writes `outputs/deliverables/Root_costs_conductivity_for_Arjun_v5_Nov2026.xlsx` (segment means, replicates, statistics, notes). If the Monte Carlo results are present in `outputs/mecha/sensitivity_mc/`, it also adds a "kr uncertainty" sheet.

```
python 4_costs/rep_costs.py && python 4_costs/rep_stats.py && python 4_costs/fig_rep.py && python 4_costs/build_xlsx.py
```
