"""v5 statistics, anatomy check against Emily Miller's measurements, and the conductivity figure with legend.
Run from anywhere after 2_mecha/batch.py.
Writes <outputs>/mecha/stats/ (ANOVA/Tukey tables, anova_results.json), <outputs>/mecha/anatomy_vs_emily.csv, anatomy_by_segment.csv,
<outputs>/figures/Fig_conductivity_barriers_v5[_withLegend].png/.pdf and Fig_conductivity_barriers_v5_legend.md."""
import os, sys; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from paths import *
import pandas as pd, numpy as np, json, os, sys
sys.path.insert(0, os.path.join(ROOT, 'shared'))
from stats_conductivity import analyse, ORDER
from fig_conductivity import draw
F = 1.03214286e-7
SEGS = ['MP', 'MM', 'MT', 'LP', 'LM', 'LT']
CSV = COND
df = pd.read_csv(CSV).rename(columns={'Kx_cm4_hPa_d': 'Kx', 'kr_cm_hPa_d': 'kr'})
R = {}; S = {}
for var, groups in (('kr', ORDER), ('Kx', [g for g in ORDER if g != 'MT0'])):
    s, t, r = analyse(df, var, groups, 'v5_AllSections')
    s.to_csv(os.path.join(STATS, f'{var}_summary_letters.csv')); t.to_csv(os.path.join(STATS, f'{var}_tukey.csv'), index=False)
    R[var] = r; S[var] = s
    print(f"== v5 {var}: F({r['df_between']},{r['df_within']}) = {r['F']:.1f}, p = {r['p']:.2g}; Levene p = {r['levene_p']:.2g}; Welch p = {r['welch_p']:.2g}")
    print(s.to_string(float_format=lambda x: f'{x:.3g}'))
json.dump(R, open(os.path.join(STATS, 'anova_results.json'), 'w'), indent=1)

# ---- anatomy vs Emily (per anatomy and per segment)
G = pd.read_csv(GRANAR_SUMMARY); G = G[G.ok]
V = pd.read_csv(REPLICATES)
X = pd.read_excel(SHEET); X.insert(0, 'MPn', range(1, len(X) + 1))
em = pd.DataFrame({'MPn': X.MPn, 'em_root_mm2': pd.to_numeric(X['Root Area (pixels)'], errors='coerce') * F,
                   'em_cortex_mm2': pd.to_numeric(X['Cortex Total (pixels)'], errors='coerce') * F})
A = G.merge(V, on='id').merge(em, on='MPn')
A['sim_cortex_mm2'] = A.area_cortex + A.area_exo
A['root_ratio'] = A.area_all / A.em_root_mm2; A['cortex_ratio'] = A.sim_cortex_mm2 / A.em_cortex_mm2
A['root_d_um'] = 2 * np.sqrt(A.area_all / np.pi) * 1000
A.to_csv(os.path.join(MECHA, 'anatomy_vs_emily.csv'), index=False)
AS = A.groupby('segment').agg(n_anat=('id', 'size'), n_sections=('MPn', 'nunique'), root_mm2=('area_all', 'mean'), root_d_um=('root_d_um', 'mean'),
                              em_root_mm2=('em_root_mm2', 'mean'), root_ratio_median=('root_ratio', 'median'), cortex_ratio_median=('cortex_ratio', 'median'),
                              layers=('emily_layers', 'mean'), cortex_cells=('n_cortex_cells', 'mean'), exo_cells=('n_exo_cells', 'mean'),
                              cortex_cell_d_um=('cortex_mm', lambda x: x.mean() * 1000), stele_d_um=('stele_layer_mm', lambda x: x.mean() * 1000)).reindex(SEGS)
AS.to_csv(os.path.join(MECHA, 'anatomy_by_segment.csv'))
print(AS.to_string(float_format=lambda x: f'{x:.3g}'))
print('overall median root ratio', A.root_ratio.median().round(2), 'cortex ratio', A.cortex_ratio.median().round(2))

# ---- figure + legend
SUP = str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹')
def sci(x): m, e = f'{x:.2e}'.split('e'); return f'{float(m):.2f} × 10' + str(int(e)).translate(SUP)
def pf(p): return 'p < 1 × 10⁻¹⁵' if p < 1e-15 else f'p = {sci(p)}'
def means(v): return '; '.join(f"{f} {sci(r['mean'])} ± {sci(r['sd'])} (n = {int(r['n'])}, {r['tukey_letter']})" for f, r in S[v].iterrows())
anat = '; '.join(f"{sg}: {int(r.n_sections)} sections / {int(r.n_anat)} anatomies, root {r.root_mm2:.3f} mm² ({r.root_d_um:.0f} µm; Emily {r.em_root_mm2:.3f} mm², median ratio {r.root_ratio_median:.2f}), "
                 f"{r.layers:.1f} cortical layers incl. exodermis, {r.cortex_cells:.0f} cortex + {r.exo_cells:.0f} exodermis cells, cortex cell input {r.cortex_cell_d_um:.1f} µm, stele {r.stele_d_um:.0f} µm"
                 for sg, r in AS.iterrows())
failed = ', '.join(pd.read_csv(GRANAR_SUMMARY).query('not ok').id)
excl = V[~V.used]
a, b = R['kr'], R['Kx']
L = ("Figure. Root hydraulic conductivity of tomato (S. pimpinellifolium) root fragments under the assigned apoplastic barriers — all usable sections of Emily Miller's measurements.\n\n"
     "(A) Schematic of a main root with one lateral root. Sampled positions are 3 cm apart: main true tip (MT0), tip zone 3 cm from the tip (MT), middle 3 cm above MT (MM), proximal 3 cm above MM (MP); "
     "lateral true tip (LT), lateral middle (LM), lateral proximal (LP). Outer ring = exodermis, inner ring = endodermis; coloured = barrier present, dashed grey = absent. "
     "Barriers (Magda Julkowska, 6 Nov 2026, from Emily Miller's staining): MP, MM endodermal lignin + exodermal suberin and lignin (MECHA B6); MT endodermal lignin only (B1); "
     "LP, LM exodermal suberin and lignin only (B10); MT0, LT none (mature-xylem stage with every barrier wall conductivity set to primary-wall kw, after V. Couvreur).\n\n"
     "(B) Radial conductivity kr (cm hPa⁻¹ d⁻¹) and (C) xylem axial conductance Kx (cm⁴ hPa⁻¹ d⁻¹), log scales. Points are individual anatomies; bars are medians; letters are Tukey HSD groups. "
     "Each anatomy is one measured cell (replicate 1–4) of one of Emily Miller's cross-sections.\n\n"
     "Anatomies (GRANAR 0.0.1.2, built directly from Emily Miller's sheet AllDatammAndpixelsNov2026_2.xlsx at 1.032 × 10⁻⁷ mm² px⁻¹; diameters = 2√(area/π)): for replicate k, the cortex, endodermis and epidermis cell "
     "diameters are Emily's cell #k of that section, and the stele cell is her vasculature cell #k (segment median where not measured); stele diameter from her vasculature area; "
     "cortex layers = her count with the exodermis counted as one of them. Pericycle 21 µm, exodermis 20 µm, max xylem 35 µm for all (template values). "
     ""
     f"Seed = 10 × section + cell (+100 for MP8_cell3, MP10_cell2 and MP16_cell3 after singular hydraulic solves). Excluded: rows 4 (DSC_00425PR, torn) and 60 (DSC_0030RT, no stele visible), and {len(excl[~excl.MPn.isin([4,60])])} further replicates lacking a measurement; "
     f"GRANAR failed to tessellate {failed} with 11 seeds. Mean anatomy properties: {anat}.\n\n"
     f"Means ± SD (Tukey group). kr: {means('kr')}. Kx: {means('Kx')}. MT0 Kx equals MT.\n\n"
     f"Statistics: one-way ANOVA on log10 values, Tukey HSD (α = 0.05). kr: F({a['df_between']},{a['df_within']}) = {a['F']:.1f}, {pf(a['p'])}; Levene p = {sci(a['levene_p'])}; Welch F = {a['welch_F']:.1f}, {pf(a['welch_p'])}. "
     f"Kx: F({b['df_between']},{b['df_within']}) = {b['F']:.1f}, {pf(b['p'])}; Levene p = {b['levene_p']:.2g}; Welch F = {b['welch_F']:.1f}, {pf(b['welch_p'])}. "
     "MT and MT0 are the same anatomies with and without barriers (not independent); MT0 is omitted from the Kx test (grey italic letter repeats MT).\n\n"
     "Model: MECHA (openalea.mecha, MECHA-dev develop 945b2be with a one-line merge fix), defaults: kw 2.4 × 10⁻⁴ cm² hPa⁻¹ d⁻¹, barrier walls 1 × 10⁻¹⁶, septa 1.2 × 10⁻⁴, kAQP 4.3 × 10⁻⁴ cm hPa⁻¹ d⁻¹, "
     "Kpl 5.3 × 10⁻¹² cm³ hPa⁻¹ d⁻¹, cell length 200 µm, im_scale 1000.")
open(os.path.join(FIG, 'Fig_conductivity_barriers_v5_legend.md'), 'w').write(L)
lk, lx = S['kr'].tukey_letter.to_dict(), S['Kx'].tukey_letter.to_dict()
ttl = "Root hydraulic conductivity by segment and barrier "
draw(CSV, lk, lx, ttl, os.path.join(FIG, 'Fig_conductivity_barriers_v5'))
draw(CSV, lk, lx, ttl, os.path.join(FIG, 'Fig_conductivity_barriers_v5_withLegend'), legend_text=L)
print('done')
