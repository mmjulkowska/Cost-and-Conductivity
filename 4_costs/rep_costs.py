"""Costs per replicate (RootSlice), paired with MECHA conductivities; segment cost statistics.
Run from anywhere after RootSlice (3_rootslice/) and MECHA (2_mecha/batch.py). Outputs in <outputs>/costs/.
Writes replicate_costs.csv (all RootSlice replicates), replicate_level_costs_conductivity.csv (paired),
section_costs.csv (replicate means per section) and costs_<var>_summary_letters.csv / _tukey.csv / costs_anova.json
(ANOVA on log10 section means, n = sections)."""
import os, sys; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from paths import *
import pandas as pd, numpy as np, json, sys, os
sys.path.insert(0, os.path.join(ROOT, 'shared')); from stats_conductivity import cld
os.chdir(COSTS)  # cost tables live in <outputs>/costs/

from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.oneway import anova_oneway
ORDER = ['MP', 'MM', 'MT', 'LP', 'LM', 'LT']
seg = lambda n: 'MP' if n <= 9 else 'MM' if n <= 20 else 'MT' if n <= 26 else 'LP' if n <= 40 else 'LM' if n <= 50 else 'LT'
v = pd.read_csv(RS_VOL)
# 13 of the 19 replicates that hung/crashed were recovered with a simplified stele interior (run_small_stele.sh; controls gave identical cortical volumes)
_ss = RS_VOL_SS
if os.path.exists(_ss):
    _s = pd.read_csv(_ss); _s = _s[_s.model.str.startswith('R5_') & ~_s.model.isin(v.model)]; v = pd.concat([v, _s], ignore_index=True)
o = pd.DataFrame({'anatomy': v.model.str.replace('R5_', '', regex=False)})
o['section'] = o.anatomy.str.extract(r'MP(\d+)_')[0].astype(int); o['segment'] = o.section.map(seg)
k = 1000 / v['TotalHeight(um)'] / v['cortex2SectionRatio']
o['C_ng_mm'] = v['corticalCellWallVolume(um3)'] * 1e-12 * 1.51 * 0.45 * 1e9 * k
o['Ncyt_ng_mm'] = 1.8348777e-15 * v['corticalCytoplasmVolume(um3)'] * 1e9 * k
o['Nvac_ng_mm'] = 4.20201e-15 * v['corticalVacuoleVolume(um3)'] * 1e9 * k
o['N_ng_mm'] = o.Ncyt_ng_mm + o.Nvac_ng_mm
o['maintC_ng_mm_d'] = v['corticalCytoplasmVolume(um3)'] * 1e-12 * 50e-9 * 12.011 * 86400 * 1e9 * k
o['cortex_fraction_model'] = v['cortex2SectionRatio']
cfg = pd.read_csv(RS_VALUES)[['anatomy', 'stele_r', 'cortex_cell', 'layers', 'epidermis', 'endodermis', 'stele_cell']]
o = o.merge(cfg, on='anatomy').sort_values(['section', 'anatomy']); o.to_csv('replicate_costs.csv', index=False)
allcfg = pd.read_csv(RS_VALUES)
print(len(o), 'of', len(allcfg), 'RootSlice replicates ran; missing:', ' '.join(sorted(set(allcfg.anatomy) - set(o.anatomy))))
# pair with conductivities
kk = pd.read_csv(COND); rows = []
for a, r in kk.groupby('anatomy'):
    main = r[r.fragment != 'MT0'].iloc[0]; t0 = r[r.fragment == 'MT0']
    rows.append(dict(anatomy=a, barrier=main.mecha_barrier, kr=main.kr_cm_hPa_d, Kx=main.Kx_cm4_hPa_d, kr_trueTip_MT0=(t0.kr_cm_hPa_d.iloc[0] if len(t0) else np.nan)))
P = o.merge(pd.DataFrame(rows), on='anatomy', how='inner').sort_values(['section', 'anatomy'])
P.to_csv('replicate_level_costs_conductivity.csv', index=False)
print(len(P), 'paired replicates;', P.section.nunique(), 'sections')
print(P.groupby('segment').agg(reps=('anatomy', 'count'), sections=('section', 'nunique')).reindex(ORDER).T)
# section means (all RootSlice replicates) and stats on sections
VARS = ['C_ng_mm', 'maintC_ng_mm_d', 'Ncyt_ng_mm', 'Nvac_ng_mm', 'N_ng_mm']
S = o.groupby(['segment', 'section'])[VARS + ['stele_r', 'cortex_cell', 'layers']].mean().reset_index(); S['n_reps'] = o.groupby(['segment', 'section']).size().values
S.to_csv('section_costs.csv', index=False); res = {}
for var in VARS:
    d = S.copy(); d['y'] = np.log10(d[var]); smp = [d.y[d.segment == g].values for g in ORDER]
    F, p = stats.f_oneway(*smp); lev = stats.levene(*smp, center='median'); w = anova_oneway(smp, use_var='unequal')
    tk = pairwise_tukeyhsd(d.y, d.segment); t = pd.DataFrame(tk._results_table.data[1:], columns=tk._results_table.data[0])
    t['fold_change(group2/group1)'] = 10 ** t.meandiff.astype(float); t.to_csv(f'costs_{var}_tukey.csv', index=False)
    L = cld(ORDER, {(r.group1, r.group2): bool(r.reject) for r in t.itertuples()})
    s = d.groupby('segment')[var].agg(n='count', mean='mean', sd='std', median='median').reindex(ORDER); s.index.name = 'fragment'; s['tukey_letter'] = [L[g] for g in ORDER]
    s.to_csv(f'costs_{var}_summary_letters.csv')
    res[var] = dict(F=float(F), dfb=len(ORDER) - 1, dfw=len(d) - len(ORDER), p=float(p), levene_p=float(lev.pvalue), welch_F=float(w.statistic), welch_p=float(w.pvalue))
    print(var, f"F({len(ORDER)-1},{len(d)-len(ORDER)})={F:.1f} p={p:.2g}", ' '.join(f"{g}:{s.loc[g,'mean']:.0f}{L[g]}" for g in ORDER))
json.dump(res, open('costs_anova.json', 'w'), indent=1)
