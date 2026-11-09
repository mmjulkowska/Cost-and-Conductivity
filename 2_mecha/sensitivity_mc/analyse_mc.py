"""Analyse the kr Monte Carlo (mc_kr.py). Run after the three mc_kr.py designs finished (files in <outputs>/mecha/sensitivity_mc/).
Writes mc_summary_by_barrier.csv, mc_barrier_ratios.csv, mc_elasticities.csv, mc_segment_means.csv, mc_rank_checks.json,
mc_cost_kr_correlation_per_draw.csv; figure via fig_mc.py."""
import os, sys; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from paths import *
import pandas as pd, numpy as np, json, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(ROOT, 'shared')); from stats_conductivity import cld
os.chdir(MC)  # all Monte Carlo files live in <outputs>/mecha/sensitivity_mc/

from statsmodels.stats.multicomp import pairwise_tukeyhsd
import statsmodels.formula.api as smf
from scipy import stats
P = ['kw', 'kmb', 'kaqp', 'kpl']; ORDER = ['MP', 'MM', 'MT', 'MT0', 'LP', 'LM', 'LT']; TR = ['none', 'B1', 'B10', 'B6']
REPSEG = {}
out_sum, out_rat, out_el = [], [], []
for r in (2, 10):
    d = pd.read_csv(f'mc_rep_x{r}.csv'); d['range'] = f'x{r}'
    for t, g in d.groupby('treatment'):
        X = np.column_stack([np.ones(len(g))] + [np.log(g['f_' + p]) for p in P] + [pd.get_dummies(g.anatomy).values[:, 1:].astype(float)])
        y = np.log(g.kr.values); b = np.linalg.lstsq(X, y, rcond=None)[0]; R2 = 1 - ((y - X @ b) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        out_el.append(dict(range=f'x{r}', treatment=t, **{f'elasticity_{p}': b[i + 1] for i, p in enumerate(P)}, R2=R2))
        for (a, frag), h in g.groupby(['anatomy', 'fragment']):
            out_sum.append(dict(range=f'x{r}', treatment=t, fragment=frag, anatomy=a, kr_median=h.kr.median(), kr_q025=h.kr.quantile(.025), kr_q975=h.kr.quantile(.975),
                                fold_95=h.kr.quantile(.975) / h.kr.quantile(.025)))
    w = d.pivot_table(index=['draw', 'anatomy', 'fragment'], columns='treatment', values='kr')
    for a, b_ in [('B10', 'B1'), ('B6', 'B10'), ('B6', 'B1'), ('none', 'B6'), ('none', 'B1'), ('none', 'B10')]:
        x = w[a] / w[b_]; out_rat.append(dict(range=f'x{r}', ratio=f'{a}/{b_}', median=x.median(), q025=x.quantile(.025), q975=x.quantile(.975), frac_gt1=(x > 1).mean()))
    # what decides B10 vs B1: correlation with log(Kpl/kAQP)
    x = np.log((w['B10'] / w['B1']).values); z = np.log(d.drop_duplicates('draw').set_index('draw').loc[w.index.get_level_values('draw'), 'f_kpl'].values /
                                                     d.drop_duplicates('draw').set_index('draw').loc[w.index.get_level_values('draw'), 'f_kaqp'].values)
    out_rat.append(dict(range=f'x{r}', ratio='B10/B1 vs Kpl/kAQP (pearson r of logs)', median=np.corrcoef(x, z)[0, 1]))
pd.DataFrame(out_sum).to_csv('mc_summary_by_barrier.csv', index=False); R = pd.DataFrame(out_rat); R.to_csv('mc_barrier_ratios.csv', index=False)
E = pd.DataFrame(out_el); E.to_csv('mc_elasticities.csv', index=False)
print(R.round(3).to_string(index=False)); print(E.round(3).to_string(index=False))
# ---------- design 1: all anatomies, assigned barriers
A = pd.read_csv('mc_all_x2.csv')
base = pd.read_csv(COND).rename(columns={'kr_cm_hPa_d': 'kr'})
dl = pd.read_csv(os.path.join(STATS, 'kr_summary_letters.csv')).set_index('fragment').tukey_letter.reindex(ORDER)
M = A.groupby(['draw', 'fragment']).kr.mean().unstack()[ORDER]
SM = pd.DataFrame({'default_mean': base.groupby('fragment').kr.mean().reindex(ORDER), 'mc_median_of_means': M.median(), 'mc_q025': M.quantile(.025), 'mc_q975': M.quantile(.975)})
SM['fold_95'] = SM.mc_q975 / SM.mc_q025
rel = M.div(M['LT'], axis=0); SM['rel_LT_default'] = SM.default_mean / SM.default_mean['LT']; SM['rel_LT_q025'] = rel.quantile(.025); SM['rel_LT_q975'] = rel.quantile(.975)
checks = {'n_draws': len(M), 'MP,MM < MT,LP,LM (every draw)': float(((M[['MP', 'MM']].max(axis=1)) < M[['MT', 'LP', 'LM']].min(axis=1)).mean()),
          'MT,LP,LM < MT0 < LT': float(((M[['MT', 'LP', 'LM']].max(axis=1) < M.MT0) & (M.MT0 < M.LT)).mean()),
          'MT > mean(LP,LM)': float((M.MT > M[['LP', 'LM']].mean(axis=1)).mean()),
          'main/lateral kr ratio (MM/LM) median, 95%': [float((M.MM / M.LM).median()), float((M.MM / M.LM).quantile(.025)), float((M.MM / M.LM).quantile(.975))]}
same = []; pats = {}
for dr, g in A.groupby('draw'):
    g = g.copy(); g['y'] = np.log10(g.kr); t = pairwise_tukeyhsd(g.y, g.fragment)
    tt = pd.DataFrame(t._results_table.data[1:], columns=t._results_table.data[0])
    L = cld(ORDER, {(r.group1, r.group2): bool(r.reject) for r in tt.itertuples()}); pat = ' '.join(L[o] for o in ORDER)
    pats[pat] = pats.get(pat, 0) + 1; same.append(all(L[o] == dl[o] for o in ORDER))
checks['Tukey letters identical to default'] = float(np.mean(same)); checks['letter patterns (count)'] = dict(sorted(pats.items(), key=lambda x: -x[1]))
checks['default letters'] = ' '.join(dl)
SM.to_csv('mc_segment_means.csv'); print(SM.to_string(float_format=lambda x: f'{x:.3g}'))
# cost-kr correlation per draw
COST = pd.read_csv(os.path.join(COSTS, 'replicate_costs.csv'))[['anatomy', 'section', 'segment', 'C_ng_mm', 'N_ng_mm']]
rows = []
for dr, g in A[A.fragment != 'MT0'].groupby('draw'):
    x = g.merge(COST, on='anatomy'); x['lc'] = np.log10(x.C_ng_mm); x['lk'] = np.log10(x.kr); x['ln'] = np.log10(x.N_ng_mm)
    o = dict(draw=dr, n=len(x), spearman_C=stats.spearmanr(x.C_ng_mm, x.kr)[0], spearman_N=stats.spearmanr(x.N_ng_mm, x.kr)[0])
    for c in ('lc', 'ln'):
        a = smf.ols(f'lk ~ C(section) + {c}', x).fit(); tv = a.tvalues[c]; o[f'rm_{c}'] = np.sign(tv) * np.sqrt(tv ** 2 / (tv ** 2 + a.df_resid)); o[f'rm_p_{c}'] = a.pvalues[c]
    rows.append(o)
CR = pd.DataFrame(rows); CR.to_csv('mc_cost_kr_correlation_per_draw.csv', index=False)
checks['cost-kr within-section r (building C) median, 95%, frac p<0.05'] = [float(CR.rm_lc.median()), float(CR.rm_lc.quantile(.025)), float(CR.rm_lc.quantile(.975)), float((CR.rm_p_lc < 0.05).mean())]
checks['cost-kr within-section r (N) median, 95%, frac p<0.05'] = [float(CR.rm_ln.median()), float(CR.rm_ln.quantile(.025)), float(CR.rm_ln.quantile(.975)), float((CR.rm_p_ln < 0.05).mean())]
checks['spearman C vs kr median, 95%'] = [float(CR.spearman_C.median()), float(CR.spearman_C.quantile(.025)), float(CR.spearman_C.quantile(.975))]
json.dump(checks, open('mc_rank_checks.json', 'w'), indent=1); print(json.dumps(checks, indent=1))
