"""Bracket lateral kr (Nov 2026): LP, LM with current B10 (exodermis only) vs B6 (exodermis + endodermal Casparian strip).
Endodermal strip in laterals was 'not detected' with toluidine blue brightfield, which cannot resolve a CS. Default hydraulics, all v5 anatomies."""
import os, sys; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from paths import *
import os, sys, numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lignin_cap import run, HERE
def job(t):
    a, s = t; out = []
    for tr in ('B10', 'B6'):
        try: kr = run(a, tr)
        except Exception: kr = np.nan
        out.append(dict(anatomy=a, fragment=s, treatment=tr, kr=kr))
    return out
if __name__ == '__main__':
    C = pd.read_csv(COND); C = C[C.kr_cm_hPa_d.notna()]
    L = C[C.fragment.isin(['LP', 'LM'])]
    with Pool(max(1, (os.cpu_count() or 2) - 1)) as p: R = [r for res in p.map(job, list(zip(L.anatomy, L.fragment))) for r in res]
    D = pd.DataFrame(R); D.to_csv(os.path.join(HERE, 'bracket_laterals_kr.csv'), index=False)
    W = D.pivot_table(index=['fragment', 'anatomy'], columns='treatment', values='kr').reset_index(); W['ratio'] = W.B6 / W.B10
    S = W.groupby('fragment').agg(n=('ratio', 'count'), kr_B10_mean=('B10', 'mean'), kr_B10_sd=('B10', 'std'), kr_B6_mean=('B6', 'mean'), kr_B6_sd=('B6', 'std'),
                                 ratio_median=('ratio', 'median'), ratio_min=('ratio', 'min'), ratio_max=('ratio', 'max'))
    ref = C.groupby('fragment').kr_cm_hPa_d.agg(['mean', 'std', 'count']); S.to_csv(os.path.join(HERE, 'bracket_laterals_summary.csv'))
    print(S.to_string(float_format=lambda x: f'{x:.3g}')); print(ref.to_string(float_format=lambda x: f'{x:.3g}'))
    # how does ranking change: all segments with laterals as B6
    A = C[['anatomy', 'fragment', 'kr_cm_hPa_d']].rename(columns={'kr_cm_hPa_d': 'kr'})
    B = A.copy(); m = W.set_index('anatomy').B6; B.loc[B.fragment.isin(['LP', 'LM']), 'kr'] = B.anatomy.map(m)
    import statsmodels.api as sm; from statsmodels.stats.multicomp import pairwise_tukeyhsd
    for nm, X in (('current', A), ('laterals_B6', B)):
        X = X.dropna(); t = pairwise_tukeyhsd(np.log10(X.kr), X.fragment)
        r = pd.DataFrame(t.summary().data[1:], columns=t.summary().data[0]); r.to_csv(os.path.join(HERE, f'bracket_tukey_{nm}.csv'), index=False)
        print(nm); print(r[r.group1.isin(['LP','LM','MP','MM','MT']) & r.group2.isin(['LP','LM','MP','MM','MT'])].to_string(index=False))
