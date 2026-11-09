"""Lignin-cap check (Nov 2026): exodermal barrier as a lignin cap (MECHA B9: outer tangential exodermis wall + exodermal radial band)
instead of full exodermal suberin + lignin (B6/B10). Main-root MP, MM keep the endodermal band (B9 as is); laterals LP, LM get
B9 with the endodermal band removed ('B9noEndo'). Default hydraulics; all v5 anatomies of these segments."""
import os, sys; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from paths import *
import os, io, sys, contextlib, numpy as np, pandas as pd
from multiprocessing import Pool
from openalea.mecha import Mecha, InData
from openalea.mecha.utils import data_loader as DL
HERE = CHECKS
cls = [c for c in vars(DL).values() if isinstance(c, type) and hasattr(c, 'get_wall_conductivities')][0]
_orig = cls.get_wall_conductivities
NOENDO = {'on': False}
def patched(self, barrier, h):
    c = _orig(self, barrier, h)
    if barrier == 9 and NOENDO['on']: c['kw_endo_endo'] = self.get_kw_value(h)
    return c
cls.get_wall_conductivities = patched
def run(a, treat):
    b = {'B1': 1, 'B6': 6, 'B10': 10, 'B9': 9, 'B9noEndo': 9}[treat]; NOENDO['on'] = treat == 'B9noEndo'
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        d = InData(cellset_file=os.path.join(ANAT, a + '_root.xml')); d.geometry.im_scale = 1000.0
        d.geometry.set_maturity_stages([b], [200.0]); m = Mecha(d); m.compute_conductivities()
    return m.root_hydraulic_properties[-1]['kr']
TR = {'MP': ['B6', 'B9', 'B1'], 'MM': ['B6', 'B9', 'B1'], 'LP': ['B10', 'B9noEndo', 'B9'], 'LM': ['B10', 'B9noEndo', 'B9']}
def job(t):
    a, s = t; out = []
    for tr in TR[s]:
        try: kr = run(a, tr)
        except Exception: kr = np.nan
        out.append(dict(anatomy=a, fragment=s, treatment=tr, kr=kr))
    return out
if __name__ == '__main__':
    C = pd.read_csv(COND)
    C = C[C.fragment.isin(TR) & C.kr_cm_hPa_d.notna()]
    with Pool(max(1, (os.cpu_count() or 2) - 1)) as p: R = [r for res in p.map(job, list(zip(C.anatomy, C.fragment))) for r in res]
    D = pd.DataFrame(R); D.to_csv(os.path.join(HERE, 'lignin_cap_kr.csv'), index=False)
    W = D.pivot_table(index=['fragment', 'anatomy'], columns='treatment', values='kr').reset_index()
    rows = []
    for s, g in W.groupby('fragment'):
        ref = 'B6' if s in ('MP', 'MM') else 'B10'
        for tr in TR[s][1:]:
            r = g[tr] / g[ref]; rows.append(dict(fragment=s, reference=ref, treatment=tr, n=r.notna().sum(), kr_ref_mean=g[ref].mean(), kr_mean=g[tr].mean(),
                                                ratio_median=r.median(), ratio_min=r.min(), ratio_max=r.max()))
    S = pd.DataFrame(rows).set_index('fragment').loc[['MP', 'MM', 'LP', 'LM']]; S.to_csv(os.path.join(HERE, 'lignin_cap_summary.csv'))
    print(S.to_string(float_format=lambda x: f'{x:.3g}'))
