"""Inflation check for the salt Krs traits: raw accession means (no prior experiment adjustment) with experiment as
EMMAX covariates, next to the 5 PCs. Control ln Krs is included as a reference (lambda was already ok).
Writes outputs/krs/krs_gwas_expcov/ (copy this folder to the server): pheno/emmax/trait.<name>.INT.pheno, exp_covariates.txt (265 rows, tfam order, share of the
accession's day-5 plants in each experiment, PimpiBig1 = reference, 0 for unphenotyped samples), run_expcov.sh."""
import pandas as pd, numpy as np, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import krs_paths as KP
from scipy.stats import norm
KEY = KP.inp('sample_key.csv'); GEO = KP.inp('01_panel_root_geometry.csv')
KP.to_out()
R = pd.read_csv('krs_per_plant_day.csv', dtype={'accession_code': str}); D = R[R.day == 5].copy()
D['code'] = D.accession_code.str.zfill(3); D['cond'] = D.condition.map({'control': 'C', 'NaCl_100mM': 'S'})
D['y'] = np.log(D.Krs)
M = D.groupby(['code', 'cond']).y.mean().unstack()
T = {'logKrsRaw.C': M['C'], 'logKrsRaw.S': M['S'], 'logKrsRaw.SminusC': M['S'] - M['C']}
key = pd.read_csv(KEY, dtype=str)
geo = pd.read_csv(GEO, dtype=str)[['accession_code', 'tgrc_code']].drop_duplicates()
geo['code'] = geo.accession_code.str.zfill(3); la2code = dict(zip(geo.tgrc_code, geo.code))
code = key.IID.map(la2code).where(key.is_replicate != 'True')
out = 'krs_gwas_expcov'; os.makedirs(f'{out}/pheno/emmax', exist_ok=True)


def INT(x):
    x = x.copy(); ok = x.notna(); r = x[ok].rank(method='average'); x[ok] = norm.ppf((r - 3 / 8) / (ok.sum() + 1 / 4)); return x


for t, v in T.items():
    vv = code.map(v)
    pd.DataFrame({'FID': key.FID, 'IID': key.IID, 'v': INT(vv)}).to_csv(f'{out}/pheno/emmax/trait.{t}.INT.pheno', sep='\t', header=False, index=False, na_rep='NA', float_format='%.10g')
    print(t, vv.notna().sum())
share = pd.crosstab(D.code, D.experiment, normalize='index')
ex = [e for e in sorted(share.columns) if e != 'PimpiBig1']
cov = share.reindex(code.values)[ex].fillna(0.).round(4)
cov.to_csv(f'{out}/exp_covariates.txt', sep='\t', header=False, index=False)
import shutil; shutil.copy(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'run_expcov.sh'), out)
print(cov.shape, 'experiment columns:', ex); print((cov > 0).sum().to_dict())
