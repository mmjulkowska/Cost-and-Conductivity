"""Heritability of Krs-derived traits (tracing day 5) and accession-level values for GWAS.

Traits per plant (day 5; d4->d5 for the rate):
  logKrs      ln Krs (cm3 hPa-1 d-1)
  logTRL      ln total root length (reference: how much of Krs is just size)
  KrsPerTRL   Krs / TRL  (x1e5)
  KrsResid    residual of ln Krs on ln TRL (pooled regression, both conditions): conductance beyond size
  KrsRGR      ln(Krs d5 / Krs d4), relative daily gain in Krs
Model per condition: y ~ experiment (fixed) + (1 | accession)  (REML, statsmodels MixedLM)
  H2_plant = Vg / (Vg + Ve);  H2_mean = Vg / (Vg + Ve / n_h), n_h = harmonic mean plants per accession
  95% CI: 300 bootstrap resamples of accessions.
G x E: y ~ condition * experiment + (1|accession) + (1|accession:condition), both conditions.
Accession values: y minus the fitted experiment effect, averaged per accession (control, salt) and
salt response = salt - control (on ln scale for logKrs = ln of the salt/control ratio).
"""
import pandas as pd, numpy as np, warnings, os, statsmodels.formula.api as smf
import krs_paths as KP
KP.to_out()
warnings.filterwarnings('ignore')
rng = np.random.default_rng(1)
R = pd.read_csv('krs_per_plant_day.csv', dtype={'accession_code': str})
k = ['experiment', 'label']
d5 = R[R.day == 5].set_index(k); d4 = R[R.day == 4].set_index(k)
D = d5.reset_index().copy()
D['logKrs'] = np.log(D.Krs); D['logTRL'] = np.log(D.TRL); D['KrsPerTRL'] = D.Krs_per_TRL * 1e5
b = np.polyfit(D.logTRL, D.logKrs, 1); D['KrsResid'] = D.logKrs - np.polyval(b, D.logTRL)
D['KrsRGR'] = np.log(d5.Krs / d4.Krs.reindex(d5.index)).values
D['cond'] = D.condition.map({'control': 'C', 'NaCl_100mM': 'S'}); D['acc'] = D.accession_code
TRAITS = ['logKrs', 'logTRL', 'KrsPerTRL', 'KrsResid', 'KrsRGR']


def fit(d, y):
    d = d.dropna(subset=[y])
    m = smf.mixedlm(f'{y} ~ C(experiment)', d, groups=d['acc']).fit(reml=True)
    vg = float(m.cov_re.iloc[0, 0]); ve = float(m.scale)
    n = d.groupby('acc').size(); nh = len(n) / np.sum(1 / n)
    return vg, ve, nh, m


def h2(d, y):
    vg, ve, nh, m = fit(d, y)
    return vg / (vg + ve), vg / (vg + ve / nh), vg, ve, nh, m


rows = []; vals = {}
for y in TRAITS:
    for c in 'CS':
        d = D[D.cond == c]
        hp, hm, vg, ve, nh, m = h2(d, y)
        accs = d.acc.unique(); bs = []
        for _ in range(int(os.environ.get('NBOOT', 300))):
            s = rng.choice(accs, len(accs)); parts = []
            for i, a in enumerate(s):
                x = d[d.acc == a].copy(); x['acc'] = f'{a}_{i}'; parts.append(x)
            try: bs.append(h2(pd.concat(parts), y)[:2])
            except Exception: pass
        bs = np.array(bs)
        rows.append(dict(trait=y, condition=c, n_plants=int(d[y].notna().sum()), n_acc=d.loc[d[y].notna(), 'acc'].nunique(),
                         plants_per_acc=round(nh, 2), Vg=vg, Ve=ve, H2_plant=hp, H2_plant_lo=np.quantile(bs[:, 0], .025),
                         H2_plant_hi=np.quantile(bs[:, 0], .975), H2_mean=hm, H2_mean_lo=np.quantile(bs[:, 1], .025),
                         H2_mean_hi=np.quantile(bs[:, 1], .975)))
        fe = m.fe_params; ex = {e: fe.get(f'C(experiment)[T.{e}]', 0.) for e in d.experiment.unique()}
        adj = d[y] - d.experiment.map(ex)
        vals[f'{y}.{c}'] = adj.groupby(d.acc).mean()
H = pd.DataFrame(rows); H.to_csv('krs_heritability.csv', index=False)
print(H[['trait', 'condition', 'n_acc', 'plants_per_acc', 'H2_plant', 'H2_plant_lo', 'H2_plant_hi', 'H2_mean', 'H2_mean_lo', 'H2_mean_hi']].round(2).to_string())

# G x E
g = []
for y in TRAITS:
    d = D.dropna(subset=[y]).copy(); d['one'] = 1
    m = smf.mixedlm(f'{y} ~ C(cond) * C(experiment)', d, groups=d['acc'], re_formula='1',
                    vc_formula={'gxe': '0 + C(cond)'}).fit(reml=True)
    vg = float(m.cov_re.iloc[0, 0]); vgxe = float(m.vcomp[0]); ve = float(m.scale); tot = vg + vgxe + ve
    g.append(dict(trait=y, Vg_share=vg / tot, VgxE_share=vgxe / tot, Ve_share=ve / tot,
                  salt_effect=m.fe_params['C(cond)[T.S]'], salt_p=m.pvalues['C(cond)[T.S]']))
G = pd.DataFrame(g); G.to_csv('krs_gxe.csv', index=False); print(G.round(3).to_string())

A = pd.DataFrame(vals)
for y in TRAITS: A[f'{y}.SminusC'] = A[f'{y}.S'] - A[f'{y}.C']
A.index.name = 'accession_code'; A.to_csv('krs_accession_values.csv')
print(A.shape, A.notna().sum().to_dict())
print('corr of accession values with Arjun panel types etc. skipped')
print('C vs S correlation, logKrs:', A[['logKrs.C', 'logKrs.S']].corr().iloc[0, 1].round(2), ' KrsResid:', A[['KrsResid.C', 'KrsResid.S']].corr().iloc[0, 1].round(2))
print('pooled ln Krs ~ ln TRL slope/intercept', b)
