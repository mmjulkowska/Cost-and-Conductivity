"""Salt-induced Krs loss against collection site and genetic group.

Joins the accession values (krs_accession_values.csv, from krs_heritability.py) to
  - S.pimpinellifolium GPS data.xlsx (latitude/longitude, 227 accessions; Dropbox Tomato_Pimpinellifolium/Pimp_SNPs)
  - pcs_with_krs_traits.tsv (5 PCs and the k-means group of each genotyped accession, written on the server by
    gwas/krs_vs_structure.R; groups other than the largest are the 'minority' groups). Optional.
Tests latitude, longitude and Ecuador (north of 3.5 S) vs Peru. Writes krs_accession_geo.csv; prints the tests.
"""
import krs_paths as KP
import os, pandas as pd, numpy as np
from scipy.stats import spearmanr, mannwhitneyu
gps = pd.read_excel(KP.inp('S.pimpinellifolium GPS data.xlsx')); gps['Accession'] = gps.Accession.str.strip()
geo = pd.read_csv(KP.inp('01_panel_root_geometry.csv'), dtype=str)[['accession_code', 'tgrc_code']].drop_duplicates()
grp_file = os.path.join(KP.KRS_IN, 'pcs_with_krs_traits.tsv')
KP.to_out()
A = pd.read_csv('krs_accession_values.csv', dtype={'accession_code': str})
geo['ac'] = geo.accession_code.astype(int).astype(str); A['ac'] = A.accession_code.astype(int).astype(str)
A = A.merge(geo[['ac', 'tgrc_code']], on='ac').merge(gps, left_on='tgrc_code', right_on='Accession', how='left')
if os.path.exists(grp_file):
    g = pd.read_csv(grp_file, sep='\t')[['IID', 'group']]
    big = g.group.value_counts().idxmax()
    A = A.merge(g.rename(columns={'IID': 'tgrc_code'}), on='tgrc_code', how='left'); A['minor'] = A.group.notna() & (A.group != big)
else:
    A['group'] = np.nan; A['minor'] = False
A['ecuador'] = A.Latitude > -3.5
print('accessions with coordinates:', A.Latitude.notna().sum(), 'of', len(A))
print(pd.crosstab(A.minor, A.ecuador, rownames=['minority group'], colnames=['north of 3.5 S']))
for t in ['logKrs.C', 'logKrs.S', 'logKrs.SminusC', 'KrsResid.C', 'KrsResid.S', 'KrsResid.SminusC']:
    d = A.dropna(subset=[t, 'Latitude'])
    r, p = spearmanr(d.Latitude, d[t]); r2, p2 = spearmanr(d.Longitude, d[t])
    e = d[d.ecuador][t]; q = d[~d.ecuador][t]
    print(f'{t:18s} n={len(d)} lat rho={r:+.2f} P={p:.2g} | lon rho={r2:+.2f} P={p2:.2g} | '
          f'Ecuador {e.median():+.3f} (n={len(e)}) vs Peru {q.median():+.3f} (n={len(q)}) P={mannwhitneyu(e, q).pvalue:.2g}')
A.to_csv('krs_accession_geo.csv', index=False)
