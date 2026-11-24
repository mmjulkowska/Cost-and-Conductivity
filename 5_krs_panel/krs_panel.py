"""Krs for every traced S. pimpinellifolium seedling and tracing day (PimpiBig1-6, control and 100 mM NaCl).

Inputs (KRS_IN, see krs_paths.py): pooled_plants_by_day.csv, pooled_laterals.csv (SmartRoot exports, laterals >= 0.1 cm),
timeseries_roots.csv (PimpiBig2 RSML, days 1-3 added from here; laterals filtered at >= 0.1 cm to match the exports).
PimpiBig5 tracing days follow the image dates (as in the v3 CPlantBox fit): 20 Oct = d2, 21 = d3, 22 = d4, 23 = d5;
19 Oct is dropped (main roots longer than on the next day).
Outputs (outputs/krs/): krs_per_plant_day.csv, plants_assembled.csv, laterals_assembled.csv
"""
import pandas as pd, numpy as np, krs_tree as T, sys
import krs_paths as KP
P = pd.read_csv(KP.inp('pooled_plants_by_day.csv')); L = pd.read_csv(KP.inp('pooled_laterals.csv'))
KP.to_out()
P['date'] = P.date.astype(str)

# PimpiBig5 day remap by date
b5 = {'20191020': 2, '20191021': 3, '20191022': 4, '20191023': 5}
L = L.merge(P[['experiment', 'label', 'day', 'date']], on=['experiment', 'label', 'day'], how='left')
P = P[~((P.experiment == 'PimpiBig5') & (P.date == '20191019'))].copy()
L = L[~((L.experiment == 'PimpiBig5') & (L.date == '20191019'))].copy()
for X in (P, L):
    m = X.experiment == 'PimpiBig5'; X.loc[m, 'day'] = X.loc[m, 'date'].map(b5)
assert not P.duplicated(['experiment', 'label', 'day']).any()

# PimpiBig2 days 1-3 from the RSML time series
TS = pd.read_csv(KP.inp('timeseries_roots.csv')); TS = TS[TS.day <= 3]
mr = TS[TS.order == 0].sort_values('length').drop_duplicates(['label', 'day'], keep='last')
lat = TS[(TS.order == 1) & (TS.length >= 0.1)].merge(mr[['day', 'id', 'length']].rename(columns={'id': 'parent', 'length': 'MR'}), on=['day', 'parent'])
meta = P[P.experiment == 'PimpiBig2'][['label', 'architectural_type', 'condition', 'imaging_day']].drop_duplicates('label')
g = lat.groupby(['label', 'day']); agg = pd.DataFrame({'nLR': g.size(), 'totLR': g.length.sum()}).reset_index()
p2 = mr[['label', 'day', 'length']].rename(columns={'length': 'MR'}).merge(agg, on=['label', 'day'], how='left').fillna({'nLR': 0, 'totLR': 0})
p2 = p2.merge(meta, on='label'); p2['experiment'] = 'PimpiBig2'
P = pd.concat([P, p2], ignore_index=True)
L = pd.concat([L[['experiment', 'label', 'day', 'length', 'insertion_position']],
               lat.assign(experiment='PimpiBig2').rename(columns={'ins': 'insertion_position'})[['experiment', 'label', 'day', 'length', 'insertion_position']]],
              ignore_index=True)

P.to_csv('plants_assembled.csv', index=False); L.to_csv('laterals_assembled.csv', index=False)
LG = {k: list(zip(v.insertion_position, v.length)) for k, v in L.groupby(['experiment', 'label', 'day'])}
rows = []
for r in P.itertuples():
    lats = LG.get((r.experiment, r.label, r.day), [])
    k, up, ln = T.krs(r.MR, lats)
    k_endo, _, _ = T.krs(r.MR, lats, 'lat_endo')
    k_t05, _, _ = T.krs(r.MR, lats, TIP=0.5)
    k_t3, _, _ = T.krs(r.MR, lats, TIP=3.0)
    trl = r.MR + sum(l for _, l in lats)
    d = dict(experiment=r.experiment, label=r.label, day=int(r.day), condition=r.condition,
             architectural_type=r.architectural_type, panel_day=r.imaging_day, MR=r.MR, nLR=len(lats), TRL=trl,
             Krs=k, Krs_lat_endo=k_endo, Krs_tip0p5=k_t05, Krs_tip3=k_t3, Krs_per_TRL=k / trl)
    d.update({f'up_{c}': up[c] / k for c in T.CLASSES})
    rows.append(d)
R = pd.DataFrame(rows)
R['accession_code'] = R.label.str.split('_').str[0]; R['replicate'] = R.label.str.split('_').str[1].astype(int)
R.to_csv('krs_per_plant_day.csv', index=False)
print(R.groupby(['experiment', 'day']).size().unstack())
print(R.groupby(['condition', 'day']).Krs.median().unstack() * 1e4)
