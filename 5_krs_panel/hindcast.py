"""Hindcast check of krs_project.py: plant parameters from tracing days <= 3 only (the panel rate-by-age curve
from all days), predict days 4-5 and compare with the traced values. Writes hindcast_<LNMODE>.csv."""
import os, sys, pandas as pd, numpy as np
import krs_project as K   # changes into outputs/krs/
obs = pd.read_csv('krs_per_plant_day.csv')
pl = set(K.P[K.P.experiment.isin(['PimpiBig2','PimpiBig5','PimpiBig6'])].groupby(K.k3).filter(lambda g: (g.day<=3).sum()>=2 and (g.day==5).any())[K.k3].drop_duplicates().itertuples(index=False,name=None))
print(len(pl),'plants with >=2 days up to day 3 and a day-5 tracing')
C5 = K.build_curves(K.INC)   # panel-level rate-by-age curve from all days; plant-level parameters from days <= 3 only
R,Pp,C = K.run(3,5,pl,C5)
m = R[R.day.isin([4,5])].merge(obs[['experiment','label','day','MR','nLR','TRL','Krs']], on=['experiment','label','day'], suffixes=('_pred','_obs'))
m.to_csv(f"hindcast_{os.environ.get('LNMODE','kinetic')}.csv", index=False)
def s(x):
    out={}
    for v in ['MR','nLR','TRL','Krs']:
        p=x[v+'_pred']; o=x[v+'_obs']
        out[v+'_ratio_med']=np.median(p/o.clip(lower=1e-9)) if v!='nLR' else np.median((p+1)/(o+1))
        out[v+'_r']=np.corrcoef(np.log(p+1e-9 if v!='nLR' else p+1),np.log(o+1e-9 if v!='nLR' else o+1))[0,1]
    return pd.Series(out)
print(m.groupby(['rule','cond','day']).apply(s).round(2).to_string())
