import cpb_paths as CP   # paths; changes into outputs/cplantbox/
import pandas as pd, numpy as np
d=pd.read_csv(CP.inp('01_panel_root_geometry.csv'))
d=d[d.condition=='control']
L=pd.read_csv('rsml_laterals.csv'); L=L[L.cond=='control']
BINS=np.arange(0,1.0001,0.1)
def profile(rel,length):
    h,_=np.histogram(np.clip(rel,0,0.9999),BINS,weights=length); return h/h.sum()
T={}
for t in ['B','C','dT']:
    g=d[d.architectural_type==t]; l=L[L.type==t]
    prof=np.mean([profile(x.rel,x.len) for _,x in l.groupby('arbor')],axis=0)
    last=l.groupby('arbor').rel.max().median()
    T[t]=dict(MR=g.MR_length.median(),nLR=g.n_lateral_roots.median(),totLR=g.total_lateral_root_length.median(),
              ins=g.LR_insertion_angle_deg.median(),vert=g.LR_angle_from_vertical_deg.median(),last=last,prof=prof,
              n=len(g),n_rsml=l.arbor.nunique(),
              cvMR=g.MR_length.std()/g.MR_length.mean(),cvN=g.n_lateral_roots.std()/g.n_lateral_roots.mean())
