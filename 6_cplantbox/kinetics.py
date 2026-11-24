"""Growth kinetics per architectural type from the PimpiBig2 day 1-5 tracings (roots keep their ID across days)."""
import cpb_paths as CP   # paths; changes into outputs/cplantbox/
import pandas as pd, numpy as np
pd.set_option('display.width',200)
D=pd.read_csv(CP.inp('timeseries_roots.csv'),dtype={'num':str})
p=pd.read_csv(CP.inp('01_panel_root_geometry.csv')); p=p[p.arbor_id.str.contains('Big2')].copy()
p['num']=p.arbor_id.str.split('_').str[6]; p['label']=p.arbor_id.str.split('_').str[7:].str.join('_')
D=D.merge(p[['num','label','architectural_type','condition']],on=['num','label'],how='inner')
D['plant']=D.num+'_'+D.label
M=D[D.order==0][['plant','day','length','architectural_type','condition']].rename(columns={'length':'MR'})
L=D[D.order==1].merge(M[['plant','day','MR']],on=['plant','day'])
# --- main root
mr=M.pivot_table(index=['architectural_type','condition','plant'],columns='day',values='MR')
print('MAIN ROOT LENGTH (cm), median'); print(mr.groupby(level=[0,1]).median().round(2))
rate=mr.diff(axis=1).iloc[:,1:]; print('MAIN ROOT GROWTH (cm/day), median'); print(rate.groupby(level=[0,1]).median().round(2))
# --- laterals: number, apical unbranched zone (MR length - deepest insertion)
g=L.groupby(['architectural_type','condition','plant','day'])
z=pd.DataFrame({'nLR':g.size(),'la':g.MR.first()-g.ins.max(),'tot':g.length.sum()})
print('LATERALS PER PLANT'); print(z.nLR.unstack().groupby(level=[0,1]).median())
print('APICAL UNBRANCHED ZONE (cm) = MR - deepest visible lateral'); print(z.la.unstack().groupby(level=[0,1]).median().round(2))
# --- each lateral: emergence day, length trajectory, elongation rate
first=L.groupby('id').day.min().rename('emerg'); L=L.join(first,on='id')
L['age']=L.day-L.emerg
s=L.sort_values('day').groupby('id')
L['dlen']=s.length.diff(); L['dday']=s.day.diff(); L['rate']=L.dlen/L.dday
r=L[L.dday==1]
print('LATERAL ELONGATION (cm/day), median over lateral-days, by lateral age (days since first seen)')
print(r.pivot_table(index=['architectural_type','condition'],columns='age',values='rate',aggfunc='median').round(2))
print('LATERAL ELONGATION by insertion third along the day-5 MR (control, cm/day)')
L5=L[L.day==5][['id','ins','MR']]; L5['third']=pd.cut(L5.ins/L5.MR,[0,1/3,2/3,1.01],labels=['top','mid','bottom'])
r2=r.merge(L5[['id','third']],on='id'); print(r2[r2.condition=='control'].pivot_table(index='architectural_type',columns='third',values='rate',aggfunc='median',observed=False).round(2))
# distance of new laterals from MR tip at the day they are first seen
new=L[(L.age==0)&(L.day>1)]; new=new.assign(from_tip=new.MR-new.ins)
print('NEW LATERALS: distance from MR tip when first seen (cm), median'); print(new.groupby(['architectural_type','condition']).from_tip.median().round(2))
print('n plants'); print(M[M.day==5].groupby(['architectural_type','condition']).size())
L.to_csv('timeseries_laterals_tracked.csv',index=False); M.to_csv('timeseries_mainroot.csv',index=False)
# ---- targets for the v2 fit: control medians per type and day, lateral rate by age (Big2 only)
import json
M=M[M.condition=='control']; L=L[L.condition=='control']; out={}
for t in ['B','C','dT']:
    m_=M[M.architectural_type==t]; l=L[L.architectural_type==t]
    g=l.groupby(['plant','day']); per=pd.DataFrame({'n':g.size(),'tot':g.length.sum()}).reset_index()
    per=m_[['plant','day']].merge(per,on=['plant','day'],how='left').fillna(0)
    r_=l[l.dday==1].groupby('age').rate.median()
    out[t]=dict(MR=m_.groupby('day').MR.median().round(3).tolist(),nLR=per.groupby('day').n.median().tolist(),totLR=per.groupby('day').tot.median().round(3).tolist(),
                rate_by_age=r_.round(3).to_dict(),n_plants=int(m_.plant.nunique()))
json.dump(out,open('kinetics_targets.json','w'),indent=1)
