"""Pooled targets (control) per type for the v3 fit.
Time courses: Big2 (RSML, days 1-5), Big5 (days by date: 20 Oct = 2, 22 Oct = 4, 23 Oct = 5; 19 Oct dropped: main roots shorter on the next day,
so the day-1 export is not the same tracing series), Big6 (days 1-5). Each experiment's median trajectory is divided by its day-5 value,
the experiments are averaged, and the result is scaled to the panel day-5 median (all six experiments).
Lateral position profiles: every experiment at its panel imaging day. Lateral rate by age: Big2 (RSML IDs) + Big5/6 (position-tracked)."""
import cpb_paths as CP   # paths; changes into outputs/cplantbox/
import pandas as pd, numpy as np, json
P=pd.read_csv(CP.inp('pooled_plants_by_day.csv'),dtype={'date':str}); L=pd.read_csv(CP.inp('pooled_laterals.csv')); R=pd.read_csv('pooled_lateral_rates.csv')
K2=json.load(open('kinetics_targets.json')); pan=pd.read_csv(CP.inp('01_panel_root_geometry.csv')); pan=pan[pan.condition=='control']
B5={'20191020':2,'20191022':4,'20191023':5}
P['eday']=np.where(P.experiment=='PimpiBig5',P.date.map(B5),P.day)
BINS=np.arange(0,1.0001,0.1); out={}
for t in ['B','C','dT']:
    traj=[]
    k=K2[t]; traj.append(pd.DataFrame({'MR':k['MR'],'nLR':k['nLR'],'tot':k['totLR']},index=range(1,6)))
    for ex in ['PimpiBig5','PimpiBig6']:
        c=P[(P.experiment==ex)&(P.architectural_type==t)&(P.condition=='control')&P.eday.notna()]
        traj.append(c.groupby('eday')[['MR','nLR','totLR']].median().rename(columns={'totLR':'tot'}).reindex(range(1,6)))
    norm=[x/x.loc[5] for x in traj]; avg=pd.concat(norm).groupby(level=0).mean()
    g=pan[pan.architectural_type==t]; anchor={'MR':g.MR_length.median(),'nLR':g.n_lateral_roots.median(),'tot':g.total_lateral_root_length.median()}
    tg={k_:(avg[k_]*anchor[k_]).round(3).tolist() for k_ in ['MR','nLR','tot']}
    # profiles at panel day
    l=L[(L.architectural_type==t)&(L.condition=='control')&(L.day==L.imaging_day)].copy(); l['rel']=(l.insertion_position/l.MR).clip(0,0.9999)
    cnt=np.mean([np.histogram(x.rel,BINS)[0]/len(x) for _,x in l.groupby(['experiment','label'])],axis=0)
    prof=np.mean([np.histogram(x.rel,BINS,weights=x.length)[0]/x.length.sum() for _,x in l.groupby(['experiment','label'])],axis=0)
    r=R[(R.type==t)&(R.condition=='control')&(R.experiment.isin(['PimpiBig5','PimpiBig6']))].groupby('age').rate.agg(['median','size'])
    rate=[]
    for a in (1,2,3,4):
        v=[(k['rate_by_age'][str(a)],k['n_plants'])]
        if a in r.index: v.append((r.loc[a,'median'],r.loc[a,'size']/5))
        rate.append(round(sum(x*w for x,w in v)/sum(w for _,w in v),3))
    out[t]=dict(**tg,cnt=cnt.round(4).tolist(),prof=prof.round(4).tolist(),rate=rate,n_profile_plants=int(l.groupby(['experiment','label']).ngroups),
                n_timecourse=dict(Big2=k['n_plants'],**{ex:int(P[(P.experiment==ex)&(P.architectural_type==t)&(P.condition=='control')].label.nunique()) for ex in ['PimpiBig5','PimpiBig6']}),
                ins_angle_median=float(l.insertion_angle.median()))
    print(t,json.dumps({k_:out[t][k_] for k_ in ['MR','nLR','tot','rate','n_profile_plants','n_timecourse']}))
json.dump(out,open('kinetics_targets_v3.json','w'),indent=1)
