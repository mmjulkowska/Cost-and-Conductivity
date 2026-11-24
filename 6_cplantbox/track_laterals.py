"""Follow laterals from day to day by insertion position (SmartRoot IDs are not kept across all days in Big5/Big6).
A lateral on day d+1 is the same root as the lateral on day d of the same plant whose insertion lies within 0.15 cm and that is not longer
than it (+0.05 cm); nearest first. Checked against Big2 RSML IDs (days 4->5). Writes pooled_lateral_rates.csv."""
import cpb_paths as CP   # paths; changes into outputs/cplantbox/
import pandas as pd, numpy as np
L=pd.read_csv(CP.inp('pooled_laterals.csv')); P=pd.read_csv(CP.inp('pooled_plants_by_day.csv'),dtype={'date':str})
L=L.merge(P[['experiment','label','day','date']],on=['experiment','label','day'])
L['t']=pd.to_datetime(L.date.astype(str)).map(pd.Timestamp.toordinal)
def match(a,b,tol=0.15):
    pairs=[]; used=set()
    cand=sorted(((abs(x.insertion_position-y.insertion_position),i,j) for i,x in a.iterrows() for j,y in b.iterrows()
                 if abs(x.insertion_position-y.insertion_position)<tol and y.length>=x.length-0.05))
    ua=set()
    for dd,i,j in cand:
        if i in ua or j in used: continue
        ua.add(i); used.add(j); pairs.append((i,j))
    return pairs
rows=[]; check=[]
for (ex,lab),g in L.groupby(['experiment','label']):
    ts=sorted(g.t.unique()); age={}
    for k,t in enumerate(ts):
        b=g[g.t==t]
        if k==0 or t-ts[k-1]!=1:
            for j in b.index: age[j]=None if k==0 else None
            continue
        a=g[g.t==ts[k-1]]; pr=match(a,b); matched=set(j for _,j in pr)
        for j in b.index:
            if j not in matched: age[j]=0        # first seen today
        for i,j in pr:
            age[j]=None if age.get(i) is None else age[i]+1
            if age[j] is not None:
                rows.append(dict(experiment=ex,label=lab,type=g.architectural_type.iloc[0],condition=g.condition.iloc[0],age=age[j],
                                 rate=b.at[j,'length']-a.at[i,'length'],ins=b.at[j,'insertion_position'],MR=b.at[j,'MR']))
            check.append(a.at[i,'root']==b.at[j,'root'])
R=pd.DataFrame(rows); R.to_csv('pooled_lateral_rates.csv',index=False)
print('matches where SmartRoot ID also agrees (only meaningful where IDs persist):', np.mean(check).round(3), len(check))
c=R[R.condition=='control']
print(c.pivot_table(index=['type'],columns='age',values='rate',aggfunc='median').round(2))
print(c.pivot_table(index=['type','experiment'],columns='age',values='rate',aggfunc='size'))
