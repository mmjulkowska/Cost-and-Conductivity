"""Krs over time and at equal size (30 cm total root length) for the v4 architectures, 12 replicates per topology (3D).
Variants: v5 (laterals without endodermal barrier, B10), lat_endo (laterals with one, B6), and tip-zone length 0.5 / 1.5 / 3 cm."""
import cpb_paths as CP   # paths; changes into outputs/cplantbox/
import numpy as np, pandas as pd, json, os, plantbox as pb, rsa_model as m, conductance as K
P=json.load(open(os.environ.get('PARAMS','params_v4.json'))); OUT=os.environ.get('KOUT','krs'); os.makedirs(OUT,exist_ok=True)
END=12.; NREP=12; rows=[]; at30=[]
for t in ['C','B','dT']:
    t0=P[t]['T_day5']-5
    for s in range(NREP):
        p=m.build(P[t]['3D'],2000+s,stochastic=True); tt=0.; done30=False
        marks=np.round(np.arange(0.5,t0+END+1e-9,0.5),3)
        for mk in marks:
            while tt<mk-1e-9:
                trl=sum(o.getLength() for o in p.getOrgans(pb.root)) if False else None
                st=min(0.1,mk-tt)
                if not done30:
                    L=sum(o.getLength() for o in p.getOrgans(m.pb.root))
                    if L>27: st=min(st,0.01)
                    if L>=30:
                        segs=K.segments(p); k,up,_=K.krs(p,segs=segs); k2,_,_=K.krs(p,'lat_endo',segs=segs)
                        at30.append(dict(type=t,rep=s+1,age=tt,day=tt-t0,TRL=sum(x[2] for x in segs),Krs=k,Krs_lat_endo=k2,**{f'up_{c}':up[c]/k for c in K.CLASSES},
                                         **{f'len_{c}':sum(x[2] for x in segs if x[3]==c) for c in K.CLASSES}))
                        done30=True
                p.simulate(st,False); tt+=st
            for tip in (0.5,1.5,3.0):
                segs=K.segments(p,TIP=tip)
                for sc in ('v5','lat_endo'):
                    if tip!=1.5 and sc!='v5': continue
                    k,up,_=K.krs(p,sc,segs=segs)
                    rows.append(dict(type=t,rep=s+1,age=mk,day=mk-t0,tip=tip,scenario=sc,TRL=sum(x[2] for x in segs),Krs=k,**{f'up_{c}':(up[c]/k if k else 0) for c in K.CLASSES},
                                     **{f'len_{c}':sum(x[2] for x in segs if x[3]==c) for c in K.CLASSES}))
    print(t,'done')
R=pd.DataFrame(rows); R.to_csv(f'{OUT}/krs_timecourse.csv',index=False); A=pd.DataFrame(at30); A.to_csv(f'{OUT}/krs_at_30cm.csv',index=False)
pd.set_option('display.width',250)
print(A.groupby('type')[['day','TRL','Krs','Krs_lat_endo']+[f'up_{c}' for c in K.CLASSES]].median().round(4).to_string())
q=R[(R.tip==1.5)&(R.scenario=='v5')]; print(q[q.day.isin([3,5,8,12])].groupby(['type','day'])[['TRL','Krs']].median().unstack().round(4).to_string())
