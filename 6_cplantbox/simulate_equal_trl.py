"""Grow stochastic replicates of each topology until total root length (main root + laterals) reaches the same target."""
import cpb_paths as CP   # paths; changes into outputs/cplantbox/
import numpy as np, json, os, pandas as pd, plantbox as pb, rsa_model as m
TARGET=float(os.environ.get('TRL',30.)); NREP=int(os.environ.get('NREP',12)); OUT=os.environ.get('OUT','sims')
PAR=json.load(open(os.environ.get('PARAMS','params_fitted.json'))); os.makedirs(OUT,exist_ok=True)
NAME={'C':'Christmas tree','B':'Broomstick','dT':'Droopy telephone pole'}
def trl(p): return sum(o.getLength() for o in p.getOrgans(pb.root))
def grow_to(P,seed,plate):
    p=m.build(P,seed,stochastic=True,plate=plate); t=0.; dt=0.1
    while True:
        L0=trl(p)
        if L0>=TARGET-1e-3 or t>40: break
        step=dt if TARGET-L0>1.0 else 0.01
        p.simulate(step,False); t+=step
    return p,t
rows=[]
for geo in ['plate','3D']:
    for t in ['C','B','dT']:
        for s in range(NREP):
            p,T=grow_to(PAR[t][geo],1000+s,geo=='plate'); o=m.measure(p); L=o['lats']
            segs=np.array([[n.x,n.y,n.z] for n in p.getNodes()])
            fn=f'{OUT}/{geo}_{t}_rep{s+1:02d}'
            p.write(fn+'.rsml'); p.write(fn+'.vtp')
            np.save(fn+'_polylines.npy',np.array([np.array([[n.x,n.y,n.z] for n in r.getNodes()]) for r in p.getOrgans(pb.root)],dtype=object),allow_pickle=True)
            depth=-segs[:,2].min(); width=np.ptp(segs[:,0]) if geo=='plate' else 2*np.max(np.hypot(segs[:,0],segs[:,1]))
            rows.append(dict(geometry=geo,type=t,topology=NAME[t],rep=s+1,age_d=round(T,2),TRL=trl(p),MR=o['MR'],nLR=len(L),totLR=L[:,1].sum() if len(L) else 0,
                             meanLR=L[:,1].mean() if len(L) else 0,maxLR=L[:,1].max() if len(L) else 0,branched_zone_frac=L[:,0].max()/o['MR'] if len(L) else 0,
                             LR_ins_angle=np.median(L[:,3]),LR_angle_from_vertical=np.median(L[:,2]),depth=depth,width=width,width_depth=width/depth))
    pb_par=None
R=pd.DataFrame(rows); R.to_csv(f'{OUT}/replicate_metrics.csv',index=False)
pd.set_option('display.width',250)
print(R.groupby(['geometry','topology']).median(numeric_only=True).drop(columns='rep').round(2).to_string())
