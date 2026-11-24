"""Angle fit (as fit_angle.py / fit_angle_plate.py) on the kinetics-refit shapes, simulated to tracing day 5 (t0 + 5 d)."""
import cpb_paths as CP   # paths; changes into outputs/cplantbox/
import numpy as np, json, rsa_model as m, targets as tg
def shape(t):
    x=json.load(open(f'fit_v3_{t}.json'))['x']
    return dict(r=x['r'],mr_lmax=x['mr_lmax'],la=x['la'],lb=0.1,ln=x['ln'],lnk=x['lnk'],rL=x['rL'],sL=x['sL'],lmaxL=15.), x['t0']+5
def ang(P,T,plate,n):
    a=np.vstack([m.measure(m.simulate(P,T,seed=s+1,stochastic=True,plate=plate))['lats'] for s in range(n)]); return np.median(a[:,3]),np.median(a[:,2])
out={}
for t in ['C','B','dT']:
    S,T=shape(t); ti,tv=tg.T[t]['ins'],tg.T[t]['vert']
    g3=[(abs(ang(dict(S,theta=ti,tN=1.,tS=tS),T,False,12)[1]-tv),tS) for tS in [0.02,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6]]
    P3=dict(S,theta=ti,tN=1.,tS=min(g3)[1])
    g2=[]
    for th in range(72,92,2):
        for tS in [0.05,0.1,0.2,0.3,0.4,0.5]:
            for tN in [1,2]:
                i,v=ang(dict(S,theta=th,tN=tN,tS=tS),T,True,20); g2.append(((i-ti)**2+(v-tv)**2,th,tS,tN,i,v))
    b=min(g2); P2=dict(S,theta=b[1],tS=b[2],tN=b[3])
    out[t]={'3D':P3,'plate':P2,'T_day5':T}
    print(t,'target',round(ti,1),round(tv,1),'| 3D tS',P3['tS'],np.round(ang(P3,T,False,20),1),'| plate',b[1:4],round(b[4],1),round(b[5],1))
json.dump(out,open('params_v3.json','w'),indent=1)
