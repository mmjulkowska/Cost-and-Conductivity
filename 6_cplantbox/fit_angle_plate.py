"""Refine the 2D plate angle fit with 30 seeds and tropism number tN in {1,2,3}."""
import cpb_paths as CP   # paths; changes into outputs/cplantbox/
import numpy as np, json, rsa_model as m, targets as tg
F=json.load(open('params_fitted.json')); S=json.load(open('fit_shape.json'))
def P2(t,th,tS,tN):
    r,la,ln,lnk,rL,lmaxL=S[t]['x']; return dict(r=r,la=la,lb=0.1,ln=ln,lnk=lnk,rL=rL,lmaxL=lmaxL,theta=th,tN=tN,tS=tS,mr_lmax=15.)
def ang(P,n=30):
    a=np.vstack([m.measure(m.simulate(P,5,seed=s+1,stochastic=True,plate=True))['lats'] for s in range(n)]); return np.median(a[:,3]),np.median(a[:,2])
for t in ['C','B','dT']:
    ti,tv=tg.T[t]['ins'],tg.T[t]['vert']; g=[]
    for th in range(72,92,2):
        for tS in [0.05,0.1,0.2,0.3,0.4,0.5]:
            for tN in [1,2]:
                i,v=ang(P2(t,th,tS,tN)); g.append(((i-ti)**2+(v-tv)**2,th,tS,tN,i,v))
    b=min(g); F[t]['plate']=P2(t,b[1],b[2],b[3]); print(t,round(ti,1),round(tv,1),'->',b[1:4],round(b[4],1),round(b[5],1))
json.dump(F,open('params_fitted.json','w'),indent=1)
