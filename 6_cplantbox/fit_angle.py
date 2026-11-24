"""Lateral angles. 3D: theta = measured insertion angle; tropism strength tS fitted so the median whole-lateral angle from vertical matches.
2D plate: the thin agar slab deflects laterals, so theta and tS are both fitted so that the simulated plate measurements
(insertion angle and angle from vertical) match the measured ones."""
import cpb_paths as CP   # paths; changes into outputs/cplantbox/
import numpy as np, json, rsa_model as m, targets as tg
F=json.load(open('fit_shape.json'))
def params(t, theta, tS, tN=1.):
    r,la,ln,lnk,rL,lmaxL=F[t]['x']
    return dict(r=r,la=la,lb=0.1,ln=ln,lnk=lnk,rL=rL,lmaxL=lmaxL,theta=theta,tN=tN,tS=tS,mr_lmax=15.)
def ang(P,plate,n=10):
    a=np.vstack([m.measure(m.simulate(P,5,seed=s+1,stochastic=True,plate=plate))['lats'] for s in range(n)])
    return np.median(a[:,3]), np.median(a[:,2])
out={}
for t in ['B','C','dT']:
    ti,tv=tg.T[t]['ins'],tg.T[t]['vert']
    best3=min(((abs(ang(params(t,ti,tS),False)[1]-tv),tS) for tS in [0.02,0.05,0.1,0.15,0.2,0.3,0.4,0.5]))
    g=[]
    for th in range(70,95,3):
        for tS in [0.02,0.05,0.1,0.2,0.3,0.4]:
            i,v=ang(params(t,th,tS),True); g.append(((i-ti)**2+(v-tv)**2,th,tS,i,v))
    b=min(g)
    P3=params(t,ti,best3[1]); P2=params(t,b[1],b[2])
    out[t]={'3D':P3,'plate':P2}
    print(t,'target ins %.1f vert %.1f'%(ti,tv),'| 3D tS',best3[1],'->',np.round(ang(P3,False),1),'| plate theta',b[1],'tS',b[2],'->',round(b[3],1),round(b[4],1))
json.dump(out,open('params_fitted.json','w'),indent=1)
