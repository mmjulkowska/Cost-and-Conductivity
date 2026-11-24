import cpb_paths as CP   # paths; changes into outputs/cplantbox/
import numpy as np, json, plantbox as pb, rsa_model as m, targets as tg
from scipy.optimize import differential_evolution
BINS=np.arange(0,1.0001,0.1)
def sim_metrics(x, T=5.0):
    r,la,ln,lnk,rL,lmaxL=x
    P=dict(r=r,la=la,lb=0.1,ln=ln,lnk=lnk,rL=rL,lmaxL=lmaxL,theta=70,tN=1.)
    p=m.build(P); mr=p.getOrganRandomParameter(pb.root)[1]; mr.lmax=15.; p.initialize(False)
    for _ in range(int(T/0.1)): p.simulate(0.1,False)
    o=m.measure(p); L=o['lats']
    if len(L)==0: return None
    rel=L[:,0]/o['MR']
    cnt,_=np.histogram(np.clip(rel,0,.9999),BINS); lenp,_=np.histogram(np.clip(rel,0,.9999),BINS,weights=L[:,1])
    return dict(MR=o['MR'],nLR=len(L),totLR=L[:,1].sum(),cnt=cnt/cnt.sum(),prof=lenp/lenp.sum(),last=rel.max())
def make_obj(t):
    import pandas as pd
    Ls=pd.read_csv('rsml_laterals.csv'); Ls=Ls[(Ls.cond=='control')&(Ls.type==t)]
    cnt_t=np.mean([np.histogram(np.clip(g.rel,0,.9999),BINS)[0]/len(g) for _,g in Ls.groupby('arbor')],axis=0)
    tt=tg.T[t]; tt['cnt']=cnt_t
    def obj(x):
        r,la,ln,lnk,rL,lmaxL=x
        if ln+lnk*(0-(15-la-0.1)/2) < 0.05: return 10.   # spacing at the top must stay positive
        s=sim_metrics(x)
        if s is None: return 10.
        e=((s['MR']/tt['MR']-1)**2+(s['nLR']/tt['nLR']-1)**2+(s['totLR']/tt['totLR']-1)**2)*3
        e+=np.sum((s['cnt']-tt['cnt'])**2)+np.sum((s['prof']-tt['prof'])**2)
        return e
    return obj
if __name__=='__main__':
    B=[(0.8,3.),(0.3,6.),(0.05,1.5),(0.,0.15),(0.2,1.5),(0.5,6.)]
    res={}
    for t in ['B','C','dT']:
        o=make_obj(t); r=differential_evolution(o,B,seed=1,maxiter=150,popsize=15,tol=1e-8,polish=True)
        s=sim_metrics(r.x); res[t]=dict(x=list(r.x),err=r.fun)
        print(t,np.round(r.x,3),round(r.fun,4),{k:np.round(v,2) for k,v in s.items()})
        print('   target',{k:np.round(tg.T[t][k],2) for k in ['MR','nLR','totLR','last']}, np.round(tg.T[t]['cnt'],2), np.round(tg.T[t]['prof'],2))
    json.dump(res,open('fit_shape.json','w'),indent=1)
