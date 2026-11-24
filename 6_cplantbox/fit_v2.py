"""Refit with growth kinetics: PimpiBig2 day 1-5 tracings (same plants tracked daily) give the time course of main-root length,
lateral number and total lateral length, and lateral elongation rate by lateral age. Sizes are anchored to the panel day-5 medians
(target(d) = Big2(d) x panel(5) / Big2(5)); the day-5 position profiles come from the RSML as before.
Free parameters: t0 (root age at tracing day 0), MR r and lmax, la, ln, lnk, lateral rL and acceleration sL."""
import cpb_paths as CP   # paths; changes into outputs/cplantbox/
import numpy as np, json, sys, pandas as pd, plantbox as pb, rsa_model as m, targets as tg
from scipy.optimize import differential_evolution
K=json.load(open('kinetics_targets.json')); BINS=np.arange(0,1.0001,0.1)
NAMES=['t0','r','mr_lmax','la','ln','lnk','rL','sL']
BOUNDS=[(0,4),(0.5,4),(8,40),(0.3,6),(0.05,1.5),(0,0.15),(0.05,1.5),(0,1.5)]
def P_of(x,t='C'):
    d=dict(zip(NAMES,x)); return dict(r=d['r'],mr_lmax=d['mr_lmax'],la=d['la'],lb=0.1,ln=d['ln'],lnk=d['lnk'],rL=d['rL'],sL=d['sL'],lmaxL=15.,theta=70,tN=1.,tS=0.05), d['t0']
def run(x):
    P,t0=P_of(x); p=m.build(P); out=[]; t=0.
    for day in range(1,6):
        while t < t0+day-1e-9: st=min(0.1,t0+day-t); p.simulate(st,False); t+=st
        o=m.measure(p); L=o['lats']; L=L[L[:,1]>0.1] if len(L) else L
        out.append((o['MR'],len(L),L[:,1].sum() if len(L) else 0.,L))
    return out
def targets_for(t):
    k=K[t]; T=tg.T[t]
    sc=lambda v,a: np.array(v)*a/v[-1]
    Ls=pd.read_csv('rsml_laterals.csv'); Ls=Ls[(Ls.cond=='control')&(Ls.type==t)]
    cnt=np.mean([np.histogram(np.clip(g.rel,0,.9999),BINS)[0]/len(g) for _,g in Ls.groupby('arbor')],axis=0)
    return dict(MR=sc(k['MR'],T['MR']),nLR=sc(k['nLR'],T['nLR']),tot=sc(k['totLR'],T['totLR']),cnt=cnt,prof=T['prof'],
                rate=np.array([k['rate_by_age'][str(a)] for a in (1,2,3,4)]))
class Obj:
    def __init__(self,t): self.T=targets_for(t)
    def __call__(self,x):
        P,_=P_of(x)
        if P['ln']+P['lnk']*(0-(P['mr_lmax']-P['la']-0.1)/2)<0.05: return 50.
        try: S=run(x)
        except Exception: return 50.
        T=self.T; e=0.
        for i,(mr,n,tot,L) in enumerate(S):
            e+=((mr-T['MR'][i])/T['MR'][-1])**2+((n-T['nLR'][i])/T['nLR'][-1])**2+((tot-T['tot'][i])/T['tot'][-1])**2
        e+=2*(((S[-1][0]/T['MR'][-1])-1)**2+((S[-1][1]/T['nLR'][-1])-1)**2+((S[-1][2]/T['tot'][-1])-1)**2)
        L=S[-1][3]
        if len(L)==0: return 50.
        rel=L[:,0]/S[-1][0]; c,_=np.histogram(np.clip(rel,0,.9999),BINS); pr,_=np.histogram(np.clip(rel,0,.9999),BINS,weights=L[:,1])
        e+=np.sum((c/c.sum()-T['cnt'])**2)+np.sum((pr/pr.sum()-T['prof'])**2)
        rate=P['rL']+P['sL']*np.arange(1,5); e+=0.5*np.mean(((rate-T['rate'])/T['rate'])**2)
        return e
if __name__=='__main__':
    t=sys.argv[1]
    r=differential_evolution(Obj(t),BOUNDS,seed=2,maxiter=int(sys.argv[2]) if len(sys.argv)>2 else 120,popsize=15,tol=1e-7,workers=2,updating='deferred',polish=False)
    S=run(r.x); T=targets_for(t)
    res=dict(x=dict(zip(NAMES,map(float,r.x))),err=float(r.fun),sim=dict(MR=[s[0] for s in S],nLR=[s[1] for s in S],tot=[s[2] for s in S]),
             target=dict(MR=T['MR'].tolist(),nLR=T['nLR'].tolist(),tot=T['tot'].tolist(),rate=T['rate'].tolist()))
    json.dump(res,open(f'fit_v2_{t}.json','w'),indent=1)
    print(t,{k:round(v,3) for k,v in res['x'].items()},'err',round(r.fun,4))
    for k in ('MR','nLR','tot'): print('  ',k,'sim',np.round(res['sim'][k],1),'target',np.round(res['target'][k],1))
