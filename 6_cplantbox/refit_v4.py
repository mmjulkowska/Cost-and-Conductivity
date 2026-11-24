"""After switching to a constant-rate main root and type-specific lateral maxima, refit only la, rL and sL so that days 1-5 still match
(same objective as fit_v3: time courses, day-5 position profiles, lateral rate by age). Everything else is kept from v3."""
import cpb_paths as CP   # paths; changes into outputs/cplantbox/
import json, sys, numpy as np, rsa_model as m
from scipy.optimize import differential_evolution
import fit_v3 as F
P4=json.load(open('params_v4.json'))
def run(P,t0):
    p=m.build(P); out=[]; t=0.
    for day in range(1,6):
        while t<t0+day-1e-9: st=min(0.1,t0+day-t); p.simulate(st,False); t+=st
        o=m.measure(p); L=o['lats']; L=L[L[:,1]>0.1] if len(L) else L; out.append((o['MR'],len(L),L[:,1].sum() if len(L) else 0.,L))
    return out
def obj_for(t):
    T=F.targets_for(t); t0=json.load(open(f'fit_v3_{t}.json'))['x']['t0']; base=P4[t]['3D']
    def obj(x):
        P=dict(base,la=x[0],rL=x[1],sL=x[2])
        S=run(P,t0); e=0.
        for i,(mr,n,tot,L) in enumerate(S):
            e+=((mr-T['MR'][i])/T['MR'][-1])**2+((n-T['nLR'][i])/T['nLR'][-1])**2+((tot-T['tot'][i])/T['tot'][-1])**2
        e+=2*(((S[-1][1]/T['nLR'][-1])-1)**2+((S[-1][2]/T['tot'][-1])-1)**2)
        L=S[-1][3]; rel=L[:,0]/S[-1][0]; B=F.BINS
        c,_=np.histogram(np.clip(rel,0,.9999),B); pr,_=np.histogram(np.clip(rel,0,.9999),B,weights=L[:,1])
        e+=np.sum((c/c.sum()-T['cnt'])**2)+np.sum((pr/pr.sum()-T['prof'])**2)
        e+=0.5*np.mean(((x[1]+x[2]*np.arange(1,5)-T['rate'])/T['rate'])**2)
        return e
    return obj,t0
if __name__=='__main__':
    for t in ['C','B','dT']:
        o,t0=obj_for(t); r=differential_evolution(o,[(0.3,6),(0.05,1.5),(0,1.5)],seed=3,maxiter=60,popsize=12,tol=1e-7,polish=False)
        for g in ('plate','3D'): P4[t][g].update(la=float(r.x[0]),rL=float(r.x[1]),sL=float(r.x[2]))
        S=run(P4[t]['3D'],t0); T=F.targets_for(t)
        print(t,'la %.2f rL %.2f sL %.2f err %.3f'%(*r.x,r.fun))
        for i,k in enumerate(['MR','nLR','tot']): print('   ',k,'sim',np.round([s[i] for s in S],1),'target',np.round(T[k],1))
    json.dump(P4,open('params_v4.json','w'),indent=1)
