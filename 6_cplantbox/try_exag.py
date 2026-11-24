import cpb_paths as CP   # paths; changes into outputs/cplantbox/
import json, sys, numpy as np, plantbox as pb, rsa_model as m
F=json.load(open('params_fitted.json'))
def trl(p): return sum(o.getLength() for o in p.getOrgans(pb.root))
def run(P,plate=False,seed=1,stoch=False,target=30.):
    p=m.build(P,seed,stochastic=stoch,plate=plate); t=0
    while trl(p)<target-1e-3 and t<40:
        st=0.1 if target-trl(p)>1 else 0.01; p.simulate(st,False); t+=st
    o=m.measure(p); L=o['lats']; o['t']=t; return o
def show(name,P):
    o=run(P); L=o['lats']; L=L[np.argsort(L[:,0])]
    print(f"{name}: t {o['t']:.1f} d MR {o['MR']:.1f} nLR {len(L)} mean {L[:,1].mean():.2f} CV {L[:,1].std()/L[:,1].mean():.2f} branched {L[:,0].max()/o['MR']:.2f}")
    print('   pos ',np.round(L[:,0],1)); print('   len ',np.round(L[:,1],1))
