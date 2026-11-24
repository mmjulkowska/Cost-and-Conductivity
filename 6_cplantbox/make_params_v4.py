"""v4 = v3 with (1) the main root growing at a constant rate (least-squares fit to the day 1-5 main-root lengths, through zero at germination)
and (2) a type-specific maximal lateral length: Broomstick 8 cm > Christmas tree 6 cm > Droopy telephone pole 4 cm.
Lateral spacing is unchanged within the old main-root lmax and continues with the same widening beyond it (rsa_model: extend_mr)."""
import cpb_paths as CP   # paths; changes into outputs/cplantbox/
import json, numpy as np, rsa_model as m
K=json.load(open('kinetics_targets_v3.json')); P3=json.load(open('params_v3.json'))
LMAX={'B':8.,'C':6.,'dT':4.}; MRMAX=200.
out={}
for t in ['C','B','dT']:
    x=json.load(open(f'fit_v3_{t}.json'))['x']; t0=x['t0']; d=np.arange(1,6); y=np.array(K[t]['MR'])
    v=float(np.sum(y*(t0+d))/np.sum((t0+d)**2))
    out[t]={'T_day5':P3[t]['T_day5'],'MR_rate':v}
    for g in ('plate','3D'):
        P=dict(P3[t][g]); P.update(r=v,gfMR=2,extend_mr=True,lmaxL=LMAX[t]); out[t][g]=P
json.dump(out,open('params_v4.json','w'),indent=1)
if __name__=='__main__':
    for t in ['C','B','dT']:
        x=json.load(open(f'fit_v3_{t}.json'))['x']; t0=x['t0']; p=m.build(out[t]['3D']); tt=0.; S=[]
        for dd in range(1,6):
            while tt<t0+dd-1e-9: st=min(0.1,t0+dd-tt); p.simulate(st,False); tt+=st
            o=m.measure(p); L=o['lats']; L=L[L[:,1]>0.1]; S.append((o['MR'],len(L),L[:,1].sum()))
        S=np.array(S); print(t,'MR rate %.2f cm/d, LR max %.0f cm'%(out[t]['MR_rate'],LMAX[t]))
        for i,k in enumerate(['MR','nLR','tot']): print('   ',k,'sim',np.round(S[:,i],1),'target',np.round(K[t][k],1))
