"""Lateral length vs insertion position: data-fitted vs exaggerated topologies (all at 30 cm total root length)."""
import cpb_paths as CP   # paths; changes into outputs/cplantbox/
import numpy as np, pandas as pd, os, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
ORDER=['C','B','dT']; NAME={'C':'Christmas tree','B':'Broomstick','dT':'Droopy telephone pole'}; COL={'C':'#1b7a4a','B':'#b5651d','dT':'#2a5fa8'}
def lats(d,g,t):
    out=[]
    for f in sorted(__import__('glob').glob(f'{d}/{g}_{t}_rep*_polylines.npy')):
        pl=np.load(f,allow_pickle=True); mr=pl[0]; cum=np.r_[0,np.cumsum(np.linalg.norm(np.diff(mr,axis=0),axis=1))]
        for q in pl[1:]: out.append((cum[np.argmin(np.linalg.norm(mr-q[0],axis=1))],np.linalg.norm(np.diff(q,axis=0),axis=1).sum()))
    return np.array(out)
fig,axs=plt.subplots(1,3,figsize=(12,3.8),sharey=True)
for ax,t in zip(axs,ORDER):
    for d,lab,a,mk in (('sims','fitted to the panel',0.35,'o'),('sims_exaggerated','exaggerated',0.9,'o')):
        L=lats(d,'3D',t); R=pd.read_csv(f'{d}/replicate_metrics.csv'); r=R[(R.geometry=='3D')&(R.type==t)]
        ax.scatter(L[:,0],L[:,1],s=10,color=COL[t],alpha=a,lw=0,label=f'{lab}: {r.nLR.median():.0f} laterals, mean {L[:,1].mean():.1f} cm, CV {L[:,1].std()/L[:,1].mean():.2f}')
    ax.set_title(NAME[t],color=COL[t],fontweight='bold'); ax.set_xlabel('Lateral insertion, distance from main-root base (cm)'); ax.legend(fontsize=7,frameon=False,loc='upper right')
axs[0].set_ylabel('Lateral root length (cm)')
fig.suptitle('Lateral length along the main root, 12 replicates per topology (3D), all at 30 cm total root length. Pale: fitted to the panel; dark: exaggerated.',fontsize=9.5)
fig.tight_layout(rect=(0,0,1,0.93))
os.makedirs('figs_exaggerated',exist_ok=True)
for e in ('png','pdf'): fig.savefig(f'figs_exaggerated/lateral_profiles_fitted_vs_exaggerated.{e}',dpi=200)
