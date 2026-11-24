import cpb_paths as CP   # paths; changes into outputs/cplantbox/
import numpy as np, pandas as pd, json, os, glob, xml.etree.ElementTree as ET, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib.collections import LineCollection
OUT=os.environ.get('OUT','sims'); FIG=os.environ.get('FIG','figs'); os.makedirs(FIG,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42})
ORDER=['C','B','dT']; NAME={'C':'Christmas tree','B':'Broomstick','dT':'Droopy telephone pole'}
COL={'C':'#1b7a4a','B':'#b5651d','dT':'#2a5fa8'}; MRC='#222222'
R=pd.read_csv(f'{OUT}/replicate_metrics.csv'); TRL=R.TRL.median(); TAG=os.environ.get('TAG','')
def load(geo,t,s):
    pl=list(np.load(f'{OUT}/{geo}_{t}_rep{s:02d}_polylines.npy',allow_pickle=True)); o=pl[0][0].copy(); return [q-o for q in pl]   # seed at origin
def draw2d(ax,pls,c,x=0,y=1):
    ax.add_collection(LineCollection([q[:,[0,2]] for q in pls[1:]],colors=c,linewidths=0.9))
    ax.add_collection(LineCollection([pls[0][:,[0,2]]],colors=MRC,linewidths=1.6))
# ---------- measured examples (SmartRoot RSML, day 5, control)
PX=2.54/300; D=CP.RSML_DIR+'/'
sel=pd.read_csv(CP.RSML_SELECTED); P=pd.read_csv('rsml_plants.csv')
def measured(arbor,rsml):
    t=ET.parse(D+rsml).getroot()
    for mr in t.iter('root'):
        if mr.get('label') and arbor.endswith(mr.get('label')) and mr.find('geometry') is not None:
            pl=[np.array([[float(q.get('x')),float(q.get('y'))] for q in mr.find('geometry/polyline')])*PX]
            pl+=[np.array([[float(q.get('x')),float(q.get('y'))] for q in lr.find('geometry/polyline')])*PX for lr in mr.findall('root')]
            o=pl[0][0]; return [np.c_[q[:,0]-o[0],np.zeros(len(q)),-(q[:,1]-o[1])] for q in pl]
def pick(t):   # measured plant closest to the type's median MR, n LR and total LR length
    p=P[(P.type==t)&(P.cond=='control')].copy(); L=pd.read_csv('rsml_laterals.csv').groupby('arbor').len.sum(); p['tot']=p.arbor.map(L)
    z=sum(((p[c]-p[c].median())/p[c].std())**2 for c in ['MR','nLR','tot']); a=p.loc[z.idxmin(),'arbor']
    return a, sel.set_index('arbor_id').loc[a,'rsml']
NREP=int(R.rep.max()); ncol=6
# ---------- figure 1: 2D plate gallery
fig,axs=plt.subplots(3,ncol+1,figsize=(2.0*(ncol+1),3*3.1),sharex=True,sharey=True)
for i,t in enumerate(ORDER):
    a,f=pick(t); pls=measured(a,f); draw2d(axs[i,0],pls,'#888888')
    axs[i,0].set_title(f'measured (day 5)\n{a.split("_")[-3]}_{a.split("_")[-2]}',fontsize=7,color='#555')
    for j in range(ncol):
        draw2d(axs[i,j+1],load('plate',t,j+1),COL[t]); r=R[(R.geometry=='plate')&(R.type==t)&(R.rep==j+1)].iloc[0]
        axs[i,j+1].set_title(f'rep {j+1} · {r.age_d:.1f} d',fontsize=7)
    axs[i,0].set_ylabel(NAME[t],fontsize=11,color=COL[t],fontweight='bold')
for ax in axs.flat:
    ax.set_aspect('equal'); ax.set_xlim(-5,5); ax.set_ylim(-11,0.5); ax.tick_params(labelsize=7)
    for s in ['top','right']: ax.spines[s].set_visible(False)
for ax in axs[-1]: ax.set_xlabel('cm')
fig.suptitle(f'CPlantBox{TAG}, 2D plate: every simulated root system has total root length {TRL:.0f} cm (main root + laterals). Left column: a measured plant of that type.',fontsize=10)
fig.tight_layout(rect=(0,0,1,0.96))
for e in ('png','pdf'): fig.savefig(f'{FIG}/gallery_2D_plate.{e}',dpi=200)
plt.close(fig)
# ---------- figure 2: 3D gallery (same view for all; a 2 cm grid on the soil surface and a depth scale)
ncol3=6
fig=plt.figure(figsize=(2.1*ncol3,3*3.0))
for i,t in enumerate(ORDER):
    for j in range(ncol3):
        ax=fig.add_subplot(3,ncol3,i*ncol3+j+1,projection='3d'); pls=load('3D',t,j+1)
        ax.add_collection3d(Line3DCollection(pls[1:],colors=COL[t],linewidths=0.8)); ax.add_collection3d(Line3DCollection([pls[0]],colors=MRC,linewidths=1.5))
        g=np.arange(-3,3.1,2)
        for v in g: ax.plot([v,v],[-3,3],[0,0],color='#cccccc',lw=0.5); ax.plot([-3,3],[v,v],[0,0],color='#cccccc',lw=0.5)
        ax.plot([3,3],[3,3],[0,-10],color='#999999',lw=0.8)
        for z in (0,-5,-10): ax.plot([3,3],[2.6,3],[z,z],color='#999999',lw=0.8)
        if j==0: ax.text(3,3.2,-10.6,'10 cm',fontsize=6,color='#777')
        ax.set_xlim(-3,3); ax.set_ylim(-3,3); ax.set_zlim(-10.5,0.3); ax.set_box_aspect((1,1,1.75),zoom=1.15); ax.view_init(elev=18,azim=-55); ax.set_axis_off()
        r=R[(R.geometry=='3D')&(R.type==t)&(R.rep==j+1)].iloc[0]
        ax.set_title(f'rep {j+1} · {r.age_d:.1f} d',fontsize=7,pad=2)
    fig.text(0.008,0.95-(i+0.5)/3*0.94,NAME[t],rotation=90,va='center',fontsize=11,color=COL[t],fontweight='bold')
fig.suptitle(f'CPlantBox{TAG}, 3D soil: total root length {TRL:.0f} cm in every root system (grid on the soil surface: 2 cm squares; grey bar: 10 cm depth)',fontsize=10)
fig.subplots_adjust(left=0.03,right=0.99,top=0.95,bottom=0.01,wspace=0.0,hspace=0.02)
for e in ('png','pdf'): fig.savefig(f'{FIG}/gallery_3D.{e}',dpi=200)
plt.close(fig)
print('ok')
