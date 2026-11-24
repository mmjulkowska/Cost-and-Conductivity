import cpb_paths as CP   # paths; changes into outputs/cplantbox/
import numpy as np, pandas as pd, json, os, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import plantbox as pb, rsa_model as m, conductance as K
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
P=json.load(open('params_v4.json')); OUT='krs'; ORDER=['C','B','dT']
NAME={'C':'Christmas tree','B':'Broomstick','dT':'Droopy telephone pole'}; COL={'C':'#1b7a4a','B':'#b5651d','dT':'#2a5fa8'}
R=pd.read_csv(f'{OUT}/krs_timecourse.csv'); A=pd.read_csv(f'{OUT}/krs_at_30cm.csv')
WIDTH={'MP':2.4,'MM':2.0,'MT':1.8,'MT0':1.8,'LP':1.2,'LM':1.1,'LT':1.0}
def draw(ax,p,lim=None):
    S=K.segments(p); o=np.array([0,0,-0.5])
    for c in K.CLASSES:
        segs=[np.array([(s[4]-o)[[0,2]],(s[5]-o)[[0,2]]]) for s in S if s[3]==c]
        if segs: ax.add_collection(LineCollection(segs,colors=K.COLOR[c],linewidths=WIDTH[c],capstyle='round'))
    k,up,_=K.krs(p,segs=S); return k,up,sum(s[2] for s in S)
def grow_plate(t,seed,until_trl=None,until_day=None):
    pp=m.build(P[t]['plate'],seed,stochastic=True,plate=True); t0=P[t]['T_day5']-5; tt=0.
    while True:
        L=sum(o.getLength() for o in pp.getOrgans(pb.root))
        if until_trl and L>=until_trl: break
        if until_day is not None and tt>=t0+until_day-1e-9: break
        st=0.01 if (until_trl and L>until_trl-3) else 0.1; pp.simulate(st,False); tt+=st
    return pp,tt-t0
# ---------- Figure 1: segment maps (2D plate), at 30 cm and at tracing days 5 and 12
cols=[('30 cm total root length',dict(until_trl=30)),('tracing day 5 (last measured day)',dict(until_day=5)),('tracing day 12 (extrapolated)',dict(until_day=12))]
fig,axs=plt.subplots(3,3,figsize=(13,13.5),gridspec_kw=dict(width_ratios=[1,1,1.6]))
for i,t in enumerate(ORDER):
    for j,(lab,kw) in enumerate(cols):
        ax=axs[i,j]; pp,d=grow_plate(t,7,**kw); k,up,trl=draw(ax,pp)
        ax.set_aspect('equal'); ax.autoscale_view()
        xl=ax.get_xlim(); yl=ax.get_ylim()
        ax.set_title(f'{lab}\nday {d:.1f} · {trl:.0f} cm · Krs = {k*1e4:.1f}×10⁻⁴ cm³ hPa⁻¹ d⁻¹',fontsize=8)
        if j==0: ax.set_ylabel(f'{NAME[t]}\ndepth (cm)',fontsize=10,color=COL[t],fontweight='bold')
        ax.set_xlabel('cm')
lims=[(-4,4,-9,0.5),(-5,5,-11,0.5),(-9,9,-27,0.5)]
for j,(x0,x1,y0,y1) in enumerate(lims):
    for i in range(3): axs[i,j].set_xlim(x0,x1); axs[i,j].set_ylim(y0,y1)
h=[Line2D([],[],color=K.COLOR[c],lw=WIDTH[c]+1,label=K.LABEL[c]) for c in K.CLASSES]
fig.legend(handles=h,loc='lower center',ncol=4,frameon=False,fontsize=8.5,bbox_to_anchor=(0.5,0.0))
fig.suptitle('Which v5 conductivity each root piece gets (by distance from its own tip: true tip < 1.5 cm, then 3 cm steps), 2D plate view',fontsize=10.5)
fig.tight_layout(rect=(0,0.05,1,0.97))
for e in ('png','pdf'): fig.savefig(f'{OUT}/segment_map.{e}',dpi=170)
plt.close(fig)
# ---------- Figure 2: Krs
fig,axs=plt.subplots(2,3,figsize=(14,8.4))
base=R[(R.tip==1.5)&(R.scenario=='v5')]; endo=R[(R.tip==1.5)&(R.scenario=='lat_endo')]
for t in ORDER:
    for D,ls,lab in ((base,'-','laterals without endodermal barrier (B10)'),(endo,'--','laterals with endodermal barrier (B6)')):
        g=D[D.type==t].groupby('day'); x=g.Krs.median().index
        axs[0,0].plot(x,g.Krs.median()*1e4,ls,color=COL[t],lw=1.8,label=f'{NAME[t]}' if ls=='-' else None)
        if ls=='-': axs[0,0].fill_between(x,g.Krs.quantile(.1)*1e4,g.Krs.quantile(.9)*1e4,color=COL[t],alpha=0.12,lw=0)
        g2=D[D.type==t]; axs[0,1].plot(g2.groupby('day').TRL.median(),g.Krs.median()*1e4,ls,color=COL[t],lw=1.8)
        axs[0,2].plot(x,(g.Krs.median()/g.TRL.median())*1e6,ls,color=COL[t],lw=1.8)
for ax in axs[0]: ax.axvspan(-1,5,color='#eeeeee',zorder=0) if ax is not axs[0,1] else None
axs[0,0].set_xlim(0,12); axs[0,2].set_xlim(0,12); axs[0,0].set_xlabel('Tracing day (grey: measured window)'); axs[0,2].set_xlabel('Tracing day (grey: measured window)')
axs[0,0].set_ylabel('Krs (×10⁻⁴ cm³ hPa⁻¹ d⁻¹)'); axs[0,1].set_xlabel('Total root length (cm)'); axs[0,1].set_ylabel('Krs (×10⁻⁴ cm³ hPa⁻¹ d⁻¹)')
axs[0,1].axvline(30,color='#999',lw=0.8,ls=':'); axs[0,2].set_ylabel('Krs per root length (×10⁻⁶ cm² hPa⁻¹ d⁻¹)')
axs[0,0].legend(handles=axs[0,0].get_legend_handles_labels()[0]+[Line2D([],[],color='k',ls='-',label='laterals without endodermal barrier (B10)'),Line2D([],[],color='k',ls='--',label='laterals with endodermal barrier (B6)')],fontsize=7,frameon=False)
# at 30 cm
ax=axs[1,0]
for i,t in enumerate(ORDER):
    a=A[A.type==t]
    for off,col,mk in ((-0.17,'Krs','o'),(0.17,'Krs_lat_endo','^')):
        ax.scatter(i+off+np.random.default_rng(i).uniform(-.06,.06,len(a)),a[col]*1e4,s=16,color=COL[t],marker=mk,alpha=.8,lw=0)
        ax.plot([i+off-.12,i+off+.12],[a[col].median()*1e4]*2,color='k',lw=1.5)
ax.set_xticks(range(3)); ax.set_xticklabels([NAME[t].replace(' telephone','\ntelephone') for t in ORDER]); ax.set_ylabel('Krs at 30 cm total root length (×10⁻⁴)')
ax.set_title('Equal size (30 cm): ● laterals without, ▲ with endodermal barrier',fontsize=8.5)
# uptake share
ax=axs[1,1]; bottom=np.zeros(3)
for c in K.CLASSES:
    v=np.array([A[A.type==t][f'up_{c}'].median() for t in ORDER]); v=v
    ax.bar(range(3),v*100,bottom=bottom*100,color=K.COLOR[c],label=c,width=0.6); bottom+=v
ax.set_xticks(range(3)); ax.set_xticklabels([NAME[t].replace(' telephone','\ntelephone') for t in ORDER]); ax.set_ylabel('Share of water uptake at 30 cm (%)'); ax.legend(fontsize=7,frameon=False,ncol=2,loc='lower right')
ax.set_title('Where the water enters (laterals without endodermal barrier)',fontsize=8.5)
# tip zone sensitivity
ax=axs[1,2]
for t in ORDER:
    for tip,ls in ((0.5,':'),(1.5,'-'),(3.0,'--')):
        g=R[(R.type==t)&(R.tip==tip)&(R.scenario=='v5')].groupby('day').Krs.median(); ax.plot(g.index,g*1e4,ls,color=COL[t],lw=1.5)
ax.axvspan(-1,5,color='#eeeeee',zorder=0); ax.set_xlim(0,12); ax.set_xlabel('Tracing day'); ax.set_ylabel('Krs (×10⁻⁴ cm³ hPa⁻¹ d⁻¹)')
ax.legend(handles=[Line2D([],[],color='k',ls=ls,label=f'barrier-free tip zone {tip} cm') for tip,ls in ((0.5,':'),(1.5,'-'),(3.0,'--'))],fontsize=7,frameon=False)
ax.set_title('Sensitivity to how long the barrier-free tip is',fontsize=8.5)
fig.suptitle('Root-system hydraulic conductance (Krs) with the v5 segment conductivities, 12 replicates per topology (3D; median, band = 10–90%)',fontsize=10.5)
fig.tight_layout(rect=(0,0,1,0.96))
for e in ('png','pdf'): fig.savefig(f'{OUT}/krs_summary.{e}',dpi=170)
# numbers
s=A.groupby('type')[['day','TRL','Krs','Krs_lat_endo']].median(); s['ratio_endo']=s.Krs_lat_endo/s.Krs; print(s.round(6))
for d in (5.0,12.0):
    q=R[(R.day.round(2)==round(d-0.0,2))|(abs(R.day-d)<0.26)]
    print(d, q[(q.tip==1.5)].groupby(['type','scenario'])[['TRL','Krs']].median().round(6).to_string())
print(R[(R.scenario=='v5')&(abs(R.day-12)<0.26)].groupby(['type','tip']).Krs.median().unstack().round(6))
