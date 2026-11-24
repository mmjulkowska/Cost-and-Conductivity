"""Growth animation with every root piece coloured by the v5 conductivity segment it is given, and Krs as it develops (v4 architectures, 2D plate view)."""
import cpb_paths as CP   # paths; changes into outputs/cplantbox/
import numpy as np, json, os, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.animation as an, plantbox as pb, rsa_model as m, conductance as K
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9})
P=json.load(open('params_v4.json')); ORDER=['C','B','dT']; END=12.; STEP=0.1; SEED=7
NAME={'C':'Christmas tree','B':'Broomstick','dT':'Droopy telephone pole'}; COL={'C':'#1b7a4a','B':'#b5651d','dT':'#2a5fa8'}
WIDTH={'MP':2.4,'MM':2.0,'MT':1.8,'MT0':1.8,'LP':1.2,'LM':1.1,'LT':1.0}
F={}
for t in ORDER:
    p=m.build(P[t]['plate'],SEED,stochastic=True,plate=True); t0=P[t]['T_day5']-5; tt=0.; fr=[]
    days=np.arange(0,END+1e-9,STEP)
    for d in days:
        while tt<t0+d-1e-9: st=min(STEP,t0+d-tt); p.simulate(st,False); tt+=st
        S=K.segments(p); k=K.krs(p,segs=S)[0] if S else 0.
        o=np.array([0,0,-0.5]); fr.append(dict(day=d,krs=k,trl=sum(s[2] for s in S),segs={c:[np.array([(s[4]-o)[[0,2]],(s[5]-o)[[0,2]]]) for s in S if s[3]==c] for c in K.CLASSES}))
    F[t]=fr
n=len(F['C'])
fig=plt.figure(figsize=(13,7)); axs=[fig.add_axes([0.06+0.235*i,0.16,0.22,0.72]) for i in range(3)]; axk=fig.add_axes([0.8,0.3,0.18,0.5])
ttl=fig.suptitle('',fontsize=11)
fig.legend(handles=[Line2D([],[],color=K.COLOR[c],lw=WIDTH[c]+1.5,label=K.LABEL[c]) for c in K.CLASSES],loc='lower center',ncol=4,frameon=False,fontsize=8.5,bbox_to_anchor=(0.45,0.0))
def upd(k):
    d=F['C'][k]['day']
    for ax,t in zip(axs,ORDER):
        ax.cla(); f=F[t][k]
        for c in K.CLASSES:
            if f['segs'][c]: ax.add_collection(LineCollection(f['segs'][c],colors=K.COLOR[c],linewidths=WIDTH[c],capstyle='round'))
        ax.set_xlim(-9,9); ax.set_ylim(-26,0.5); ax.set_aspect('equal'); ax.set_title(f"{NAME[t]}\n{f['trl']:.0f} cm · Krs {f['krs']*1e4:.1f}×10⁻⁴",fontsize=9,color=COL[t])
        for s in ('top','right'): ax.spines[s].set_visible(False)
        ax.set_xlabel('cm')
    axs[0].set_ylabel('depth (cm)')
    axk.cla(); axk.axvspan(0,5,color='#eeeeee')
    for t in ORDER:
        x=[f['day'] for f in F[t][:k+1]]; y=[f['krs']*1e4 for f in F[t][:k+1]]; axk.plot(x,y,color=COL[t],lw=2)
    axk.set_xlim(0,END); axk.set_ylim(0,max(f['krs'] for t in ORDER for f in F[t])*1e4*1.05); axk.set_xlabel('Tracing day'); axk.set_ylabel('Krs (×10⁻⁴ cm³ hPa⁻¹ d⁻¹)')
    for s in ('top','right'): axk.spines[s].set_visible(False)
    inside=d<=5; ttl.set_text(f"tracing day {d:4.1f}   ·   {'measured window (days 1–5)' if inside else 'EXTRAPOLATION beyond the measured days'}"); ttl.set_color('#333' if inside else '#c0392b')
    return []
a=an.FuncAnimation(fig,upd,frames=n,blit=False)
a.save('krs/conductance_segments_growth.mp4',writer=an.FFMpegWriter(fps=12,bitrate=3000),dpi=110)
a.save('krs/conductance_segments_growth.gif',writer=an.PillowWriter(fps=12),dpi=65)
