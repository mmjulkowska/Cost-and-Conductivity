"""Animation of a pot drying out: three topologies side by side (root system fixed at 30 cm, replicate 0), soil water content as background,
roots coloured by their conductivity segment, transpiration (% of demand) below. Usage: python animate_soil.py "loam" drydown"""
import cpb_paths as CP   # paths; changes into outputs/cplantbox/
import json, sys, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, matplotlib.animation as an
from matplotlib.collections import LineCollection
import plantbox as pb, rsa_model as m, conductance as K, soil_sim as S
soil=sys.argv[1] if len(sys.argv)>1 else 'loam'; sc=sys.argv[2] if len(sys.argv)>2 else 'drydown'; arch='frozen'
ORDER=['C','B','dT']; NAME={'C':'Christmas tree','B':'Broomstick','dT':'Droopy telephone pole'}; COL={'C':'#1b7a4a','B':'#b5651d','dT':'#2a5fa8'}
WIDTH={'MP':2.2,'MM':1.9,'MT':1.7,'MT0':1.7,'LP':1.1,'LM':1.0,'LT':0.9}
P=json.load(open('params_v4.json')); so=S.Soil(soil)
runs={}; geo={}
for t in ORDER:
    runs[t]=json.load(open(f"soil_runs/{t}_{soil.replace(' ','-')}_{sc}_{arch}_r0.json"))['rows']
    p=m.build(dict(P[t]['3D'],pot=(S.POT_R,S.POT_D)),3000,stochastic=True)
    while sum(o.getLength() for o in p.getOrgans(pb.root))<30: p.simulate(0.01,False)
    Sg=K.segments(p); o=np.array([0,0,-0.5]); geo[t]={c:[np.array([(s[4]-o)[[0,2]],(s[5]-o)[[0,2]]]) for s in Sg if s[3]==c] for c in K.CLASSES}
n=min(len(r) for r in runs.values()); z=-(np.arange(int(S.POT_D/S.DZ))+0.5)*S.DZ
th_s,th_r=so.ts,so.tr
fig=plt.figure(figsize=(11,8.2)); axs=[fig.add_axes([0.06+0.25*i,0.33,0.2,0.58]) for i in range(3)]; axT=fig.add_axes([0.08,0.07,0.6,0.17]); cax=fig.add_axes([0.83,0.4,0.02,0.45])
ttl=fig.suptitle('',fontsize=11)
sm=plt.cm.ScalarMappable(cmap='YlGnBu',norm=plt.Normalize(th_r,th_s)); cb=fig.colorbar(sm,cax=cax); cb.set_label('soil water content (cm³ cm⁻³)')
def upd(k):
    for ax,t in zip(axs,ORDER):
        ax.cla(); r=runs[t][k]; th=so.theta(np.array(r['h_profile']))
        ax.imshow(np.tile(th[:,None],(1,2)),extent=(-S.POT_R,S.POT_R,-S.POT_D,0),aspect='auto',cmap='YlGnBu',vmin=th_r,vmax=th_s)
        for c in K.CLASSES:
            if geo[t][c]: ax.add_collection(LineCollection(geo[t][c],colors=K.COLOR[c],linewidths=WIDTH[c]))
        ax.set_xlim(-S.POT_R,S.POT_R); ax.set_ylim(-S.POT_D,0.3); ax.set_xlabel('cm')
        ax.set_title(f"{NAME[t]}\ntranspiring {100*r['Tact']/max(r['Tpot'],1e-9):.0f}% of demand" if r['Tpot']>0.05 else f"{NAME[t]}\nnight",fontsize=9,color=COL[t])
    axs[0].set_ylabel('depth (cm)')
    axT.cla()
    for t in ORDER:
        rr=runs[t][:k+1]; axT.plot([q['t'] for q in rr],[q['Tact'] for q in rr],color=COL[t],lw=1.6)
    rr=runs['C'][:k+1]; axT.plot([q['t'] for q in rr],[q['Tpot'] for q in rr],color='#999',lw=1,ls=':')
    axT.set_xlim(0,10); axT.set_ylim(0,S.TPOT*np.pi*1.05); axT.set_xlabel('Day in the pot'); axT.set_ylabel('transpiration\n(cm³ d⁻¹)')
    for s_ in ('top','right'): axT.spines[s_].set_visible(False)
    ttl.set_text(f"{soil}, {'dry-down from container capacity' if sc=='drydown' else 'dry top, wet below' if sc=='drytop' else 'well watered'} · day {runs['C'][k]['t']:.2f} · dotted = demand")
    return []
a=an.FuncAnimation(fig,upd,frames=n,blit=False)
stem=f"soil/anim_{soil.replace(' ','-')}_{sc}"
a.save(stem+'.mp4',writer=an.FFMpegWriter(fps=8,bitrate=2500),dpi=110); a.save(stem+'.gif',writer=an.PillowWriter(fps=8),dpi=65)
print(stem)
