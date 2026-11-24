"""Growth animations from the v3 model: each topology from germination to tracing day END (default 12).
Days 1-5 are inside the measured window; later days are extrapolation (shaded red label).
Writes anim/<type>_growth.mp4 / .gif (2D plate + 3D) and anim/all_topologies_growth.mp4 / .gif (3D, side by side)."""
import cpb_paths as CP   # paths; changes into outputs/cplantbox/
import numpy as np, json, os, sys, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.animation as an, plantbox as pb, rsa_model as m
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9})
END=float(os.environ.get('END',12)); STEP=0.1; FPS=12; SEED=int(os.environ.get('SEED',3))
PAR=json.load(open(os.environ.get('PARAMS','params_v3.json'))); AD=os.environ.get('ANIM','anim'); ORDER=['C','B','dT']
NAME={'C':'Christmas tree','B':'Broomstick','dT':'Droopy telephone pole'}; COL={'C':'#1b7a4a','B':'#b5651d','dT':'#2a5fa8'}
os.makedirs(AD,exist_ok=True)
def series(P,t0,plate):
    p=m.build(P,SEED,stochastic=True,plate=plate); t=0.; frames=[]
    for k in range(int(round((t0+END)/STEP))+1):
        if k: p.simulate(STEP,False); t+=STEP
        pl=[np.array([[n.x,n.y,n.z] for n in r.getNodes()]) for r in p.getOrgans(pb.root)]
        o=np.array([0,0,-0.5]); pl=[q-o for q in pl if len(q)] or [np.zeros((1,3))]
        L=[np.linalg.norm(np.diff(q,axis=0),axis=1).sum() for q in pl]
        frames.append(dict(t=t,day=t-t0,pl=pl,TRL=sum(L),MR=L[0] if L else 0,nLR=sum(1 for x in L[1:] if x>0.1)))
    return frames
S={t:{g:series(PAR[t][g],PAR[t]['T_day5']-5,g=='plate') for g in ('plate','3D')} for t in ORDER}
def ext(fr):
    a=np.vstack([q for q in fr[-1]['pl'] if len(q)]); return a
def label(ax_or_fig,f,x,y,**kw):
    d=f['day']; inside=d<=5.0
    s=f"tracing day {max(d,0):4.1f}   (root age {f['t']:.1f} d)\n{'measured window (days 1–5)' if inside else 'EXTRAPOLATION beyond the measured days'}"
    return s, ('#333333' if inside else '#c0392b')
def draw2d(ax,f,c,lim):
    ax.cla(); pl=f['pl']
    if len(pl)>1: ax.add_collection(LineCollection([q[:,[0,2]] for q in pl[1:] if len(q)>1],colors=c,linewidths=0.9))
    if len(pl[0])>1: ax.add_collection(LineCollection([pl[0][:,[0,2]]],colors='#222',linewidths=1.6))
    ax.set_xlim(-lim[0],lim[0]); ax.set_ylim(-lim[1],0.5); ax.set_aspect('equal'); ax.set_xlabel('cm'); ax.set_ylabel('depth (cm)')
    for s in ('top','right'): ax.spines[s].set_visible(False)
def draw3d(ax,f,c,lim,az,zoom=1.45):
    ax.cla(); pl=f['pl']
    if len(pl)>1: ax.add_collection3d(Line3DCollection([q for q in pl[1:] if len(q)>1],colors=c,linewidths=0.8))
    if len(pl[0])>1: ax.add_collection3d(Line3DCollection([pl[0]],colors='#222',linewidths=1.5))
    R=lim[0]; g=np.arange(-R,R+0.01,2)
    for v in g: ax.plot([v,v],[-R,R],[0,0],color='#ddd',lw=0.5); ax.plot([-R,R],[v,v],[0,0],color='#ddd',lw=0.5)
    ax.plot([R,R],[R,R],[0,-10],color='#999',lw=1); ax.text(R,R,-10.8,'10 cm',fontsize=7,color='#777')
    ax.set_xlim(-R,R); ax.set_ylim(-R,R); ax.set_zlim(-lim[1],0.3); ax.set_box_aspect((1,1,lim[1]/R),zoom=zoom); ax.view_init(elev=14,azim=az); ax.set_axis_off()
def save(anim,stem):
    anim.save(f'{AD}/{stem}.mp4',writer=an.FFMpegWriter(fps=FPS,bitrate=2500),dpi=110)
    anim.save(f'{AD}/{stem}.gif',writer=an.PillowWriter(fps=FPS),dpi=70)
TRL30={t:next((f['day'] for f in S[t]['3D'] if f['TRL']>=30),None) for t in ORDER}
# ---- one animation per topology: plate (left) + 3D (right)
for t in ORDER:
    F2,F3=S[t]['plate'],S[t]['3D']
    a2=ext(F2); a3=ext(F3); lim2=(max(4,np.abs(a2[:,0]).max()*1.05),max(10,-a2[:,2].min()*1.05)); lim3=(max(4,np.ceil(np.abs(a3[:,:2]).max())),max(10,-a3[:,2].min()*1.05))
    fig=plt.figure(figsize=(11,6.4)); ax1=fig.add_axes([0.06,0.1,0.42,0.72]); ax2=fig.add_axes([0.5,0.06,0.5,0.8],projection='3d')
    ttl=fig.suptitle('',fontsize=11); info=fig.text(0.5,0.02,'',ha='center',fontsize=9)
    def upd(k):
        f2,f3=F2[k],F3[k]; draw2d(ax1,f2,COL[t],lim2); draw3d(ax2,f3,COL[t],lim3,-60+0.6*k,zoom=1.15 if lim3[1]/lim3[0]<2 else 0.85)
        ax1.set_title('2D plate',fontsize=9); ax2.set_title('3D soil',fontsize=9)
        s,c=label(fig,f3,0,0); ttl.set_text(f'{NAME[t]}\n{s}'); ttl.set_color(c)
        info.set_text(f"3D: total root length {f3['TRL']:.0f} cm · main root {f3['MR']:.1f} cm · {f3['nLR']} laterals"
                      + (f"   ·   reaches 30 cm on tracing day {TRL30[t]:.1f}" if TRL30[t] is not None and f3['day']>=TRL30[t] else ''))
        return []
    save(an.FuncAnimation(fig,upd,frames=len(F3),blit=False),f'{t}_growth'); plt.close(fig); print(t,'done')
# ---- all three together (3D), same scale and clock
A=np.vstack([ext(S[t]['3D']) for t in ORDER]); lim=(max(4,np.ceil(np.abs(A[:,:2]).max())),max(10,-A[:,2].min()*1.05))
fig=plt.figure(figsize=(12,6.4)); axs=[fig.add_subplot(1,3,i+1,projection='3d') for i in range(3)]; ttl=fig.suptitle('',fontsize=11)
n=min(len(S[t]['3D']) for t in ORDER)
def upd(k):
    # same clock = same tracing day for all three
    d=S['C']['3D'][k]['day']
    for ax,t in zip(axs,ORDER):
        F=S[t]['3D']; j=int(np.argmin([abs(f['day']-d) for f in F])); f=F[j]
        draw3d(ax,f,COL[t],lim,-60+0.6*k,zoom=1.0); ax.set_title(f"{NAME[t]}\n{f['TRL']:.0f} cm total · {f['nLR']} laterals",fontsize=9,color=COL[t])
    s,c=label(fig,S['C']['3D'][k],0,0); ttl.set_text(s); ttl.set_color(c); return []
fig.subplots_adjust(left=0.01,right=0.99,top=0.84,bottom=0.02,wspace=0.0)
save(an.FuncAnimation(fig,upd,frames=n,blit=False),'all_topologies_growth'); plt.close(fig)
json.dump({t:{'day_TRL_30cm':TRL30[t],'END_tracing_day':END,'final':{k:S[t]['3D'][-1][k] for k in ('TRL','MR','nLR')}} for t in ORDER},open(f'{AD}/summary.json','w'),indent=1)
print(json.load(open(f'{AD}/summary.json')))
