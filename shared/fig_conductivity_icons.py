"""Variant with photo/GRANAR half-and-half cross-section icons (icons from make_icons.py). Root schematic with barriers + kr and Kx per fragment, Tukey letters. Body text kept minimal; details go in the legend."""
import pandas as pd, numpy as np, matplotlib, textwrap, sys; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle
from PIL import Image
import os
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,
                     'axes.edgecolor':'#8a8984','xtick.color':'#52514e','ytick.color':'#52514e','pdf.fonttype':42})
MAIN,LAT,EXO,ENDO,OFF,TXT2,TISSUE='#2a78d6','#1baf7a','#eda100','#e34948','#c9c8c2','#52514e','#f3efe6'
BAR={'MP':(1,1),'MM':(1,1),'MT':(0,1),'MT0':(0,0),'LP':(1,0),'LM':(1,0),'LT':(0,0)}

def draw(data_csv, letters_kr, letters_kx, title, out, legend_text=None, icon_dir='icons', diam_um=None):
    # icon radius proportional to the photographed root diameter (largest = 1.2 plot units)
    dmax=max(diam_um.values()) if diam_um else 1; UPU=1.2/(dmax/2)
    RAD=(lambda f: diam_um[f]/2*UPU) if diam_um else (lambda f: 1.0)
    d=pd.read_csv(data_csv)
    H=8.3; WR=185; LH=0 if legend_text is None else sum(0.158*len(textwrap.wrap(pp,WR))+0.14 for pp in legend_text.split('\n\n'))+0.35
    fig=plt.figure(figsize=(14,H+LH))
    top=1-0.85/(H+LH); bot=LH/(H+LH)+0.45/(H+LH)
    gs=fig.add_gridspec(2,2,width_ratios=[1.0,1.2],hspace=0.38,wspace=0.10,left=0.01,right=0.985,top=top,bottom=bot)
    ax=fig.add_subplot(gs[:,0]); ax.set_xlim(-0.4,11.6); ax.set_ylim(-0.6,14.6); ax.set_aspect('equal'); ax.axis('off')
    def strip(p0,p1,w0,w1,exo,endo):
        p0,p1=np.asarray(p0,float),np.asarray(p1,float); v=p1-p0; n=np.array([-v[1],v[0]])/np.linalg.norm(v)
        ax.add_patch(Polygon([p0+n*w0,p1+n*w1,p1-n*w1,p0-n*w0],fc=TISSUE,ec='none',zorder=1))
        for s in (1,-1):
            ax.plot(*np.array([p0+s*n*w0*0.92,p1+s*n*w1*0.92]).T,color=EXO if exo else OFF,lw=4 if exo else 1.1,ls='-' if exo else (0,(3,2)),solid_capstyle='butt',zorder=3)
            ax.plot(*np.array([p0+s*n*w0*0.42,p1+s*n*w1*0.42]).T,color=ENDO if endo else OFF,lw=3 if endo else 1.1,ls='-' if endo else (0,(3,2)),solid_capstyle='butt',zorder=3)
    def apex(p,u,w):   # rounded root cap: half-ellipse, slightly elongated along the root axis
        u=np.asarray(u,float)/np.linalg.norm(u); n=np.array([-u[1],u[0]]); p=np.asarray(p,float)
        t=np.linspace(-np.pi/2,np.pi/2,60)
        pts=[p+n*w*np.sin(a)+u*w*1.35*np.cos(a) for a in t]
        ax.add_patch(Polygon([p+n*w*1.0005]+pts[::-1]+[p-n*w*1.0005],closed=True,fc=TISSUE,ec='none',zorder=1))
    def icon(c,r,frag):
        im=np.asarray(Image.open(os.path.join(icon_dir,f'icon_{frag}.png')))
        ax.imshow(im,extent=(c[0]-r,c[0]+r,c[1]-r,c[1]+r),zorder=7,interpolation='lanczos')
    width=lambda y:0.30+0.025*y
    mx=3.2
    for y0,y1,b in [(0.6,1.5,(0,0)),(1.5,4.5,(0,1)),(4.5,10.6,(1,1))]: strip((mx,y1),(mx,y0),width(y1),width(y0),*b)
    apex(np.array([mx,0.6]),(0,-1),width(0.6))
    ax.text(mx,11.2,'Main root',ha='center',fontsize=10.5,fontweight='bold',color=MAIN)
    u=np.array([1,-0.30]); u/=np.linalg.norm(u); nrm=np.array([-u[1],u[0]])
    Lbase=np.array([mx+width(10.4)*0.95,10.4]); Llen=6.9; P=lambda s: Lbase+u*(Llen-s)
    for s0,s1,b in [(0.0,1.5,(0,0)),(1.5,Llen,(1,0))]: strip(P(s1),P(s0),0.22+0.012*s1,0.22+0.012*s0,*b)
    apex(P(0),u,0.22)
    lt=P(1.9)-nrm*0.62; ax.text(*lt,'Lateral root',ha='center',va='center',fontsize=10.5,fontweight='bold',color=LAT,rotation=np.degrees(np.arctan2(u[1],u[0])))
    for f,y in [('MP',9.0),('MM',6.0),('MT',3.0),('MT0',0.6)]:
        ax.plot([mx-width(y)-0.12,mx+width(y)+0.12],[y,y],color='#0b0b0b',lw=1.2,zorder=5)
        r=RAD(f); c=np.array([2.05-r,y]); ax.plot([c[0]+r,mx-width(y)-0.12],[y,y],color='#b9b6ad',lw=0.7,zorder=2); icon(c,r,f)
        ax.text(mx+width(y)+0.25,y,f,fontsize=9.5,fontweight='bold',va='center',zorder=8)
    for f,s in [('LP',6.0),('LM',3.0),('LT',0.0)]:
        p=P(s); ax.plot(*np.array([p+nrm*0.36,p-nrm*0.36]).T,color='#0b0b0b',lw=1.2,zorder=5)
        r=RAD(f); c=P(s)+np.array([0,1.15+r]); ax.plot([p[0],c[0]],[p[1],c[1]-r],color='#b9b6ad',lw=0.7,zorder=2); icon(c,r,f)
        ax.text(c[0],c[1]+r+0.08,f,ha='center',va='bottom',fontsize=9.5,fontweight='bold')
    lx,ly=6.6,5.6
    for k,(c,lw,ls,t) in enumerate([(EXO,4,'-','Exodermis: suberin + lignin'),(ENDO,3,'-','Endodermis: lignin'),(OFF,1.2,(0,(3,2)),'No barrier')]):
        yy=ly-0.5*k; ax.plot([lx,lx+0.6],[yy,yy],color=c,lw=lw,ls=ls); ax.text(lx+0.75,yy,t,va='center',fontsize=8.5)
    if diam_um:
        L=200*UPU; ax.plot([lx,lx+L],[ly-1.35,ly-1.35],color='#333',lw=2.5,solid_capstyle='butt'); ax.text(lx+L+0.15,ly-1.35,'200 µm',va='center',fontsize=7.8,color=TXT2)
    ax.text(lx,ly-1.75,'Circle size = measured root diameter.\nCross-sections: left, Emily Miller\'s section;\nright, GRANAR anatomy (blue xylem,\ngreen phloem, lilac stele, pink endodermis,\nyellow exodermis; barrier walls coloured)',va='top',fontsize=7.8,color=TXT2,linespacing=1.3)
    ax.text(-0.3,13.9,'A',fontsize=12,fontweight='bold')
    cats=['Proximal','Middle','Tip zone','True tip']; X={'MP':0,'MM':1,'MT':2,'MT0':3,'LP':0,'LM':1,'LT':3}
    rng=np.random.default_rng(3)
    for i,(col,yl,ttl,L) in enumerate([('kr_cm_hPa_d','kr (cm hPa$^{-1}$ d$^{-1}$)','B  Radial conductivity',letters_kr),
                                       ('Kx_cm4_hPa_d','Kx (cm$^4$ hPa$^{-1}$ d$^{-1}$)','C  Axial conductance',letters_kx)]):
        a=fig.add_subplot(gs[i,1])
        for frags,c,dx,name in ((['MP','MM','MT','MT0'],MAIN,-0.16,'Main root'),(['LP','LM','LT'],LAT,0.16,'Lateral root')):
            for f in frags:
                s=d[d.fragment==f]; x=X[f]+dx
                a.scatter(x+rng.uniform(-0.06,0.06,len(s)),s[col],s=20,c=c,alpha=0.45,lw=0,zorder=2)
                v=s[col].median(); a.plot([x-0.12,x+0.12],[v,v],color=c,lw=3,solid_capstyle='round',zorder=3)
                if f in L: a.text(x,s[col].max()*1.15,L[f],ha='center',va='bottom',fontsize=9,color='#0b0b0b')
                elif f=='MT0' and 'MT' in L: a.text(x,s[col].max()*1.15,L['MT'],ha='center',va='bottom',fontsize=9,color='#9a9993',style='italic')
            a.scatter([],[],s=30,c=c,label=name)
        a.set_yscale('log'); a.set_xticks(range(4)); a.set_xticklabels(cats); a.set_xlim(-0.55,3.55)
        a.set_ylabel(yl); a.set_title(ttl,loc='left',fontsize=10.5,fontweight='bold')
        a.grid(axis='y',which='major',color='#e6e5e0',lw=0.6); a.set_axisbelow(True)
        lo,hi=a.get_ylim(); a.set_ylim(lo,hi*1.9)
        if i==0: a.legend(frameon=False,loc='upper left',fontsize=8.5,ncol=2)
    fig.suptitle(title,x=0.012,y=1-0.15/(H+LH),ha='left',va='top',fontsize=12.5,fontweight='bold')
    if legend_text:
        paras=legend_text.split('\n\n'); y=(LH-0.1)/(H+LH)
        for k,ptxt in enumerate(paras):
            lines=textwrap.wrap(ptxt,WR)
            fig.text(0.012,y,'\n'.join(lines),fontsize=8.2,va='top',ha='left',color='#1f1f1d',fontweight='bold' if k==0 else 'normal',linespacing=1.35)
            y-=(0.158*len(lines)+0.14)/(H+LH)
    fig.savefig(out+'.png',dpi=220); fig.savefig(out+'.pdf'); plt.close(fig)
