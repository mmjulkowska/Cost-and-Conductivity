"""Figures for the soil runs (soil_runs/*.json)."""
import cpb_paths as CP   # paths; changes into outputs/cplantbox/
import json, glob, numpy as np, pandas as pd, os, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import soil_sim as S, conductance as K
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
OUT='soil'; os.makedirs(OUT,exist_ok=True)
ORDER=['C','B','dT']; NAME={'C':'Christmas tree','B':'Broomstick','dT':'Droopy telephone pole'}; COL={'C':'#1b7a4a','B':'#b5651d','dT':'#2a5fa8'}
SOILS=['sandy loam','loam','silty clay loam']; SCN={'wet':'Well watered','drydown':'Dry-down from container capacity','drytop':f'Dry top {S.DRYTOP_D:.0f} cm, wet below'}
rows=[]; prof=[]
for f in glob.glob('soil_runs/*.json'):
    d=json.load(open(f)); m=d['meta']
    for r in d['rows']:
        rows.append({k:v for k,v in r.items() if k not in ('h_profile','uptake_profile')})
    prof.append((m,d['rows']))
R=pd.DataFrame(rows); R.to_csv(f'{OUT}/soil_runs_timeseries.csv',index=False)
R['dayi']=np.floor(R.t+1e-9).astype(int)
# daily transpiration
D=[]
for key,g in R.groupby(['topology','soil','scenario','arch','rep']):
    g=g.sort_values('t'); c=np.interp(np.arange(0,11),np.r_[0,g.t],np.r_[0,g.cumT]); cp=np.interp(np.arange(0,11),np.r_[0,g.t],np.r_[0,g.cumTpot])
    for dd in range(10): D.append(dict(zip(['topology','soil','scenario','arch','rep'],key),day=dd+1,T=c[dd+1]-c[dd],Tpot=cp[dd+1]-cp[dd],cum=c[dd+1],cumpot=cp[dd+1]))
D=pd.DataFrame(D); D['rel']=D['T']/D.Tpot; D.to_csv(f'{OUT}/daily_transpiration.csv',index=False)
# ---- Figure 1: daily transpiration relative to potential
fig,axs=plt.subplots(3,3,figsize=(13,10),sharex=True,sharey=True)
for i,so in enumerate(SOILS):
    for j,sc in enumerate(['wet','drydown','drytop']):
        ax=axs[i,j]
        for t in ORDER:
            for arch,ls in (('frozen','-'),('growing','--')):
                g=D[(D.soil==so)&(D.scenario==sc)&(D.topology==t)&(D.arch==arch)].groupby('day').rel
                if len(g.size()): ax.plot(g.median().index,g.median()*100,ls,color=COL[t],lw=1.8,marker='o' if ls=='-' else None,ms=3)
        ax.set_ylim(0,105); ax.axhline(100,color='#bbb',lw=0.6)
        if i==0: ax.set_title(SCN[sc],fontsize=9.5)
        if j==0: ax.set_ylabel(f'{so}\ntranspiration (% of potential)')
        if i==2: ax.set_xlabel('Day in the pot')
fig.legend(handles=[Line2D([],[],color=COL[t],lw=2,label=NAME[t]) for t in ORDER]+[Line2D([],[],color='k',ls='-',label='root system fixed at 30 cm'),Line2D([],[],color='k',ls='--',label='root system keeps growing')],
           loc='lower center',ncol=5,frameon=False,bbox_to_anchor=(0.5,0.0))
fig.suptitle(f'Daily transpiration as % of the shoot\'s demand ({S.TPOT:g} cm³ d⁻¹ for every plant); pot {2*S.POT_R:g} cm wide × {S.POT_D:g} cm deep; median of replicates',fontsize=10.5)
fig.tight_layout(rect=(0,0.04,1,0.96))
for e in ('png','pdf'): fig.savefig(f'{OUT}/daily_transpiration.{e}',dpi=170)
plt.close(fig)
# ---- Figure 2: total water used over 10 days + day stress starts
C=D[D.day==10].groupby(['topology','soil','scenario','arch']).apply(lambda g:(g.cum/g.cumpot).median()*100,include_groups=False).rename('pct').reset_index()
def onset(g):
    x=g.sort_values('day'); b=x[x.rel<0.9]; return b.day.min() if len(b) else np.nan
O=D.groupby(['topology','soil','scenario','arch','rep']).apply(onset,include_groups=False).rename('onset').reset_index()
O=O.groupby(['topology','soil','scenario','arch']).onset.median().reset_index()
fig,axs=plt.subplots(2,2,figsize=(13,8))
for col,arch in enumerate(['frozen','growing']):
    for row,(df,val,lab) in enumerate(((C,'pct','Water transpired over 10 days (% of demand)'),(O,'onset','First day transpiration < 90% of demand'))):
        ax=axs[row,col]; x=0; ticks=[]
        for so in SOILS:
            for sc in ['drydown','drytop']:
                for k,t in enumerate(ORDER):
                    v=df[(df.soil==so)&(df.scenario==sc)&(df.topology==t)&(df.arch==arch)][val]
                    if len(v): ax.bar(x+k*0.27,v.iloc[0],width=0.25,color=COL[t])
                ticks.append((x+0.27,f'{so}\n{"dry-down" if sc=="drydown" else "dry top"}')); x+=1.1
            x+=0.4
        ax.set_xticks([a for a,_ in ticks]); ax.set_xticklabels([b for _,b in ticks],fontsize=7.5); ax.set_ylabel(lab)
        ax.set_title('root system fixed at 30 cm' if arch=='frozen' else 'root system keeps growing',fontsize=9.5)
        if row==0: ax.set_ylim(0,100)
fig.legend(handles=[Line2D([],[],color=COL[t],lw=6,label=NAME[t]) for t in ORDER],loc='upper center',ncol=3,frameon=False,bbox_to_anchor=(0.5,0.99))
fig.tight_layout(rect=(0,0,1,0.95))
for e in ('png','pdf'): fig.savefig(f'{OUT}/drought_summary.{e}',dpi=170)
plt.close(fig)
# ---- Figure 3: where water comes from, loam, rep 0, frozen: soil pressure head and uptake by depth on day 2 noon and day 8 noon
fig,axs=plt.subplots(2,4,figsize=(14,7.5),sharey=True)
z=-(np.arange(int(S.POT_D/S.DZ))+0.5)*S.DZ
for i,sc in enumerate(['drydown','drytop']):
    for t in ORDER:
        run=[r for (m_,r) in prof if m_['topology']==t and m_['soil']=='loam' and m_['scenario']==sc and m_['arch']=='frozen' and m_['rep']==0]
        if not run: continue
        rr=run[0]
        for k,day in enumerate((1.5,7.5)):
            r=min(rr,key=lambda q:abs(q['t']-day))
            axs[i,2*k].plot(-np.array(r['h_profile']),z,color=COL[t],lw=1.8); axs[i,2*k].set_xscale('log')
            axs[i,2*k+1].plot(np.array(r['uptake_profile']),z,color=COL[t],lw=1.8)
    for k,day in enumerate((2,8)):
        axs[i,2*k].set_title(f'{SCN[sc]}\nday {day}, noon: soil suction',fontsize=8.5); axs[i,2*k].set_xlabel('−h (hPa)')
        axs[i,2*k+1].set_title(f'day {day}, noon: root water uptake',fontsize=8.5); axs[i,2*k+1].set_xlabel('uptake (cm³ cm⁻³ d⁻¹)')
    axs[i,0].set_ylabel('depth (cm)')
fig.legend(handles=[Line2D([],[],color=COL[t],lw=2,label=NAME[t]) for t in ORDER],loc='lower center',ncol=3,frameon=False)
fig.suptitle('Loam, root system fixed at 30 cm (one replicate): where the soil dries and where the roots take up water',fontsize=10.5)
fig.tight_layout(rect=(0,0.04,1,0.95))
for e in ('png','pdf'): fig.savefig(f'{OUT}/profiles_loam.{e}',dpi=170)
plt.close(fig)
print(C.pivot_table(index=['soil','scenario','arch'],columns='topology',values='pct').round(1).to_string())
print(O.pivot_table(index=['soil','scenario','arch'],columns='topology',values='onset').to_string())
