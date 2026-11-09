"""Cost vs conductivity figure (run after rep_stats.py). Outputs in <outputs>/costs/."""
import os, sys; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from paths import *
os.chdir(COSTS)  # cost tables live in <outputs>/costs/
import pandas as pd, numpy as np, textwrap, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.ticker as mt
from matplotlib.lines import Line2D
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':'#8a8984','xtick.color':'#52514e','ytick.color':'#52514e','pdf.fonttype':42})
MAIN,LAT='#2a78d6','#1baf7a'; COL={'M':MAIN,'L':LAT}; MK={'P':'o','M':'s','T':'^'}
P=pd.read_csv('replicate_level_costs_conductivity.csv'); R=pd.read_csv('correlations_replicate_level.csv')
COST=[('C_ng_mm','Building cost (ng C mm$^{-1}$)','Building C'),('maintC_ng_mm_d','Maintenance (ng C mm$^{-1}$ d$^{-1}$)','Maintenance C'),('N_ng_mm','Nitrogen (ng N mm$^{-1}$)','Total N')]
COND=[('kr','kr (cm hPa$^{-1}$ d$^{-1}$)'),('Kx','Kx (cm$^4$ hPa$^{-1}$ d$^{-1}$)')]
SUP=str.maketrans('-0123456789','⁻⁰¹²³⁴⁵⁶⁷⁸⁹')
def pf(p):
    if p>=0.001: return f'p = {p:.3f}'
    m,e=f'{p:.1e}'.split('e'); return f'p = {float(m):.1f} × 10'+str(int(e)).translate(SUP)
st='; '.join(f"{r.cost} vs {r.conductivity}: Spearman ρ = {r.spearman_rho:.2f} ({pf(r.spearman_p)}); mixed model slope {r.mixed_slope:.2f} ± {r.mixed_se:.2f} ({pf(r.mixed_p)}), with segment as fixed effect {r.mixed_seg_slope:.2f} ± {r.mixed_seg_se:.2f} ({pf(r.mixed_seg_p)}); within-section r_rm = {r.rm_corr:.2f} (df = {r.rm_df}, {pf(r.rm_p)})" for r in R.itertuples())
segn=P.groupby('segment').agg(r=('anatomy','count'),s=('section','nunique')).reindex(['MP','MM','MT','LP','LM','LT'])
nn=', '.join(f"{g} {int(x.r)} replicates / {int(x.s)} sections" for g,x in segn.iterrows())
LEG=("Figure. Modelled root costs versus hydraulic conductivity, paired replicate by replicate, across tomato (S. pimpinellifolium) root cross-sections.\n\n"
 f"Each point is one replicate: one of Emily Miller's individually measured cells (cell #1–#4 of cortex, endodermis, epidermis and vasculature) within one cross-section, modelled twice from the same measurements — RootSlice for costs and GRANAR + MECHA for conductivity (n = {len(P)} replicates from {P.section.nunique()} cross-sections: {nn}). "
 "x-axes: RootSlice costs per mm of root, whole-root basis — (A, D) construction carbon in cell walls, (B, E) maintenance respiration carbon, (C, F) total nitrogen. "
 "y-axes: MECHA conductivities — (A–C) radial (cross-section) conductivity kr, (D–F) xylem axial conductance Kx. Colour: blue, main root; green, lateral root. Marker: circle, proximal; square, middle; triangle, main tip zone or lateral true tip. All axes log scale.\n\n"
 "Barriers in kr: main proximal and middle, endodermal lignin + exodermal suberin and lignin (MECHA B6); main tip zone, endodermal lignin only (B1); lateral proximal and middle, exodermal suberin and lignin (B10); lateral true tip, none. "
 "Barriers are identical within a cross-section, so within-section relationships are free of barrier effects. The main-root true tip (MT0; MT anatomies without barriers) is not plotted because its costs equal MT; its kr values are in the table.\n\n"
 "Statistics on log10 values. Replicates within a cross-section are not independent, so three analyses are reported: (i) Spearman correlation across all replicates (descriptive); "
 "(ii) linear mixed model log10(conductivity) ~ log10(cost) with cross-section as random intercept (REML), alone and with segment as an additional fixed effect (removes segment and barrier differences; the section variance converged to the boundary in some fits); "
 "(iii) repeated-measures correlation within cross-sections (Bakdash & Marusich 2017; ANCOVA with section as factor). "
 f"Results: {st}.\n\n"
 "Interpretation: radial conductivity falls as cost rises at every level, including between replicates of the same root — replicates with larger cortical cells have a thicker cortex, which costs more and lengthens the radial path. "
 "The weak positive cost–Kx association across roots is carried by the main proximal segment (wider stele and xylem); it turns slightly negative once segment is accounted for and is absent within roots (r_rm ≈ 0).\n\n"
 "Costs: RootSlice (Sidhu et al. 2023; J. Sidhu macOS build), one run per replicate (configs_v5rep), built directly from Emily Miller's sheet (1.032 × 10⁻⁷ mm² px⁻¹; d = 2√(A/π)): cortex, endodermis and epidermis cell = her cell #k, innermost stele cell = her vasculature cell #k (segment median where missing), stele radius from her vasculature area, cortical layers = her count; "
 "cylindrical cells, cell length 100 µm, RCA 0; C = cell-wall volume × 1.51 g cm⁻³ × 0.45; N = 1.835 × 10⁻¹⁵ g µm⁻³ × cytoplasm + 4.202 × 10⁻¹⁵ g µm⁻³ × vacuole; maintenance = cytoplasm × 50 nmol CO₂ s⁻¹ cm⁻³; whole root = cortex / cortex-to-section ratio. "
 "19 replicates with small steles (radius 19–31 µm) hung or crashed in RootSlice with the full stele template; 12 of them were rerun with a simplified stele interior (fewer meta- and protoxylem vessels, smaller inner stele cells; stele radius and all cortical, endodermal and epidermal dimensions unchanged), which left cortical volumes identical in 4 control replicates. Still missing costs: MP34_cell2, MP57_cell1–3, MP59_cell1, MP66_cell2. GRANAR could not tessellate MP19_cell4, MP37_cell2, MP57_cell1–3, MP59_cell1 and MP66_cell1/2. Rows 4 (torn) and 60 (no stele visible) were excluded. "
 "Conductivities: MECHA-dev develop 945b2be, default hydraulic parameters (kw 2.4 × 10⁻⁴ cm² hPa⁻¹ d⁻¹, kAQP 4.3 × 10⁻⁴ cm hPa⁻¹ d⁻¹, Kpl 5.3 × 10⁻¹² cm³ hPa⁻¹ d⁻¹, cell length 200 µm); GRANAR 0.0.1.2 anatomies (v5) built from the same cells as the RootSlice runs, with Emily Miller's cortical layer counts (exodermis counted as one layer), pericycle 21 µm.")
open('Fig_cost_vs_conductivity_replicates_legend.md','w').write(LEG)
def draw(out,leg=None):
    WR=175; H=7.6; LH=0 if leg is None else sum(0.158*len(textwrap.wrap(p,WR))+0.14 for p in leg.split('\n\n'))+0.35
    fig,axs=plt.subplots(2,3,figsize=(13,H+LH)); fig.subplots_adjust(left=0.075,right=0.985,top=1-1.0/(H+LH),bottom=(LH+0.6)/(H+LH),hspace=0.36,wspace=0.3)
    for i,(kv,kl) in enumerate(COND):
        for j,(cv,cl,cn) in enumerate(COST):
            a=axs[i,j]
            for (sec,sg),g in P.groupby(['section','segment']):
                g=g.sort_values(cv); c=COL[sg[0]]
                a.scatter(g[cv],g[kv],c=c,marker=MK[sg[1]],s=30,alpha=0.8,lw=0.5,edgecolors='white',zorder=3)
            a.set_xscale('log'); a.set_yscale('log'); a.grid(color='#eeede8',lw=0.6); a.set_axisbelow(True)
            a.xaxis.set_major_formatter(mt.FuncFormatter(lambda v,_: f'{v:g}')); a.xaxis.set_minor_formatter(mt.NullFormatter())
            a.xaxis.set_major_locator(mt.LogLocator(base=10,subs=(1,2,5)))
            if i==1: a.set_xlabel(cl)
            if j==0: a.set_ylabel(kl)
            a.set_title('ABCDEF'[i*3+j],loc='left',fontweight='bold',fontsize=11)
    h=[Line2D([],[],ls='',marker='o',color=MAIN,label='Main root'),Line2D([],[],ls='',marker='o',color=LAT,label='Lateral root'),
       Line2D([],[],ls='',marker='o',color='#6b6a65',mfc='none',label='Proximal'),Line2D([],[],ls='',marker='s',color='#6b6a65',mfc='none',label='Middle'),
       Line2D([],[],ls='',marker='^',color='#6b6a65',mfc='none',label='Tip')]
    fig.legend(handles=h,loc='upper right',ncol=5,frameon=False,fontsize=8.5,bbox_to_anchor=(0.985,1-0.38/(H+LH)))
    fig.suptitle('Root costs vs hydraulic conductivity, paired by replicate (v5, all sections)',x=0.012,y=1-0.12/(H+LH),ha='left',va='top',fontsize=12.5,fontweight='bold')
    if leg:
        y=(LH-0.1)/(H+LH)
        for k,p in enumerate(leg.split('\n\n')):
            L=textwrap.wrap(p,WR); fig.text(0.012,y,'\n'.join(L),fontsize=8.2,va='top',fontweight='bold' if k==0 else 'normal',color='#1f1f1d',linespacing=1.35); y-=(0.158*len(L)+0.14)/(H+LH)
    fig.savefig(out+'.png',dpi=220); fig.savefig(out+'.pdf'); plt.close(fig)
draw('Fig_cost_vs_conductivity_replicates'); draw('Fig_cost_vs_conductivity_replicates_withLegend',LEG)
