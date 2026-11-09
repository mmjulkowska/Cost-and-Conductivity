"""Builds the workbook for the Pareto model (Arjun Chandrasekhar) into <outputs>/deliverables/."""
import os, sys; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from paths import *
os.chdir(COSTS)  # cost tables live in <outputs>/costs/
import pandas as pd, numpy as np
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter as L
P=pd.read_csv('replicate_level_costs_conductivity.csv'); R=pd.read_csv('correlations_replicate_level.csv')
TK={'kr':pd.read_csv(os.path.join(STATS, 'kr_summary_letters.csv')).set_index('fragment').tukey_letter.to_dict(),
    'Kx':pd.read_csv(os.path.join(STATS, 'Kx_summary_letters.csv')).set_index('fragment').tukey_letter.to_dict()}
TKc=pd.read_csv('costs_C_ng_mm_summary_letters.csv').set_index('fragment').tukey_letter.to_dict()
META={'MP':('Main proximal','3 cm above MM','lignin','lignin + suberin','B6'),'MM':('Main middle','3 cm above MT','lignin','lignin + suberin','B6'),
      'MT':('Main tip zone','3 cm from apex','lignin','none','B1'),'MT0':('Main true tip','apex (MT anatomy, barriers removed)','none','none','none'),
      'LP':('Lateral proximal','3 cm above LM','none','lignin + suberin','B10'),'LM':('Lateral middle','3 cm above LT','none','lignin + suberin','B10'),'LT':('Lateral true tip','apex','none','none','none')}
ORDER=['MP','MM','MT','MT0','LP','LM','LT']
F=Font(name='Arial',size=10); FB=Font(name='Arial',size=10,bold=True); FH=Font(name='Arial',size=10,bold=True,color='FFFFFF'); FT=Font(name='Arial',size=13,bold=True); FN=Font(name='Arial',size=9,italic=True,color='555555')
HF=PatternFill('solid',fgColor='2A78D6'); thin=Side(style='thin',color='CCCCCC')
def header(ws,row,cols):
    for j,c in enumerate(cols,1):
        x=ws.cell(row,j,c); x.font=FH; x.fill=HF; x.alignment=Alignment(wrap_text=True,vertical='center',horizontal='center')
    ws.row_dimensions[row].height=44
wb=Workbook()
# ---------------- Replicates ----------------
ws=wb.active; ws.title='Replicates'
cols=['Fragment','Section','Anatomy (replicate)','Endodermis barrier','Exodermis barrier','MECHA barrier',
      'Building C (ng C/mm)','Maintenance C (ng C/mm/d)','Cytoplasmic N (ng N/mm)','Vacuolar N (ng N/mm)','Total N (ng N/mm)',
      'kr radial (cm/hPa/d)','Kx axial (cm4/hPa/d)','Stele radius (um)','Cortex cell diameter (um)','Cortical layers','Epidermis cell (um)','Endodermis cell (um)']
ws['A1']='Replicate-level paired costs (RootSlice) and conductivities (MECHA) — one row per measured cortex cell replicate'; ws['A1'].font=FT
header(ws,3,cols); r0=4; rows=[]
for frag in ORDER:
    src='MT' if frag=='MT0' else frag
    for _,p in P[P.segment==src].iterrows():
        kr=p.kr_trueTip_MT0 if frag=='MT0' else p.kr
        rows.append([frag,f"MP{p.section}",p.anatomy,META[frag][2],META[frag][3],META[frag][4],p.C_ng_mm,p.maintC_ng_mm_d,p.Ncyt_ng_mm,p.Nvac_ng_mm,None,kr,p.Kx,p.stele_r,p.cortex_cell,p.layers,p.epidermis,p.endodermis])
for i,rw in enumerate(rows):
    r=r0+i
    for j,v in enumerate(rw,1):
        c=ws.cell(r,j,v); c.font=F
    ws.cell(r,11,f'=I{r}+J{r}').font=F
    for j in (7,8,9,10,11): ws.cell(r,j).number_format='#,##0.0'
    for j in (12,13): ws.cell(r,j).number_format='0.000E+00'
last=r0+len(rows)-1
for j,w in enumerate([9,9,16,12,16,9,12,13,13,13,12,13,13,11,12,9,11,11],1): ws.column_dimensions[L(j)].width=w
ws.freeze_panes='D4'
ws.cell(last+2,1,'MT0 rows repeat the MT replicates (same anatomy and costs) with all apoplastic barriers removed in MECHA; only kr differs. Absent: 6 replicates RootSlice could not run even with a simplified stele (MP34_cell2, MP57_cell1–3, MP59_cell1, MP66_cell2; 12 other small-stele replicates were recovered that way) and 8 GRANAR failures (MP19_cell4, MP37_cell2, MP57_cell1–3, MP59_cell1, MP66_cell1/2).').font=FN
# ---------------- Segment summary ----------------
S=wb.create_sheet('Segment summary',0)
S['A1']='Root segment costs and hydraulic conductivities — tomato (S. pimpinellifolium), for the Pareto-front model'; S['A1'].font=FT
S['A2']='Means and SDs are live formulas over the Replicates sheet. Tukey letters (hardcoded) come from one-way ANOVA on log10 values; groups sharing a letter do not differ (α = 0.05).'; S['A2'].font=FN
scol=['Fragment','Segment','Position','Endodermis','Exodermis','MECHA barrier','Sections (n)','Replicates (n)',
      'Building C mean (ng C/mm)','Building C SD','Building C Tukey','Maintenance mean (ng C/mm/d)','Maintenance SD','Total N mean (ng N/mm)','Total N SD',
      'kr mean (cm/hPa/d)','kr SD','kr Tukey','Kx mean (cm4/hPa/d)','Kx SD','Kx Tukey','C payback (d) = C / maintenance','Building C relative to LT','kr relative to LT']
header(S,4,scol)
RG=lambda col: f"Replicates!${col}${r0}:${col}${last}"
rowof={}
for i,frag in enumerate(ORDER):
    r=5+i; rowof[frag]=r; m=META[frag]
    vals=[frag,m[0],m[1],m[2],m[3],m[4]]
    for j,v in enumerate(vals,1): S.cell(r,j,v).font=F
    S.cell(r,7,f'=SUMPRODUCT(({RG("A")}=A{r})/COUNTIFS({RG("A")},{RG("A")},{RG("B")},{RG("B")}))')
    S.cell(r,8,f'=COUNTIF({RG("A")},A{r})')
    def mean(col): return f'=AVERAGEIF({RG("A")},$A{r},{RG(col)})'
    def sd(col,mc): return f'=SQRT(SUMPRODUCT(({RG("A")}=$A{r})*({RG(col)}-{mc}{r})^2)/(H{r}-1))'
    S.cell(r,9,mean('G')); S.cell(r,10,sd('G','I')); S.cell(r,11,TKc.get('MT' if frag=='MT0' else frag,'') + (' (= MT)' if frag=='MT0' else ''))
    S.cell(r,12,mean('H')); S.cell(r,13,sd('H','L'))
    S.cell(r,14,mean('K')); S.cell(r,15,sd('K','N'))
    S.cell(r,16,mean('L')); S.cell(r,17,sd('L','P')); S.cell(r,18,TK['kr'][frag])
    S.cell(r,19,mean('M')); S.cell(r,20,sd('M','S')); S.cell(r,21,TK['Kx'].get(frag,TK['Kx']['MT']+' (= MT)'))
    S.cell(r,22,f'=I{r}/L{r}')
    for j in range(7,25): S.cell(r,j).font=F
    for j in (9,10,12,13,14,15): S.cell(r,j).number_format='#,##0.0'
    for j in (16,17,19,20): S.cell(r,j).number_format='0.00E+00'
    S.cell(r,22).number_format='0.0'
for frag in ORDER:
    r=rowof[frag]; S.cell(r,23,f'=I{r}/$I${rowof["LT"]}').number_format='0.00'; S.cell(r,24,f'=P{r}/$P${rowof["LT"]}').number_format='0.00'
    S.cell(r,23).font=F; S.cell(r,24).font=F
for j,w in enumerate([8,16,20,10,15,9,9,10,12,10,9,13,11,12,10,12,10,8,12,10,9,13,11,10],1): S.column_dimensions[L(j)].width=w
S.freeze_panes='C5'
notes=["Tukey letters: building C from section means of the replicate RootSlice runs (n = sections; maintenance C has the same letters; total N differs only in MT = b); kr and Kx from the replicate-level MECHA v5 run (n = anatomies). MT0 costs and Kx are MT's by construction and were not tested separately.",
       "Costs are per mm of root on a whole-root basis (RootSlice cortex values divided by the cortex-to-section area ratio). kr = cross-section radial conductivity; Kx = xylem axial conductance (MECHA).",
       "Building C = cell-wall C; maintenance = respiration C from cytoplasm volume (50 nmol CO2 /s /cm3 cytoplasm); N = cytoplasmic + vacuolar N. Sources and versions: see Notes sheet."]
for k,t in enumerate(notes): c=S.cell(13+k,1,t); c.font=FN
# ---------------- Correlations ----------------
C=wb.create_sheet('Correlations')
C['A1']='Cost vs conductivity, paired by replicate (log10 values, n = '+str(len(P))+' replicates from '+str(P.section.nunique())+' cross-sections)'; C['A1'].font=FT
ccols=['Cost','Conductivity','Replicates','Sections','Spearman rho','Spearman p','Mixed model slope (section random)','SE','p',
       'Mixed slope + segment fixed effect','SE','p','Within-section r (repeated measures)','df','p']
header(C,3,ccols)
for i,r in enumerate(R.itertuples()):
    vals=[r.cost,r.conductivity,r.n_replicates,r.n_sections,r.spearman_rho,r.spearman_p,r.mixed_slope,r.mixed_se,r.mixed_p,r.mixed_seg_slope,r.mixed_seg_se,r.mixed_seg_p,r.rm_corr,r.rm_df,r.rm_p]
    for j,v in enumerate(vals,1):
        c=C.cell(4+i,j,v); c.font=F
        if j in (5,7,8,10,11,13): c.number_format='0.000'
        if j in (6,9,12,15): c.number_format='0.0E+00'
for j,w in enumerate([14,12,10,9,10,10,14,8,9,14,8,9,14,6,9],1): C.column_dimensions[L(j)].width=w
for k,t in enumerate(["Mixed model: log10(conductivity) ~ log10(cost) + (1 | cross-section), REML (statsmodels). Repeated-measures correlation: ANCOVA with cross-section as factor (Bakdash & Marusich 2017).",
    "Radial kr falls with cost at every level, including within a cross-section (barriers are identical within a section, so this is anatomical). The weak positive cost–Kx association across roots is carried by the main proximal segment; with segment included it turns slightly negative, and within cross-sections it is absent.",
    "MT0 rows are excluded (their costs duplicate MT). Values are computed in Python (scripts in the GitHub repository mmjulkowska/Cost-and-Conductivity, folder 4_costs/), not live formulas."]):
    C.cell(12+k,1,t).font=FN
# ---------------- Notes ----------------
N=wb.create_sheet('Notes')
lines=[('Notes, sources and versions',FT),('',F),
 ('Version v5, prepared 7 Nov 2026 for Arjun Chandrasekhar (Pareto-front root architecture model) — Julkowska lab, BTI.',F),
 ('Anatomy data: Emily Miller, root cross-sections of S. pimpinellifolium (pixel areas × 1.032e-7 mm² per px, from her 0.1 mm stage micrometer). v5 uses every usable section of her sheet (55 sections; rows 4 torn and 60 no stele excluded). Replicate k = her cell #k for cortex, endodermis, epidermis and vasculature (segment median where missing); stele from her vasculature area; pericycle 21 µm.',F),
 ('Segments: positions 3 cm apart along each root (main: MP–MM–MT–apex; lateral: LP–LM–apex). MT0 = main root apex, modelled with the MT anatomy and no barriers.',F),
 ('Barriers (Magda Julkowska, 6 Nov 2026, from Emily Miller\'s staining): main proximal/middle endodermal lignin + exodermal suberin & lignin (B6); main tip zone endodermal lignin only (B1); lateral proximal/middle exodermal suberin & lignin (B10); main and lateral apex none (all barrier walls set to primary-wall kw, after V. Couvreur).',F),
 ('Costs: RootSlice (Sidhu et al. 2023), J. Sidhu macOS build; cylindrical cells, 100 µm cell length, no aerenchyma; geometry per replicate (stele radius, cortex cell, layer count = Emily\'s, epidermis, endodermis). C = wall volume × 1.51 g/cm3 × 0.45; N = 1.835e-15 g/µm3 × cytoplasm + 4.202e-15 g/µm3 × vacuole; maintenance = cytoplasm × 50 nmol CO2/s/cm3.',F),
 ('Conductivities: MECHA (openalea.mecha, MECHA-dev develop 945b2be + one-line merge fix), default hydraulics (kw 2.4e-4 cm2/hPa/d, kAQP 4.3e-4 cm/hPa/d, Kpl 5.3e-12 cm3/hPa/d, cell length 200 µm); GRANAR 0.0.1.2 anatomies built from the same cells as RootSlice, with Emily\'s cortex layer counts (exodermis counted as one layer).',F),
 ('Caveats: (1) whole-root costs assume non-cortical tissue costs the same per area as cortex; (2) MT0 Kx assumes mature xylem at the apex (upper bound); (3) 13 small-stele replicates (stele radius < 31 µm) were run in RootSlice with a simplified stele interior; cortical volumes, which give the costs, were unchanged by this in controls; (4) replicates within a section are not independent — use section means or mixed models for inference.',F),
 ('Code and model inputs: GitHub mmjulkowska/Cost-and-Conductivity (1_granar/ GRANAR inputs, 3_rootslice/configs/ RootSlice inputs). Outputs: Box > Plant.Architecture.Lab > Magda > NSF Math Bio - ParetoFront > 20261109_GRANAR_MECHA_RootSlice > outputs.',F)]
for i,(t,f) in enumerate(lines,1): c=N.cell(i,1,t); c.font=f; c.alignment=Alignment(wrap_text=True,vertical='top')
N.column_dimensions['A'].width=150
# ---------------- kr uncertainty (Monte Carlo) ----------------
import os as _os
if _os.path.exists(os.path.join(MC, 'mc_segment_means.csv')):
    U=wb.create_sheet('kr uncertainty'); MC=pd.read_csv(os.path.join(MC, 'mc_segment_means.csv'))
    U['A1']='Uncertainty of segment kr from cell hydraulic properties (Monte Carlo, 200 draws)'; U['A1'].font=FT
    U['A2']='kw, kmb, kAQP and Kpl each varied independently between default/2 and default x2 (log-uniform, Latin hypercube); all v5 anatomies with assigned barriers. Kx is not affected.'; U['A2'].font=FN
    header(U,4,['Fragment','kr at defaults (mean)','kr 2.5% (MC)','kr 97.5% (MC)','Fold range','kr relative to LT, 2.5%','kr relative to LT, 97.5%'])
    for i,r in enumerate(MC.itertuples()):
        for j,v in enumerate([r.fragment,r.default_mean,r.mc_q025,r.mc_q975,r.fold_95,r.rel_LT_q025,r.rel_LT_q975],1):
            c=U.cell(5+i,j,v); c.font=F
            if j in (2,3,4): c.number_format='0.00E+00'
            if j in (5,6,7): c.number_format='0.00'
    for j,w in enumerate([10,16,14,14,10,16,16],1): U.column_dimensions[L(j)].width=w
    for k,t in enumerate(['Robust in every draw: MP and MM lowest; MT0 and LT highest (MT, LP, LM < MT0 < LT); negative cost-kr relationship within sections.',
                          'Not robust: MT (endodermal barrier, B1) vs LP/LM (exodermal barrier, B10) - their order depends on the aquaporin vs plasmodesmata balance (MT above LP/LM in 45% of draws).',
                          'Details: outputs/mecha/sensitivity_mc/ (Fig_kr_sensitivity_MC_withLegend); code in 2_mecha/sensitivity_mc/.']):
        U.cell(14+k,1,t).font=FN
wb.save(os.path.join(DELIV, 'Root_costs_conductivity_for_Arjun_v5_Nov2026.xlsx')); print('saved', last-r0+1,'replicate rows')
