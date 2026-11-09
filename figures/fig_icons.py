"""Conductivity figure with half-photo / half-GRANAR cross-section icons (run figures/make_icons.py and 2_mecha/analyse.py first)."""
import os, sys; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from paths import *
import sys, json, pandas as pd
sys.path.insert(0,os.path.join(ROOT,'shared')); from fig_conductivity_icons import draw
info=json.load(open(os.path.join(ICONS,'icons_info.json')))
L=open(os.path.join(FIG,'Fig_conductivity_barriers_v5_legend.md')).read()
icons=("Cross-section icons in (A): the left half is a representative section photographed by Emily Miller (stained with toluidine blue; chosen among sections judged usable for clarity and a root diameter near the segment median), "
 "the right half is a GRANAR v5 anatomy built from that same section's measurements (the replicate whose kr is closest to the segment median). Circle size is proportional to the photographed root diameter (scale bar 200 µm); within a circle, each half is scaled to its own root diameter so both fill it; "
 "the GRANAR roots are 1.0-1.4x wider than the photographed ones (measured / modelled diameter, µm: "
 + "; ".join(f"{k} {v['emily_root_d_um']} / {v['granar_root_d_um']} ({v['emily_row'] and 'row '+str(v['emily_row'])}, {v['anatomy']})" for k,v in info.items() if k!='MT0')
 + "). MT0 shows the MT section. GRANAR colours: blue, xylem vessel; green, phloem; lilac, stele and pericycle; pink, endodermis; yellow, exodermis; cream, cortex; tan, epidermis. "
 "Barrier walls are drawn in orange (exodermal suberin + lignin) and red (endodermal lignin, Casparian strip) only where that barrier is assigned in MECHA.")
paras=L.split('\n\n'); paras.insert(2,icons); L2='\n\n'.join(paras)
open(os.path.join(FIG,'Fig_conductivity_barriers_v5_icons_legend.md'),'w').write(L2)
lk=pd.read_csv(os.path.join(STATS,'kr_summary_letters.csv')).set_index('fragment').tukey_letter.to_dict()
lx=pd.read_csv(os.path.join(STATS,'Kx_summary_letters.csv')).set_index('fragment').tukey_letter.to_dict()
ttl="Root hydraulic conductivity by segment and barrier"
D={k:v['emily_root_d_um'] for k,v in info.items()}
draw(COND,lk,lx,ttl,os.path.join(FIG,'Fig_conductivity_barriers_v5_icons'),icon_dir=ICONS,diam_um=D)
draw(COND,lk,lx,ttl,os.path.join(FIG,'Fig_conductivity_barriers_v5_icons_withLegend'),legend_text=L2,icon_dir=ICONS,diam_um=D)
print('ok')
