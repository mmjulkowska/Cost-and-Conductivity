"""Per-replicate RootSlice configs (all usable sections of Emily Miller's sheet, Nov 2026).
One config per v5 GRANAR input (section x measured cell k), using exactly the same values as the GRANAR/MECHA anatomy
(1_granar/v5_replicates.csv, in mm):
  cortex CellDiameter = Emily cortex cell #k, Epidermis/Endodermis = her cell #k, Stele InnestCellDiameer = vasculature cell #k
  (segment median if not measured), stele radius = vasculature diameter / 2, NumFiles = her cortical layer count
  (exodermis counted as a layer, as in GRANAR), Thickness = cell x layers.
All other fields: template_rootslice.xml (RCA 0, cell length 100 um, Jagdeep's settings for xylem, exodermis, sclerenchyma).
Includes the 8 replicates GRANAR could not tessellate (they get costs but no conductivities).
Writes 3_rootslice/configs/R5_MP<n>_cell<k>.xml and 3_rootslice/configs_values.csv."""
import os, sys; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from paths import *
import os, re, csv
CDIR=RS_CONFIGS; os.makedirs(CDIR,exist_ok=True)
T=open(RS_TEMPLATE).read()
block=lambda s,t: re.search(rf"<{t}\b.*?</{t}>",s,re.S).group(0)
def put(s,blk,tag,v):
    b=block(s,blk); m=re.search(rf"(<{tag}>)([^<]+)(<)",b); return s.replace(b,b.replace(m.group(0),f"{m.group(1)}{v:g}{m.group(3)}",1))
rows=[]
for r in csv.DictReader(open(REPLICATES)):
    if r["used"]!="True": continue
    a=r["id"]; um=lambda c: round(float(r[c])*1000,3); nl=int(float(r["emily_layers"]))
    cd,ep,en,sc,sr=um("cortex_mm"),um("epi_mm"),um("endo_mm"),um("stele_cell_mm"),um("stele_layer_mm")/2
    s=re.sub(r'<Main name="[^"]*"',f'<Main name="R5_{a}"',T,count=1)
    s=put(s,"Cortical","BaseRadius",sr); s=put(s,"Cortical","CellDiameter",cd); s=put(s,"Cortical","NumFiles",nl); s=put(s,"Cortical","Thickness",round(cd*nl,2))
    s=put(s,"Epidermis","CellDiameter",ep); s=put(s,"Endodermis","CellDiameter",en); s=put(s,"Stele","InnestCellDiameer",sc)
    open(os.path.join(CDIR,f"R5_{a}.xml"),"w").write(s)
    rows.append(dict(anatomy=a,section=int(r["MPn"]),segment=r["segment"],stele_r=sr,cortex_cell=cd,layers=nl,epidermis=ep,endodermis=en,stele_cell=sc,root_r_approx=round(sr+en+cd*nl+ep,1)))
w=csv.DictWriter(open(RS_VALUES,"w",newline=""),fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print(len(rows),"configs in 3_rootslice/configs/")
