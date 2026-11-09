"""Small-stele fallback configs for the 19 replicates that hung or crashed in RootSlice (all with stele radius 19-31 um).
Geometry (stele radius, cortex cell, layers, epidermis, endodermis) is unchanged; only the stele interior template is simplified,
which does not enter the cortical cost outputs. Three variants, tried in order by 3_rootslice/run_small_stele.sh:
  A  MXNum 3, MXNumEachRing 2, SteleCellNumBetween 3, PXNum 3, protoxylem SurroundingCellRingNum 1   (template: MXNum 5, 4 per ring, 5 between, PXNum 4, ring 2)
  B  A + innermost stele cell diameter x 0.7
  C  MXNum 1, MXNumEachRing 1, SteleCellNumBetween 2, PXNum 2, protoxylem ring 1, innermost stele cell x 0.7
Controls: variant A on 4 replicates that ran with the full template (MP34_cell3, MP35_cell1, MP45_cell1, MP65_cell1),
to check that cortical volumes do not change.
Writes 3_rootslice/configs_smallstele/<id>_<variant>.xml and CTRL_<id>_A.xml."""
import os, sys; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from paths import *
import os, re
SRC = RS_CONFIGS; SDIR = RS_CONFIGS_SS; os.makedirs(SDIR, exist_ok=True)
FAILED = ['MP32_cell1', 'MP32_cell2', 'MP32_cell3', 'MP34_cell1', 'MP34_cell2', 'MP34_cell4', 'MP43_cell3', 'MP46_cell1', 'MP46_cell2', 'MP49_cell2',
          'MP52_cell1', 'MP57_cell1', 'MP57_cell2', 'MP57_cell3', 'MP59_cell1', 'MP66_cell1', 'MP66_cell2', 'MP66_cell3', 'MP66_cell4']
CTRL = ['MP34_cell3', 'MP35_cell1', 'MP45_cell1', 'MP65_cell1']
block = lambda s, t: re.search(rf"<{t}\b.*?</{t}>", s, re.S).group(0)
def put(s, blk, tag, v):
    b = block(s, blk); m = re.search(rf"(<{tag}>)([^<]+)(<)", b); assert m, (blk, tag)
    return s.replace(b, b.replace(m.group(0), f"{m.group(1)}{v:g}{m.group(3)}", 1))
def get(s, blk, tag): return float(re.search(rf"<{tag}>([^<]+)<", block(s, blk)).group(1))
def variant(s, v):
    if v in 'AB':
        s = put(s, 'MetaXylem', 'MXNum', 3); s = put(s, 'MetaXylem', 'MXNumEachRing', 2); s = put(s, 'MetaXylem', 'SteleCellNumBetween', 3)
        s = put(s, 'ProtoXylem', 'PXNum', 3); s = put(s, 'ProtoXylem', 'SurroundingCellRingNum', 1)
    if v == 'C':
        s = put(s, 'MetaXylem', 'MXNum', 1); s = put(s, 'MetaXylem', 'MXNumEachRing', 1); s = put(s, 'MetaXylem', 'SteleCellNumBetween', 2)
        s = put(s, 'ProtoXylem', 'PXNum', 2); s = put(s, 'ProtoXylem', 'SurroundingCellRingNum', 1)
    if v in 'BC': s = put(s, 'Stele', 'InnestCellDiameer', round(get(s, 'Stele', 'InnestCellDiameer') * 0.7, 3))
    return s
n = 0
for a in FAILED:
    T = open(os.path.join(SRC, f'R5_{a}.xml')).read()
    for v in 'ABC': open(os.path.join(SDIR, f'{a}_{v}.xml'), 'w').write(variant(T, v)); n += 1
for a in CTRL:
    T = open(os.path.join(SRC, f'R5_{a}.xml')).read().replace(f'<Main name="R5_{a}"', f'<Main name="CTRL_{a}"')
    open(os.path.join(SDIR, f'CTRL_{a}_A.xml'), 'w').write(variant(T, 'A')); n += 1
print(n, 'configs in 3_rootslice/configs_smallstele/')
