"""Draw every GRANAR anatomy (cells coloured by tissue, barrier walls coloured where MECHA assigns them).
Writes <outputs>/figures/granar/<id>.png (one per anatomy) and granar_gallery_<segment>.pdf / .png
(all replicates of a segment on one page, at a common µm scale so sizes are comparable).
Run after 1_granar/build_anatomies.R:  python figures/plot_granar.py [ids]"""
import os, sys, re, numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from paths import *
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, LineCollection
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from scipy.spatial import ConvexHull

OUTG = os.path.join(FIG, 'granar'); os.makedirs(OUTG, exist_ok=True)
SEG = lambda n: 'MP' if n <= 9 else 'MM' if n <= 20 else 'MT' if n <= 26 else 'LP' if n <= 40 else 'LM' if n <= 50 else 'LT'
ORDER = ['MP', 'MM', 'MT', 'LP', 'LM', 'LT']
NAME = {'MP': 'Main proximal', 'MM': 'Main middle', 'MT': 'Main tip zone', 'LP': 'Lateral proximal', 'LM': 'Lateral middle', 'LT': 'Lateral true tip'}
BAR = {'MP': (1, 1), 'MM': (1, 1), 'MT': (0, 1), 'LP': (1, 0), 'LM': (1, 0), 'LT': (0, 0)}   # (exodermis, endodermis) barrier in MECHA
CODE = {'MP': 'B6', 'MM': 'B6', 'MT': 'B1', 'LP': 'B10', 'LM': 'B10', 'LT': 'none'}
GROUP = {1: 'exodermis', 2: 'epidermis', 3: 'endodermis', 4: 'cortex', 5: 'stele', 11: 'phloem', 12: 'companion_cell', 13: 'xylem', 16: 'pericycle'}
FILL = {'epidermis': '#e9dcc6', 'exodermis': '#f0d6a0', 'cortex': '#f1e9d6', 'endodermis': '#f4cdc4', 'pericycle': '#e2dcee', 'stele': '#dad5ea',
        'phloem': '#c4d9b2', 'companion_cell': '#c4d9b2', 'xylem': '#34568c'}
EXO, ENDO = '#eda100', '#e34948'

def parse(path):
    t = open(path).read()
    cells_txt = t.split('<cells')[1].split('</cells>')[0]
    cells = [(int(g), [int(w) for w in re.findall(r'<wall id="\s*(\d+)"', body)])
             for g, body in re.findall(r'<cell id="\s*\d+" group="(\d+)"[^>]*>(.*?)</cell>', cells_txt, re.S)]
    walls = {int(wid): [(float(x), float(y)) for x, y in re.findall(r'x="\s*([\d.eE+-]+)" y="\s*([\d.eE+-]+)"', body)]
             for wid, body in re.findall(r'<wall id="\s*(\d+)"[^>]*>(.*?)</wall>', t.split('</cells>')[1], re.S)}
    return cells, walls

def draw(ax, aid, seg, scale_um=None):
    cells, walls = parse(os.path.join(ANAT, aid + '_root.xml'))
    pts = np.array([p for w in walls.values() for p in w]) * 1000          # mm -> µm
    cen = np.median(pts, 0); polys, cols, exo, endo = [], [], set(), set()
    for g, ws in cells:
        q = np.array([p for w in ws for p in walls.get(w, [])]) * 1000 - cen
        if len(q) < 3: continue
        try: h = ConvexHull(q)
        except Exception: continue
        tp = GROUP.get(g, 'cortex'); polys.append(q[h.vertices]); cols.append(FILL.get(tp, '#eeeeee'))
        if tp == 'exodermis': exo.update(ws)
        if tp == 'endodermis': endo.update(ws)
    ax.add_collection(PolyCollection(polys, facecolors=cols, edgecolors='#7a7066', linewidths=0.2))
    e, n = BAR[seg]
    for on, ws, c in ((e, exo, EXO), (n, endo, ENDO)):
        if on: ax.add_collection(LineCollection([np.array(walls[w]) * 1000 - cen for w in ws if len(walls.get(w, [])) >= 2], colors=c, linewidths=0.9))
    R = scale_um or np.abs(np.vstack(polys)).max() * 1.03
    ax.set_xlim(-R, R); ax.set_ylim(-R, R); ax.set_aspect('equal'); ax.axis('off')
    return R

def legend(fig, y=0.012):
    h = [Patch(fc=FILL[k], ec='#7a7066', label=k.replace('_', ' ')) for k in ('epidermis', 'exodermis', 'cortex', 'endodermis', 'pericycle', 'stele', 'phloem', 'xylem')]
    h += [Line2D([], [], color=EXO, lw=2.5, label='exodermis with barrier in MECHA (suberin + lignin)'), Line2D([], [], color=ENDO, lw=2.5, label='endodermis with Casparian strip in MECHA')]
    fig.legend(handles=h, loc='lower center', ncol=5, fontsize=8, frameon=False, bbox_to_anchor=(0.5, y))

if __name__ == '__main__':
    S = pd.read_csv(GRANAR_SUMMARY); S = S[S.ok]
    ids = sys.argv[1:] or sorted(S.id, key=lambda a: (int(a.split('_')[0][2:]), a))
    seg_of = lambda a: SEG(int(a.split('_')[0][2:]))
    for a in ids:                                                            # one image per anatomy
        sg = seg_of(a); fig, ax = plt.subplots(figsize=(4.2, 4.6)); R = draw(ax, a, sg)
        ax.plot([-R * 0.95, -R * 0.95 + 100], [-R * 0.97] * 2, color='k', lw=2); ax.text(-R * 0.95 + 50, -R * 0.93, '100 µm', ha='center', fontsize=8)
        ax.set_title(f'{a}  ·  {NAME[sg]} ({sg})  ·  MECHA {CODE[sg]}', fontsize=9)
        fig.savefig(os.path.join(OUTG, a + '.png'), dpi=200, bbox_inches='tight'); plt.close(fig)
    if sys.argv[1:]: sys.exit()
    Rmax = {}                                                                # galleries: one page per segment, common scale per page
    for sg in ORDER:
        A = [a for a in ids if seg_of(a) == sg]; nc = 6; nr = int(np.ceil(len(A) / nc))
        fig, axs = plt.subplots(nr, nc, figsize=(nc * 1.9, nr * 2.0 + 1.0)); axs = np.atleast_2d(axs)
        R = max(np.sqrt(S.set_index('id').loc[A, 'area_all'] / np.pi).max() * 1000 * 1.15, 1)
        for ax, a in zip(axs.flat, A): draw(ax, a, sg, R); ax.set_title(a, fontsize=7)
        for ax in axs.flat[len(A):]: ax.axis('off')
        ax0 = axs.flat[0]; ax0.plot([-R * 0.95, -R * 0.95 + 100], [-R * 0.98] * 2, color='k', lw=1.5); ax0.text(-R * 0.95 + 50, -R * 0.9, '100 µm', ha='center', fontsize=6)
        fig.suptitle(f'GRANAR anatomies: {NAME[sg]} ({sg}), n = {len(A)}, MECHA barrier {CODE[sg]}. All panels on the same scale.', fontsize=10)
        legend(fig, 0.0); fig.tight_layout(rect=(0, 0.05, 1, 0.97))
        for ext in ('pdf', 'png'): fig.savefig(os.path.join(FIG, f'granar_gallery_{sg}.{ext}'), dpi=150)
        plt.close(fig)
    print(len(ids), 'anatomies drawn into', OUTG, '+ 6 segment galleries in', FIG)
