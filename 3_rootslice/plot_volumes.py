"""Readable versions of RootSlice's per-run bar plots (the build's own bar_*.png put 200+ run names on one 15 cm axis).
For each metric from the RootSlice volume analysis, writes to <outputs>/rootslice/figures/:
  <metric>_by_segment.png/.pdf   one point per replicate, grouped by segment, with box plots
  <metric>_replicates.png/.pdf   one horizontal bar per replicate (every label legible), grouped by segment,
                                 coloured by the replicate's measured cortex cell diameter
Run after 3_rootslice/run_rootslice.sh and run_small_stele.sh:  python 3_rootslice/plot_volumes.py"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from paths import *
import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm, colors

OUTF = os.path.join(RS_OUT, 'figures'); os.makedirs(OUTF, exist_ok=True)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42})
ORDER = ['MP', 'MM', 'MT', 'LP', 'LM', 'LT']
NAME = {'MP': 'Main proximal', 'MM': 'Main middle', 'MT': 'Main tip zone', 'LP': 'Lateral proximal', 'LM': 'Lateral middle', 'LT': 'Lateral true tip'}
SEGC = {'MP': '#1f4e8c', 'MM': '#2a78d6', 'MT': '#7fb2ee', 'LP': '#117a55', 'LM': '#1baf7a', 'LT': '#86d9b6'}
METRICS = [('cytoplasm_pct', 'Cytoplasm volume per cortex volume (%)', 'cytoplasm'),
           ('vacuole_pct', 'Vacuole volume per cortex volume (%)', 'vacuole'),
           ('Mito_per_cm3_cortex', 'Mitochondria per cm³ cortex', 'mitochondria'),
           ('Resp_umol_per_s_per_cm3', 'Respiration (µmol CO₂ cm⁻³ cortex s⁻¹)', 'respiration'),
           ('N_g_per_cm3_cortex', 'Nitrogen (g cm⁻³ cortex)', 'nitrogen'),
           ('P_g_per_cm3_cortex', 'Phosphorus (g cm⁻³ cortex)', 'phosphorus')]
SEG = lambda n: 'MP' if n <= 9 else 'MM' if n <= 20 else 'MT' if n <= 26 else 'LP' if n <= 40 else 'LM' if n <= 50 else 'LT'

v = pd.read_csv(RS_VOL)
if os.path.exists(RS_VOL_SS):
    s = pd.read_csv(RS_VOL_SS); v = pd.concat([v, s[s.model.str.startswith('R5_') & ~s.model.isin(v.model)]], ignore_index=True)
v['id'] = v.model.str.replace('R5_', '', regex=False)
v['section'] = v.id.str.extract(r'MP(\d+)_')[0].astype(int); v['cell'] = v.id.str.extract(r'cell(\d+)')[0].astype(int)
v['segment'] = v.section.map(SEG)
v = v.merge(pd.read_csv(RS_VALUES)[['anatomy', 'cortex_cell']].rename(columns={'anatomy': 'id'}), on='id', how='left')
v['o'] = v.segment.map(ORDER.index); v = v.sort_values(['o', 'section', 'cell']).reset_index(drop=True)
norm = colors.Normalize(v.cortex_cell.min(), v.cortex_cell.max()); cmap = matplotlib.colormaps['Blues']

for col, lab, short in METRICS:
    if col not in v: continue
    # ---- summary by segment
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    data = [v.loc[v.segment == g, col].values for g in ORDER]
    ax.boxplot(data, positions=range(6), widths=0.55, showfliers=False, medianprops=dict(color='k'), boxprops=dict(color='#666'), whiskerprops=dict(color='#666'), capprops=dict(color='#666'))
    rng = np.random.default_rng(1)
    for i, g in enumerate(ORDER):
        y = data[i]; ax.scatter(i + rng.uniform(-0.18, 0.18, len(y)), y, s=14, color=SEGC[g], alpha=0.8, lw=0, zorder=3)
    ax.set_xticks(range(6)); ax.set_xticklabels([f'{g}\n{NAME[g]}\n(n = {len(d)})' for g, d in zip(ORDER, data)], fontsize=9)
    ax.set_ylabel(lab); ax.ticklabel_format(axis='y', style='sci', scilimits=(-3, 4))
    ax.set_title(f'RootSlice: {short} by segment (one point per replicate)', fontsize=11)
    fig.tight_layout()
    for ext in ('png', 'pdf'): fig.savefig(os.path.join(OUTF, f'{short}_by_segment.{ext}'), dpi=200)
    plt.close(fig)
    # ---- every replicate, horizontal bars
    gap = 1.5; y = []; pos = 0.0; prev = None; seg_mid = {}
    for g in v.segment:
        if prev is not None and g != prev: pos += gap
        y.append(pos); seg_mid.setdefault(g, []).append(pos); pos += 1; prev = g
    y = np.array(y); H = 0.16 * len(v) + 1.8
    fig, ax = plt.subplots(figsize=(8.5, H))
    ax.barh(y, v[col], height=0.8, color=cmap(norm(v.cortex_cell)), edgecolor='#333', lw=0.3)
    ax.set_yticks(y); ax.set_yticklabels(v.id, fontsize=7); ax.invert_yaxis(); ax.set_ylim(y.max() + 1, -1)
    xmax = v[col].max()
    for g, ps in seg_mid.items():
        ax.axhspan(min(ps) - 0.6, max(ps) + 0.6, color=SEGC[g], alpha=0.07, lw=0)
        ax.text(xmax * 1.02, np.mean(ps), f'{g}\n{NAME[g]}', va='center', ha='left', fontsize=9, color=SEGC[g], fontweight='bold')
    ax.set_xlim(0, xmax * 1.25); ax.set_xlabel(lab); ax.ticklabel_format(axis='x', style='sci', scilimits=(-3, 4))
    ax.xaxis.set_ticks_position('both'); ax.tick_params(axis='x', labeltop=True)
    sm = cm.ScalarMappable(norm=norm, cmap=cmap); cb = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.01, aspect=60)
    cb.set_label('Measured cortex cell diameter (µm)', fontsize=9)
    ax.set_title(f'RootSlice: {short}, every replicate ({len(v)}), grouped by segment', fontsize=11, pad=24)
    fig.tight_layout()
    for ext in ('png', 'pdf'): fig.savefig(os.path.join(OUTF, f'{short}_replicates.{ext}'), dpi=150)
    plt.close(fig)
print(len(v), 'replicates;', len(METRICS) * 2, 'figures in', OUTF)
