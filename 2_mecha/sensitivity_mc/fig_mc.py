"""Figure for the kr Monte Carlo (run analyse_mc.py first). Writes Fig_kr_sensitivity_MC[_withLegend].png/.pdf and _legend.md."""
import os, sys; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from paths import *
import pandas as pd, numpy as np, json, textwrap, matplotlib; matplotlib.use('Agg')
os.chdir(MC)  # all Monte Carlo files live in <outputs>/mecha/sensitivity_mc/

import matplotlib.pyplot as plt, matplotlib.ticker as mt
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': '#8a8984',
                     'xtick.color': '#52514e', 'ytick.color': '#52514e', 'pdf.fonttype': 42})
MAIN, LAT, TXT, MUTED, GRID = '#2a78d6', '#1baf7a', '#1f1f1d', '#8a8984', '#eeede8'
TC = {'none': '#8a8984', 'B1': '#e34948', 'B10': '#eda100', 'B6': '#7a4bb8'}
TL = {'none': 'No barrier', 'B1': 'B1 endodermis', 'B10': 'B10 exodermis', 'B6': 'B6 both'}
ORDER = ['MP', 'MM', 'MT', 'MT0', 'LP', 'LM', 'LT']; TR = ['none', 'B1', 'B10', 'B6']; P = ['kw', 'kmb', 'kaqp', 'kpl']
PL = {'kw': 'Wall kw', 'kmb': 'Membrane kmb', 'kaqp': 'Aquaporin kAQP', 'kpl': 'Plasmodesmata Kpl'}
SM = pd.read_csv('mc_segment_means.csv').set_index('fragment').reindex(ORDER); E = pd.read_csv('mc_elasticities.csv'); R = pd.read_csv('mc_barrier_ratios.csv')
K = json.load(open('mc_rank_checks.json')); rep = {r: pd.read_csv(f'mc_rep_x{r}.csv') for r in (2, 10)}
SUP = str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹')
def sci(x): m, e = f'{x:.2e}'.split('e'); return f'{float(m):.2f} × 10' + str(int(e)).translate(SUP)

def draw(out, leg=None):
    WR = 175; H = 8.4; LH = 0 if leg is None else sum(0.158 * len(textwrap.wrap(p, WR)) + 0.14 for p in leg.split('\n\n')) + 0.35
    fig = plt.figure(figsize=(13, H + LH)); top = 1 - 0.75 / (H + LH); bot = (LH + 0.55) / (H + LH)
    gs = fig.add_gridspec(2, 2, left=0.07, right=0.985, top=top, bottom=bot, hspace=0.42, wspace=0.22)
    # A segment means
    a = fig.add_subplot(gs[0, 0]); x = np.arange(len(ORDER))
    for i, f in enumerate(ORDER):
        c = MAIN if f.startswith('M') else LAT; r = SM.loc[f]
        a.plot([i, i], [r.mc_q025, r.mc_q975], color=c, lw=6, alpha=0.35, solid_capstyle='round')
        a.scatter(i, r.default_mean, s=46, color=c, edgecolors='white', lw=1, zorder=3)
        a.text(i, r.mc_q975 * 1.12, f'×{r.fold_95:.1f}', ha='center', va='bottom', fontsize=8, color='#52514e')
    a.set_yscale('log'); a.set_xticks(x, ORDER); a.set_xlim(-0.6, len(ORDER) - 0.4); a.set_ylim(SM.mc_q025.min() * 0.7, SM.mc_q975.max() * 2)
    a.set_ylabel('Segment mean kr (cm hPa$^{-1}$ d$^{-1}$)'); a.grid(axis='y', color=GRID, lw=0.6); a.set_axisbelow(True)
    a.set_title('A  Segment kr, default (dot) and 95% Monte Carlo range (bar)', loc='left', fontweight='bold', fontsize=10)
    # B spread by barrier treatment, both ranges
    b = fig.add_subplot(gs[0, 1]); pos = []
    for j, t in enumerate(TR):
        for k, r in enumerate((2, 10)):
            d = rep[r][rep[r].treatment == t]; v = d.kr / d.groupby('anatomy').kr.transform('median'); p = j * 3 + k
            bp = b.boxplot(v, positions=[p], widths=0.7, whis=(2.5, 97.5), showfliers=False, patch_artist=True,
                           medianprops=dict(color=TXT, lw=1.2), whiskerprops=dict(color=TC[t]), capprops=dict(color=TC[t]),
                           boxprops=dict(facecolor=TC[t], alpha=0.35 if r == 10 else 0.8, edgecolor=TC[t]))
            b.text(p, v.quantile(.975) * 1.15, f'×0.5–2' if r == 2 else '×0.1–10', ha='center', va='bottom', fontsize=7.5, color='#52514e')
    b.set_yscale('log'); b.set_xticks([j * 3 + 0.5 for j in range(4)], [TL[t] for t in TR]); b.axhline(1, color=MUTED, lw=0.8, ls=(0, (3, 2)))
    b.yaxis.set_major_formatter(mt.FuncFormatter(lambda v, _: f'{v:g}')); b.yaxis.set_minor_formatter(mt.NullFormatter()); b.set_ylim(0.12, 14); b.set_ylabel('kr relative to the anatomy\'s median'); b.grid(axis='y', color=GRID, lw=0.6); b.set_axisbelow(True)
    b.set_title('B  Spread of kr per barrier type (18 anatomies, 200 draws each)', loc='left', fontweight='bold', fontsize=10)
    # C elasticities
    c = fig.add_subplot(gs[1, 0]); e = E[E.range == 'x2'].set_index('treatment').reindex(TR); w = 0.19
    for i, p in enumerate(P):
        c.bar(np.arange(4) + (i - 1.5) * w, e[f'elasticity_{p}'], width=w * 0.9, color=['#2f4b7c', '#b3b1aa', '#3a9fbf', '#c2185b'][i], label=PL[p])
    c.set_xticks(range(4), [TL[t] for t in TR]); c.set_ylabel('Elasticity (Δlog kr / Δlog parameter)'); c.set_ylim(0, 1.05)
    c.legend(frameon=False, ncol=4, fontsize=8, loc='upper right', bbox_to_anchor=(1, 1.02)); c.grid(axis='y', color=GRID, lw=0.6); c.set_axisbelow(True)
    c.set_title('C  Which cell property drives kr (×0.5–2 set)', loc='left', fontweight='bold', fontsize=10)
    # D B10/B1 vs Kpl/kAQP
    dd = fig.add_subplot(gs[1, 1])
    for r, al, s in ((10, 0.25, 7), (2, 0.6, 9)):
        w_ = rep[r].pivot_table(index=['draw', 'anatomy'], columns='treatment', values='kr').reset_index()
        f = rep[r].drop_duplicates('draw').set_index('draw'); w_['ratio'] = (f.f_kpl / f.f_kaqp).reindex(w_.draw).values
        dd.scatter(w_.ratio, w_.B10 / w_.B1, s=s, alpha=al, color='#2f4b7c' if r == 2 else MUTED, lw=0, label=('×0.5–2' if r == 2 else '×0.1–10'))
    dd.set_xscale('log'); dd.set_yscale('log'); dd.xaxis.set_minor_formatter(mt.NullFormatter()); dd.yaxis.set_minor_formatter(mt.NullFormatter()); dd.axhline(1, color=MUTED, lw=0.8, ls=(0, (3, 2))); dd.axvline(1, color=MUTED, lw=0.8, ls=(0, (3, 2)))
    for ax in (dd,):
        ax.xaxis.set_major_formatter(mt.FuncFormatter(lambda v, _: f'{v:g}')); ax.yaxis.set_major_formatter(mt.FuncFormatter(lambda v, _: f'{v:g}'))
    dd.set_xlabel('Kpl factor / kAQP factor'); dd.set_ylabel('kr(B10) / kr(B1), same anatomy and draw'); dd.legend(frameon=False, fontsize=8, loc='upper left', markerscale=2)
    dd.grid(color=GRID, lw=0.6); dd.set_axisbelow(True)
    dd.set_title('D  B10 vs B1 depends on the Kpl / kAQP balance', loc='left', fontweight='bold', fontsize=10)
    fig.suptitle('Sensitivity of radial conductivity kr to cell hydraulic properties (Monte Carlo)', x=0.012, y=1 - 0.15 / (H + LH), ha='left', va='top', fontsize=12.5, fontweight='bold')
    if leg:
        y = (LH - 0.1) / (H + LH)
        for k, p in enumerate(leg.split('\n\n')):
            L = textwrap.wrap(p, WR); fig.text(0.012, y, '\n'.join(L), fontsize=8.2, va='top', fontweight='bold' if k == 0 else 'normal', color=TXT, linespacing=1.35)
            y -= (0.158 * len(L) + 0.14) / (H + LH)
    fig.savefig(out + '.png', dpi=220); fig.savefig(out + '.pdf'); plt.close(fig)

rr = lambda r, q: R[(R.range == r) & (R.ratio == q)].iloc[0]
el = lambda r, t: E[(E.range == r) & (E.treatment == t)].iloc[0]
segs = '; '.join(f"{f} {sci(r.default_mean)} ({sci(r.mc_q025)}–{sci(r.mc_q975)}, ×{r.fold_95:.1f})" for f, r in SM.iterrows())
LEG = ("Figure. Sensitivity of modelled radial conductivity kr of tomato (S. pimpinellifolium) root segments to cell-scale hydraulic properties (Monte Carlo, suggested by V. Couvreur).\n\n"
 "Four MECHA cell-scale parameters were varied at the same time and independently, each multiplied by a factor drawn log-uniformly between 1/2 and 2 (Latin hypercube, 200 draws): "
 "primary cell-wall conductivity kw (default 2.4 × 10⁻⁴ cm² hPa⁻¹ d⁻¹; septa scaled with it; in no-barrier runs every barrier wall was set to the sampled kw), membrane lipid permeability kmb (3 × 10⁻⁵ cm hPa⁻¹ d⁻¹), "
 "aquaporin contribution kAQP (4.3 × 10⁻⁴ cm hPa⁻¹ d⁻¹, all tissues) and conductance per plasmodesma Kpl (5.3 × 10⁻¹² cm³ hPa⁻¹ d⁻¹; equivalent to scaling plasmodesmata frequency). "
 "Barrier walls stayed at 1 × 10⁻¹⁶. Axial Kx does not depend on these parameters and is not shown. A second set used factors between 1/10 and 10 (B and D).\n\n"
 f"(A) All {int(pd.read_csv('mc_all_x2.csv').anatomy.nunique())} v5 anatomies with their assigned barriers (MP, MM B6; MT B1; LP, LM B10; MT0, LT none), 200 draws (×0.5–2). Dot, mean kr at default parameters; bar, 2.5–97.5 percentile of the segment mean across draws; number, fold range of that interval. "
 f"Default mean (95% range, fold): {segs}. Each draw is applied to all segments at once, so segment rankings can be checked within draws: "
 f"main proximal and middle had the lowest kr in {K['MP,MM < MT,LP,LM (every draw)']*100:.0f}% of draws, and MT0 and LT the highest (MT, LP, LM < MT0 < LT) in {K['MT,LP,LM < MT0 < LT']*100:.0f}%; "
 f"MT was above the lateral proximal/middle mean in {K['MT > mean(LP,LM)']*100:.0f}% of draws. Main/lateral kr ratio (MM/LM): median {K['main/lateral kr ratio (MM/LM) median, 95%'][0]:.2f} (95%: {K['main/lateral kr ratio (MM/LM) median, 95%'][1]:.2f}–{K['main/lateral kr ratio (MM/LM) median, 95%'][2]:.2f}). "
 f"Tukey HSD letters (one-way ANOVA on log10 kr) matched the default pattern ({K['default letters']}) in {K['Tukey letters identical to default']*100:.0f}% of draws; in most other draws MT separated from LP and LM, and in a few MP from MM. "
 f"The negative cost–kr relationship within cross-sections held in every draw: repeated-measures r with building C {K['cost-kr within-section r (building C) median, 95%, frac p<0.05'][0]:.2f} (95%: {K['cost-kr within-section r (building C) median, 95%, frac p<0.05'][1]:.2f} to {K['cost-kr within-section r (building C) median, 95%, frac p<0.05'][2]:.2f}), p < 0.05 in {K['cost-kr within-section r (building C) median, 95%, frac p<0.05'][3]*100:.0f}% of draws; with total N {K['cost-kr within-section r (N) median, 95%, frac p<0.05'][0]:.2f}.\n\n"
 "(B) 18 representative anatomies (the 10th, 50th and 90th percentile of default kr in each of MP, MM, MT, LP, LM and LT), each run under every barrier treatment in each draw; kr divided by that anatomy's median for the treatment. Boxes, quartiles; whiskers, 2.5–97.5 percentiles; dark boxes ×0.5–2, light boxes ×0.1–10. "
 f"Ratios within the same anatomy and draw (median, 95% range, ×0.5–2): B6/B10 {rr('x2','B6/B10')['median']:.2f} ({rr('x2','B6/B10').q025:.2f}–{rr('x2','B6/B10').q975:.2f}), B6/B1 {rr('x2','B6/B1')['median']:.2f} ({rr('x2','B6/B1').q025:.2f}–{rr('x2','B6/B1').q975:.2f}), none/B6 {rr('x2','none/B6')['median']:.1f} ({rr('x2','none/B6').q025:.1f}–{rr('x2','none/B6').q975:.1f}); "
 f"B6 was below B1 and B10, and no barrier above all barriers, in every draw of both sets.\n\n"
 f"(C) Elasticities from a regression of log kr on the four log factors with anatomy as a factor (×0.5–2; R² 0.97–1.00). Without barriers kr follows wall conductivity (elasticity {el('x2','none').elasticity_kw:.2f}); "
 f"with an endodermal barrier (B1) it follows aquaporins ({el('x2','B1').elasticity_kaqp:.2f}); with an exodermal barrier (B10) plasmodesmata ({el('x2','B10').elasticity_kpl:.2f}); with both (B6) aquaporins and plasmodesmata ({el('x2','B6').elasticity_kaqp:.2f}, {el('x2','B6').elasticity_kpl:.2f}). Membrane lipid permeability has little effect (≤ 0.08).\n\n"
 f"(D) kr under B10 relative to B1 for the same anatomy and draw, against the ratio of the Kpl and kAQP factors. B10/B1: median {rr('x2','B10/B1')['median']:.2f} (95%: {rr('x2','B10/B1').q025:.2f}–{rr('x2','B10/B1').q975:.2f}) at ×0.5–2 and {rr('x10','B10/B1')['median']:.2f} ({rr('x10','B10/B1').q025:.2f}–{rr('x10','B10/B1').q975:.2f}) at ×0.1–10; "
 f"B10 exceeded B1 in {rr('x2','B10/B1').frac_gt1*100:.0f}% and {rr('x10','B10/B1').frac_gt1*100:.0f}% of cases; log B10/B1 correlates with log(Kpl/kAQP) (r = {R[(R.range=='x2')&R.ratio.str.contains('pearson')]['median'].iloc[0]:.2f} and {R[(R.range=='x10')&R.ratio.str.contains('pearson')]['median'].iloc[0]:.2f}). "
 "Whether the main tip zone (B1) or the lateral proximal and middle segments (B10) conduct more radially is therefore not resolved by the model without tomato-specific aquaporin and plasmodesmata data.\n\n"
 "Model: MECHA (openalea.mecha, MECHA-dev develop 945b2be with a one-line merge fix), GRANAR v5 anatomies built from Emily Miller's measurements, cell length 200 µm, im_scale 1000. Scripts: sensitivity_mc/mc_kr.py, analyse_mc.py, fig_mc.py.")
open('Fig_kr_sensitivity_MC_legend.md', 'w').write(LEG)
draw('Fig_kr_sensitivity_MC'); draw('Fig_kr_sensitivity_MC_withLegend', LEG); print('done')
