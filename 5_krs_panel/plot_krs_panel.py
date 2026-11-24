import krs_paths as KP
KP.to_out()
import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
R = pd.read_csv('krs_per_plant_day.csv', dtype={'accession_code': str})
H = pd.read_csv('krs_heritability.csv')
COL = {'control': '#2a78d6', 'NaCl_100mM': '#eb6834'}; LAB = {'control': 'Control', 'NaCl_100mM': '100 mM NaCl'}
INK = '#52514e'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': INK, 'axes.labelcolor': '#0b0b0b', 'xtick.color': INK, 'ytick.color': INK,
                     'axes.spines.top': False, 'axes.spines.right': False, 'font.family': 'DejaVu Sans'})
fig, ax = plt.subplots(2, 2, figsize=(9.5, 7.5))

# a: Krs over tracing days
a = ax[0, 0]
for c in ['control', 'NaCl_100mM']:
    g = R[R.condition == c].groupby('day').Krs; q = g.quantile([.25, .5, .75]).unstack() * 1e4
    a.fill_between(q.index, q[.25], q[.75], color=COL[c], alpha=.18, lw=0)
    a.plot(q.index, q[.5], '-o', color=COL[c], lw=2, ms=6, label=LAB[c])
n = R.groupby('day').size()
a.set_xticks(range(1, 6)); a.set_xticklabels([f'{d}\n(n={n[d]})' for d in range(1, 6)])
a.set_xlabel('Tracing day'); a.set_ylabel('Krs (×10⁻⁴ cm³ hPa⁻¹ d⁻¹)')
a.set_title('a  Root-system conductance over time', loc='left', fontweight='bold'); a.legend(frameon=False)

# b: Krs vs TRL at day 5
b = ax[0, 1]; d5 = R[R.day == 5]
for c in ['control', 'NaCl_100mM']:
    x = d5[d5.condition == c]; b.scatter(x.TRL, x.Krs * 1e4, s=10, color=COL[c], alpha=.45, lw=0, label=LAB[c])
b.set_xscale('log'); b.set_yscale('log'); from matplotlib.ticker import FixedLocator, NullLocator, ScalarFormatter
b.xaxis.set_major_locator(FixedLocator([5,10,20,40,80])); b.xaxis.set_minor_locator(NullLocator()); b.xaxis.set_major_formatter(ScalarFormatter())
b.yaxis.set_major_locator(FixedLocator([1,2,5,10,20])); b.yaxis.set_minor_locator(NullLocator()); b.yaxis.set_major_formatter(ScalarFormatter()); b.set_xlabel('Total root length, day 5 (cm)'); b.set_ylabel('Krs (×10⁻⁴)')
r2 = np.corrcoef(np.log(d5.TRL), np.log(d5.Krs))[0, 1] ** 2
b.text(.03, .95, f'ln Krs ~ ln TRL: R² = {r2:.2f}', transform=b.transAxes, va='top', color=INK)
b.set_title('b  Krs is mostly root size', loc='left', fontweight='bold'); b.legend(frameon=False, loc='lower right')

# c: salt/control ratio of medians per day, overall and Krs per cm
c_ = ax[1, 0]
for v, ls, lab in [('Krs', '-', 'Krs'), ('TRL', '--', 'Total root length'), ('Krs_per_TRL', ':', 'Krs per cm root')]:
    m = R.groupby(['day', 'condition'])[v].median().unstack(); r = m['NaCl_100mM'] / m['control']
    c_.plot(r.index, r.values, ls, marker='o', color='#0b0b0b' if v == 'Krs' else INK, lw=2 if v == 'Krs' else 1.5, ms=5, label=lab)
c_.axhline(1, color='#bbb', lw=1); c_.set_xticks(range(1, 6)); c_.set_ylim(0, 1.2)
c_.set_xlabel('Tracing day'); c_.set_ylabel('Salt / control (ratio of medians)')
c_.set_title('c  Salt effect grows with exposure', loc='left', fontweight='bold'); c_.legend(frameon=False, loc='lower left')

# d: heritability
d = ax[1, 1]
NAMES = {'logKrs': 'Krs', 'logTRL': 'Total root length', 'KrsPerTRL': 'Krs per cm', 'KrsResid': 'Krs beyond size', 'KrsRGR': 'Krs gain d4→d5'}
tr = list(NAMES); y = np.arange(len(tr))
for k, c in enumerate(['C', 'S']):
    h = H[H.condition == c].set_index('trait').loc[tr]; off = -0.18 + 0.36 * k
    col = COL['control' if c == 'C' else 'NaCl_100mM']
    d.barh(y + off, h.H2_mean, height=.32, color=col, label=LAB['control' if c == 'C' else 'NaCl_100mM'])
    d.errorbar(h.H2_mean, y + off, xerr=[h.H2_mean - h.H2_mean_lo, h.H2_mean_hi - h.H2_mean], fmt='none', ecolor=INK, lw=1, capsize=2)
d.set_yticks(y); d.set_yticklabels([NAMES[t] for t in tr]); d.invert_yaxis(); d.set_xlim(0, 1)
d.set_xlabel('Broad-sense H², accession means (95% bootstrap CI)')
d.set_title('d  Heritability (day 5)', loc='left', fontweight='bold'); d.legend(frameon=False, loc='upper center', bbox_to_anchor=(0.5,-0.12), ncol=2)
fig.tight_layout(); fig.savefig('krs_panel_summary.png', dpi=200); fig.savefig('krs_panel_summary.pdf')
