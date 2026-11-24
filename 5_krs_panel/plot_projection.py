import krs_paths as KP
KP.to_out()
import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
O = pd.read_csv('krs_per_plant_day.csv'); O['cond'] = O.condition.map({'control': 'C', 'NaCl_100mM': 'S'})
R = pd.read_csv('krs_projected_fit5.csv'); Hc = pd.read_csv('hindcast_kinetic.csv'); H = pd.read_csv('krs_projected_H2.csv')
H2o = pd.read_csv('krs_heritability.csv'); H2o = H2o[H2o.trait == 'logKrs']
COL = {'C': '#2a78d6', 'S': '#eb6834'}; LAB = {'C': 'Control', 'S': '100 mM NaCl'}; INK = '#52514e'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': INK, 'xtick.color': INK, 'ytick.color': INK,
                     'axes.spines.top': False, 'axes.spines.right': False})
fig, ax = plt.subplots(2, 2, figsize=(9.5, 7.6))

a = ax[0, 0]
a.axvspan(5, 10.3, color='#f1f0ec', lw=0); a.text(7.6, 0.2, 'projected', ha='center', color=INK, transform=a.get_xaxis_transform())
for c in 'CS':
    q = O[O.cond == c].groupby('day').Krs.quantile([.25, .5, .75]).unstack() * 1e4
    a.fill_between(q.index, q[.25], q[.75], color=COL[c], alpha=.18, lw=0); a.plot(q.index, q[.5], '-o', color=COL[c], lw=2, ms=5, label=f'{LAB[c]}, traced')
    for rule, ls in [('accel', '-'), ('const', '--')]:
        x = R[(R.cond == c) & (R.rule == rule) & (R.day >= 6)].groupby('day').Krs.quantile([.25, .5, .75]).unstack() * 1e4
        x = pd.concat([q.loc[[5]], x])
        if rule == 'accel': a.fill_between(x.index, x[.25], x[.75], color=COL[c], alpha=.10, lw=0)
        a.plot(x.index, x[.5], ls, color=COL[c], lw=1.6)
a.plot([], [], '-', color=INK, label='accelerating laterals'); a.plot([], [], '--', color=INK, label='constant lateral rate')
a.set_xlabel('Tracing day'); a.set_ylabel('Krs (×10⁻⁴ cm³ hPa⁻¹ d⁻¹)'); a.set_xticks(range(1, 11))
a.set_title('a  Krs, traced and projected (median, IQR)', loc='left', fontweight='bold'); a.legend(frameon=False, fontsize=8, loc='upper left')

b = ax[0, 1]; h = Hc[(Hc.rule == 'accel') & (Hc.day == 5)]
for c in 'CS':
    x = h[h.cond == c]; b.scatter(x.Krs_obs * 1e4, x.Krs_pred * 1e4, s=9, color=COL[c], alpha=.45, lw=0, label=LAB[c])
lim = [0.5, 30]; b.plot(lim, lim, color='#999', lw=1); b.set_xscale('log'); b.set_yscale('log'); b.set_xlim(lim); b.set_ylim(lim)
from matplotlib.ticker import FixedLocator, ScalarFormatter, NullLocator
for axis in (b.xaxis, b.yaxis):
    axis.set_major_locator(FixedLocator([1, 2, 5, 10, 20])); axis.set_minor_locator(NullLocator()); axis.set_major_formatter(ScalarFormatter())
r = np.corrcoef(np.log(h.Krs_obs), np.log(h.Krs_pred))[0, 1]; bias = np.median(h.Krs_pred / h.Krs_obs)
b.text(.03, .97, f'fit on days ≤3, predict day 5\n{len(h)} plants (Big2, 5, 6)\nr(log) = {r:.2f}, median pred/obs = {bias:.2f}',
       transform=b.transAxes, va='top', color=INK)
b.set_xlabel('Traced Krs, day 5 (×10⁻⁴)'); b.set_ylabel('Predicted Krs, day 5 (×10⁻⁴)')
b.set_title('b  Hindcast check', loc='left', fontweight='bold'); b.legend(frameon=False, loc='lower right')

c_ = ax[1, 0]
for c in 'CS':
    x = R[(R.cond == c) & (R.rule == 'accel')].groupby('day')[['Krs', 'TRL']].median()
    base = O[(O.cond == c) & (O.day == 5)][['Krs', 'TRL']].median(); x = pd.concat([base.to_frame(5).T, x[x.index >= 6]])
    c_.plot(x.index, x.TRL / base.TRL, '--', color=COL[c], lw=1.6); c_.plot(x.index, x.Krs / base.Krs, '-', color=COL[c], lw=2, label=LAB[c])
c_.plot([], [], '-', color=INK, label='Krs'); c_.plot([], [], '--', color=INK, label='Total root length')
c_.set_xlabel('Tracing day'); c_.set_ylabel('Fold change from day 5 (median)'); c_.set_xticks(range(5, 11))
c_.set_title('c  Krs grows slower than the root system', loc='left', fontweight='bold'); c_.legend(frameon=False, fontsize=8)

d = ax[1, 1]
for c in 'CS':
    hh = H[(H.cond == c) & (H.rule == 'accel')].set_index('day').H2_mean
    hh = pd.concat([pd.Series({5: float(H2o[H2o.condition == c].H2_mean.iloc[0])}), hh])
    d.plot(hh.index, hh.values, '-o', color=COL[c], lw=2, ms=5, label=LAB[c])
d.axvspan(5.5, 10.3, color='#f1f0ec', lw=0); d.set_ylim(0, 1); d.set_xticks(range(5, 11))
d.set_xlabel('Tracing day (5 = traced)'); d.set_ylabel('H², accession means (ln Krs)')
d.set_title('d  Heritability of projected Krs', loc='left', fontweight='bold'); d.legend(frameon=False)
fig.tight_layout(); fig.savefig('krs_projection_summary.png', dpi=200); fig.savefig('krs_projection_summary.pdf')
