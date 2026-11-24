"""Map of the salt-induced loss of root-system conductance (Krs) across the S. pimpinellifolium panel."""
import krs_paths as KP
KP.to_out()
import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from scipy.stats import spearmanr, mannwhitneyu
A = pd.read_csv('krs_accession_geo.csv').dropna(subset=['Latitude', 'logKrs.SminusC'])
A['retained'] = 100 * np.exp(A['logKrs.SminusC'])          # salt Krs as % of control
A['loss'] = 100 - A.retained
INK = '#52514e'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': INK, 'xtick.color': INK, 'ytick.color': INK})
fig = plt.figure(figsize=(10, 7.2)); gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1], wspace=0.28)

ax = fig.add_subplot(gs[0])
main = A[A.Longitude > -85]; gal = A[A.Longitude <= -85]
m = Basemap(llcrnrlon=-82.5, llcrnrlat=-18.8, urcrnrlon=-68.5, urcrnrlat=2.2, resolution='i', projection='merc', ax=ax)
m.fillcontinents(color='#f1f0ec', lake_color='white'); m.drawmapboundary(fill_color='white', color=INK)
m.drawcoastlines(color='#8a8984', linewidth=.6); m.drawcountries(color='#8a8984', linewidth=.6)
m.drawparallels(range(-18, 3, 4), labels=[1, 0, 0, 0], color='#dddcd7', linewidth=.4, fontsize=8)
m.drawmeridians(range(-82, -67, 4), labels=[0, 0, 0, 1], color='#dddcd7', linewidth=.4, fontsize=8)
vmin, vmax = np.nanpercentile(A.loss, [2, 98])
x, y = m(main.Longitude.values, main.Latitude.values)
order = np.argsort(main.loss.values)
sc = ax.scatter(x[order], y[order], c=main.loss.values[order], cmap='Oranges', vmin=vmin, vmax=vmax, s=34,
                edgecolor='white', linewidth=.6, zorder=5)
mk = main.minor.values[order]
ax.scatter(x[order][mk], y[order][mk], s=78, facecolor='none', edgecolor='#0b0b0b', linewidth=.9, zorder=6)
lx, ly = m(-81.9, -3.5); ax.plot(*m([-82.5, -78.5], [-3.5, -3.5]), ls=':', color=INK, lw=.8)
ax.text(*m(-82.3, -3.2), '3.5°S', color=INK, fontsize=7.5)
for nm, lo, la in [('ECUADOR', -77.2, -0.8), ('PERU', -75.2, -9.5)]: ax.text(*m(lo, la), nm, color='#8a8984', fontsize=9, ha='center', fontweight='bold')
if len(gal): ax.text(.02, .02, f'Galápagos accession(s) not shown: {len(gal)}', transform=ax.transAxes, fontsize=7.5, color=INK)
cb = fig.colorbar(sc, ax=ax, shrink=.55, pad=.02); cb.set_label('Krs lost under 100 mM NaCl (% of control, day 5)')
cb.outline.set_edgecolor(INK)
ax.scatter([], [], s=78, facecolor='none', edgecolor='#0b0b0b', label='minority genetic groups (k-means on PCs)')
ax.legend(loc='lower left', frameon=False, fontsize=8, bbox_to_anchor=(0, .04))
ax.set_title(f'a  Krs lost under salt ({len(A)} accessions)', loc='left', fontweight='bold')

b = fig.add_subplot(gs[1])
ec = A.Latitude > -3.5
for sel, col, lab in [(~ec, '#2a78d6', 'Peru (south of 3.5°S)'), (ec, '#eb6834', 'Ecuador (north of 3.5°S)')]:
    d = A[sel]
    b.scatter(d.Latitude, d.retained, s=22, color=col, alpha=.75, lw=0, label=f'{lab}, n = {sel.sum()}')
    b.scatter(d.Latitude[d.minor], d.retained[d.minor], s=60, facecolor='none', edgecolor='#0b0b0b', lw=.8)
b.axvline(-3.5, ls=':', color=INK, lw=.8)
p = mannwhitneyu(A.retained[ec], A.retained[~ec]).pvalue
b.text(.03, .97, f'median retained: Ecuador {A.retained[ec].median():.0f}%, Peru {A.retained[~ec].median():.0f}%\nMann–Whitney P = {p:.2g}',
       transform=b.transAxes, va='top', color=INK)
b.set_xlabel('Latitude (°)'); b.set_ylabel('Salt Krs as % of control (accession mean, day 5)')
b.spines[['top', 'right']].set_visible(False); b.set_ylim(5, 100); b.legend(frameon=False, loc='lower left', fontsize=8)
b.set_title('b  Retention by latitude', loc='left', fontweight='bold')
fig.savefig('krs_salt_map.png', dpi=200, bbox_inches='tight'); fig.savefig('krs_salt_map.pdf', bbox_inches='tight')
