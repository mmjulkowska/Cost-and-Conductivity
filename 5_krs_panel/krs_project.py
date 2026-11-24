"""Project each traced plant's root system past its last tracing day from its own growth, and compute Krs.

Per plant (fitted on the days up to FIT_LAST, default 5; the hindcast uses 3):
  MR rate r_i       least-squares slope of main-root length over the plant's traced days, >= 0, held constant (as v4)
  branching zone    la_i = MR - most distal lateral insertion on the last fitted day (type x treatment median if no laterals)
  spacing           ln_i = median gap between successive lateral insertions on the last day (median of the group if < 2 laterals)
  new laterals      appear at last insertion + ln_i once the main root has grown la_i beyond that point; they enter
                    (ln_i = 1 / slope of lateral number on main-root length over the traced days; 5 cm, the cap, when the
                    plant made no new laterals; group median when it was traced on one day only)
                    the network once >= 0.1 cm long (the tracing threshold)
Lateral elongation, two rules:
  A 'accel'  rate(age) = k_i * (rL + sL*age), capped at 2 cm/d. rL, sL fitted per architectural type x treatment on
             laterals of known age (Big2/5/6, followed day to day by insertion position). An existing lateral's age is
             read off its length on that curve. k_i = plant's observed day-to-day lateral growth / curve prediction.
  B 'const'  every lateral of plant i grows at the plant's own mean observed daily lateral elongation, constant.
Laterals are followed between consecutive days by insertion position (within 0.15 cm, not shorter by > 0.05 cm),
the rule used in track_laterals.py (99.8% agreement with SmartRoot IDs).
No second-order laterals, no maximal lateral length, no plate edge.
"""
import numpy as np, pandas as pd, krs_tree as T, sys, os
import krs_paths as KP
KP.to_out()
FIT_LAST = int(os.environ.get('FIT_LAST', 5)); HORIZON = int(os.environ.get('HORIZON', 10))
DT = 0.05; VMAX = 2.0; MINLEN = 0.1
P = pd.read_csv('plants_assembled.csv'); L = pd.read_csv('laterals_assembled.csv')
P['type'] = P.architectural_type.fillna('NA'); P['cond'] = P.condition.map({'control': 'C', 'NaCl_100mM': 'S'})
k3 = ['experiment', 'label']


def match(a, b, tol=0.15):
    """a, b: arrays (pos, len). Returns list of (i, j)."""
    if len(a) == 0 or len(b) == 0: return []
    d = np.abs(a[:, None, 0] - b[None, :, 0]); ok = (d < tol) & (b[None, :, 1] >= a[:, None, 1] - 0.05)
    c = sorted((d[i, j], i, j) for i, j in zip(*np.nonzero(ok))); ua = set(); ub = set(); out = []
    for _, i, j in c:
        if i in ua or j in ub: continue
        ua.add(i); ub.add(j); out.append((i, j))
    return out


LG = {k: v[['insertion_position', 'length']].to_numpy() for k, v in L.groupby(['experiment', 'label', 'day'])}

# ---- lateral increments between consecutive days (all plants), with age where known -------------------------
inc = []
for (ex, lab), g in P.groupby(k3):
    days = sorted(g.day); age = {}
    for d0, d1 in zip(days[:-1], days[1:]):
        if d1 - d0 != 1: age = {}; continue
        a = LG.get((ex, lab, d0), np.zeros((0, 2))); b = LG.get((ex, lab, d1), np.zeros((0, 2)))
        pr = match(a, b); mb = {j for _, j in pr}; new = {}
        for j in range(len(b)):
            if j not in mb: new[j] = 0 if d0 >= min(days) else None
        for i, j in pr:
            ag = age.get(i); new[j] = None if ag is None else ag + 1
            inc.append(dict(experiment=ex, label=lab, day=d1, age=None if ag is None else ag + 1, L0=a[i, 1], L1=b[j, 1], rate=b[j, 1] - a[i, 1]))
        age = new
INC = pd.DataFrame(inc).merge(P.drop_duplicates(k3)[k3 + ['type', 'cond']], on=k3)
INC['age'] = INC.age.astype(float)


def curve(sub):
    s = sub.dropna(subset=['age']); s = s[s.age <= 3]
    if s.age.nunique() < 2 or len(s) < 20: return None
    m = s.groupby('age').rate.agg(['median', 'size'])
    sL, rL = np.polyfit(m.index.values - 0.5, m['median'].values, 1, w=np.sqrt(m['size'].values))   # rate over the day = curve at mid-age
    return max(rL, 0.05), max(sL, 0.0)


def build_curves(inc):
    C = {}
    for c in 'CS':
        base = curve(inc[inc.cond == c])
        for t in ['B', 'C', 'T', 'dT', 'NA']:
            C[(t, c)] = curve(inc[(inc.cond == c) & (inc.type == t)]) or base
    return C


def L_of_age(a, rL, sL, k=1.):
    r0 = k * rL; s = k * sL; ast = (VMAX - r0) / s if s > 0 and r0 < VMAX else (0. if r0 >= VMAX else 1e9)
    a = np.asarray(a, float); l1 = r0 * np.minimum(a, ast) + s * np.minimum(a, ast) ** 2 / 2
    return l1 + VMAX * np.maximum(a - ast, 0)


def age_of_L(l, rL, sL, k=1.):
    grid = np.linspace(0, 30, 6001); return np.interp(l, L_of_age(grid, rL, sL, k), grid)


def plant_params(g, inc_p, C, fit_last):
    g = g[g.day <= fit_last].sort_values('day'); last = g.iloc[-1]
    r = max(np.polyfit(g.day, g.MR, 1)[0], 0.) if len(g) >= 2 else np.nan
    lat = LG.get((last.experiment, last.label, last.day), np.zeros((0, 2)))
    pos = np.sort(lat[:, 0]) if len(lat) else np.array([])
    la = last.MR - pos.max() if len(pos) else np.nan
    LNMODE = os.environ.get('LNMODE', 'kinetic')
    ln = (np.median(np.diff(pos)) if LNMODE == 'median' else (pos.max() - pos.min()) / (len(pos) - 1)) if len(pos) >= 3 else np.nan
    if LNMODE == 'kinetic':      # spacing from the plant's own emergence kinetics: cm of main-root growth per new lateral
        sl = np.polyfit(g.MR, g.nLR, 1)[0] if len(g) >= 2 and g.MR.nunique() >= 2 else np.nan
        if np.isnan(sl): ln = np.nan                       # one traced day only -> group median
        elif sl <= 1e-6: ln = 5.                           # no new laterals over the traced days -> sparsest spacing allowed
        else: ln = 1 / sl
    rL, sL = C[(last.type, last.cond)]
    ip = inc_p[inc_p.day <= fit_last]
    if len(ip):
        a0 = age_of_L(ip.L0.values, rL, sL); pred = L_of_age(a0 + 1, rL, sL) - L_of_age(a0, rL, sL)
        k = float(np.clip(ip.rate.sum() / pred.sum(), 0.1, 5)) if pred.sum() > 0 else np.nan
        vc = max(ip.rate.mean(), 0.)
    else: k = np.nan; vc = np.nan
    return dict(day0=int(last.day), MR0=last.MR, r=r, la=la, ln=ln, rL=rL, sL=sL, k=k, vconst=vc, n_inc=len(ip), lat0=lat)


def fill(Pp):
    """group medians (type x treatment) for parameters a plant could not estimate"""
    for c in ['r', 'la', 'ln', 'k', 'vconst']:
        med = Pp.groupby(['type', 'cond'])[c].transform('median'); Pp[c + '_imputed'] = Pp[c].isna(); Pp[c] = Pp[c].fillna(med)
        Pp[c] = Pp[c].fillna(Pp.groupby('cond')[c].transform('median'))
    Pp['ln'] = Pp.ln.clip(lower=0.1, upper=5.); Pp['la'] = Pp.la.clip(lower=0.)
    return Pp


def grow(p, rule, days_out):
    """simulate from day0 to each day in days_out; returns {day: (MR, [(pos, len), ...])}"""
    MR = p['MR0']; lat = [[x, l] for x, l in p['lat0']]
    if rule == 'accel': ages = [float(age_of_L(l, p['rL'], p['sL'], p['k'])) for _, l in lat]
    nxt = max(x for x, _ in lat) + p['ln'] if lat else max(MR - p['la'], 0.)
    out = {}; t = p['day0']; steps = int(round((max(days_out) - t) / DT))
    for s in range(1, steps + 1):
        MR += p['r'] * DT
        while MR - p['la'] >= nxt:              # new lateral emerges
            lat.append([nxt, 0.]); nxt += p['ln']
            if rule == 'accel': ages.append(0.)
        if rule == 'accel':
            for q in range(len(lat)):
                v = min(p['k'] * (p['rL'] + p['sL'] * ages[q]), VMAX); lat[q][1] += v * DT; ages[q] += DT
        else:
            for q in range(len(lat)): lat[q][1] += p['vconst'] * DT
        tt = p['day0'] + s * DT
        for d in days_out:
            if abs(tt - d) < DT / 2: out[d] = (MR, [(x, l) for x, l in lat if l >= MINLEN])
    return out


def run(fit_last, horizon, plants=None, curves=None):
    inc = INC[INC.day <= fit_last]; C = curves or build_curves(inc)
    rows = []; prm = []
    groups = P.groupby(k3)
    for (ex, lab), g in groups:
        if plants is not None and (ex, lab) not in plants: continue
        if (g.day <= fit_last).sum() == 0: continue
        p = plant_params(g, inc[(inc.experiment == ex) & (inc.label == lab)], C, fit_last)
        p.update(experiment=ex, label=lab, type=g.type.iloc[0], cond=g.cond.iloc[0]); prm.append(p)
    Pp = fill(pd.DataFrame(prm))
    for p in Pp.to_dict('records'):
        days_out = list(range(p['day0'] + 1, horizon + 1))
        if not days_out: continue
        for rule in ('accel', 'const'):
            for d, (MR, lats) in grow(p, rule, days_out).items():
                k, up, _ = T.krs(MR, lats)
                rows.append(dict(experiment=p['experiment'], label=p['label'], cond=p['cond'], type=p['type'], rule=rule,
                                 day=d, MR=MR, nLR=len(lats), TRL=MR + sum(l for _, l in lats), Krs=k, up_LT=up['LT'] / k))
    return pd.DataFrame(rows), Pp.drop(columns=['lat0']), C


if __name__ == '__main__':
    R, Pp, CUR = run(FIT_LAST, HORIZON)
    tag = f'fit{FIT_LAST}'
    R.to_csv(f'krs_projected_{tag}.csv', index=False); Pp.to_csv(f'growth_params_{tag}.csv', index=False)
    pd.DataFrame([(t, c, *v) for (t, c), v in CUR.items()], columns=['type', 'cond', 'rL', 'sL']).to_csv(f'lateral_curves_{tag}.csv', index=False)
    print(pd.DataFrame([(t, c, *v) for (t, c), v in CUR.items()], columns=['type', 'cond', 'rL', 'sL']).round(3))
    print(Pp[[c for c in Pp if c.endswith('_imputed')]].mean().round(3))
    print(Pp.groupby('cond')[['r', 'la', 'ln', 'k', 'vconst']].median().round(3))
    print((R.groupby(['rule', 'cond', 'day']).Krs.median().unstack() * 1e4).round(1))
    # heritability of projected ln Krs per day (point estimates; same model as krs_heritability.py)
    import statsmodels.formula.api as smf, warnings; warnings.filterwarnings('ignore')
    R['acc'] = R.label.str.split('_').str[0]; rows = []
    for (rule, c, d), g in R[R.day >= 6].groupby(['rule', 'cond', 'day']):
        g = g.assign(y=np.log(g.Krs)); m = smf.mixedlm('y ~ C(experiment)', g, groups=g.acc).fit(reml=True)
        vg = float(m.cov_re.iloc[0, 0]); ve = m.scale; n = g.groupby('acc').size(); nh = len(n) / np.sum(1 / n)
        rows.append(dict(rule=rule, cond=c, day=d, H2_plant=vg / (vg + ve), H2_mean=vg / (vg + ve / nh)))
    H = pd.DataFrame(rows); H.to_csv('krs_projected_H2.csv', index=False); print(H.round(2).to_string())
