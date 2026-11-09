"""Monte Carlo sensitivity of radial conductivity kr to cell-scale hydraulic properties (suggested by V. Couvreur, Nov 2026).
Run from anywhere in the mecha_dev env (outputs go to <outputs>/mecha/sensitivity_mc/):
    python 2_mecha/sensitivity_mc/mc_kr.py all 2        # all v5 anatomies with their assigned barriers, factors in [1/2, 2]      (design 1)
    python 2_mecha/sensitivity_mc/mc_kr.py rep 2        # 18 representative anatomies x every barrier treatment, factors in [1/2, 2]  (design 2)
    python 2_mecha/sensitivity_mc/mc_kr.py rep 10       # same, factors in [1/10, 10]                                               (design 2, wide)
Options: --draws 200 --workers N (default: CPU count - 1). Resumable: finished draws are skipped.
Varied independently, log-uniform, Latin hypercube (seed fixed): kw (primary wall; septa scaled with it), kmb (membrane lipid),
kAQP (aquaporins, all tissues), Kpl (per plasmodesma; equivalent to scaling plasmodesmata frequency).
Barrier walls stay at 1e-16. 'none' = mature-xylem stage with every barrier wall set to the (sampled) kw (Couvreur recipe).
Output: mc_<design>_x<range>.csv (one row per draw x anatomy x treatment) and mc_draws_x<range>.csv (the factors)."""
import os, sys; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from paths import *
import sys, os, io, contextlib, argparse, time, numpy as np, pandas as pd
from multiprocessing import Pool
HERE = MC
PARAMS = ['kw', 'kmb', 'kaqp', 'kpl']
ASSIGNED = {'MP': ('B6', 6), 'MM': ('B6', 6), 'MT': ('B1', 1), 'LP': ('B10', 10), 'LM': ('B10', 10), 'LT': ('none', 0)}
TREAT = {'none': 0, 'B1': 1, 'B10': 10, 'B6': 6}
seg = lambda n: 'MP' if n <= 9 else 'MM' if n <= 20 else 'MT' if n <= 26 else 'LP' if n <= 40 else 'LM' if n <= 50 else 'LT'

def draws(n, rng_fold, seed=20261107):
    rng = np.random.default_rng(seed + int(rng_fold))
    u = np.empty((n, len(PARAMS)))
    for j in range(len(PARAMS)): u[:, j] = (rng.permutation(n) + rng.random(n)) / n   # Latin hypercube
    f = np.exp((2 * u - 1) * np.log(rng_fold))                                           # log-uniform in [1/r, r]
    D = pd.DataFrame(f, columns=PARAMS); D.insert(0, 'draw', range(n)); return D

def run(f, b, fac):
    from openalea.mecha import Mecha, InData
    none = (b == 0); bb = 10 if none else b
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        d = InData(cellset_file=f); d.geometry.im_scale = 1000.0; d.geometry.set_maturity_stages([bb], [200.0]); h = d.hydraulic
        h.kw = [x * fac['kw'] for x in h.kw]; h.kw_elems[0]['value'] *= fac['kw']
        h.kw_septa = [x * fac['kw'] for x in h.kw_septa]; h.kw_septa_elems[0]['value'] *= fac['kw']
        h.kmb *= fac['kmb']; h.kaqp[0]['value'] *= fac['kaqp']; h.kpl[0]['value'] *= fac['kpl']
        if none:
            kw = h.kw[0]
            for e in h.kw_barrier_elems:
                for k in ('value', 'Casp', 'Sub', 'Sub_in', 'Sub_out', 'Lig'): e[k] = kw
            h.kw_phi_thick = [kw]
        m = Mecha(d); m.compute_conductivities()
    return m.root_hydraulic_properties[-1]['kr']

def job(args):
    draw, fac, tasks = args; out = []
    for a, frag, tr, b in tasks:
        try: kr = run(os.path.join(ANAT, a + '_root.xml'), b, fac)
        except Exception as e: kr = np.nan
        out.append(dict(draw=draw, anatomy=a, fragment=frag, treatment=tr, kr=kr, **{'f_' + p: fac[p] for p in PARAMS}))
    return out

def tasks_for(design):
    C = pd.read_csv(COND); C = C[C.kr_cm_hPa_d.notna()]
    if design == 'all':
        return [(r.anatomy, r.fragment, 'none' if r.mecha_barrier.startswith('none') else r.mecha_barrier, 0 if r.mecha_barrier.startswith('none') else int(r.mecha_barrier[1:]))
                for r in C.itertuples()]
    reps = representative()
    return [(a, s, tr, b) for s, ids in reps.items() for a in ids for tr, b in TREAT.items()]

def representative():
    """3 anatomies per segment at the 10th, 50th and 90th percentile of v5 kr (assigned barrier)."""
    C = pd.read_csv(COND)
    C = C[(C.fragment != 'MT0') & C.kr_cm_hPa_d.notna()]; out = {}
    for s, g in C.groupby('fragment'):
        g = g.sort_values('kr_cm_hPa_d').reset_index(drop=True); n = len(g)
        out[s] = [g.anatomy[int(round(q * (n - 1)))] for q in (0.1, 0.5, 0.9)]
    return out

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('design', choices=['all', 'rep']); ap.add_argument('range', type=float)
    ap.add_argument('--draws', type=int, default=200); ap.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument('--first', type=int, default=None, help='only the first N draws (pilot)')
    a = ap.parse_args(); r = int(a.range)
    D = draws(a.draws, a.range); D.to_csv(os.path.join(HERE, f'mc_draws_x{r}.csv'), index=False)
    T = tasks_for(a.design); out = os.path.join(HERE, f'mc_{a.design}_x{r}.csv')
    done = set(pd.read_csv(out).draw) if os.path.exists(out) else set()
    todo = [(int(d.draw), {p: float(getattr(d, p)) for p in PARAMS}, T) for d in D.itertuples() if d.draw not in done]
    if a.first is not None: todo = [t for t in todo if t[0] < a.first]
    print(f'{a.design} x{r}: {len(T)} runs per draw, {len(todo)} draws to do, {a.workers} workers', flush=True)
    t0 = time.time()
    with Pool(a.workers) as pool:
        for i, res in enumerate(pool.imap_unordered(job, todo), 1):
            pd.DataFrame(res).to_csv(out, mode='a', header=not os.path.exists(out), index=False)
            el = time.time() - t0; print(f'  draw {res[0]["draw"]} done ({i}/{len(todo)}, {el/60:.1f} min, ~{el/i*(len(todo)-i)/60:.0f} min left)', flush=True)
    print('finished', out)
