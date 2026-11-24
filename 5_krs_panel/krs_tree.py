"""Root-system conductance (Krs) of a traced root system, built straight from SmartRoot measurements.

A plant is described by its main-root length MR and its first-order laterals as (insertion position from the
main-root base, length). Under the assumptions of conductance.py (uniform soil water potential, no gravity), Krs
depends only on this topology and on the segment lengths, not on angles or 3D shape, so no CPlantBox growth
simulation is needed: the network is assembled directly from the tracing.

Segment classes, conductivities (v5 MECHA medians), radii (Emily's cross-sections) and the solver are the same
as conductance.py (CPlantBox_RSA, 20 Nov 2026):
  main root: 0-TIP cm from tip = MT0 (no barrier), then MT (B1), MM (B6), MP (B6) in STEP-cm bands
  laterals : 0-TIP cm = LT (no barrier), then LM (B10), LP (B10)   ('lat_endo': LP/LM kr x B6/B10 ratio)
Discretisation DX = 0.1 cm (as in the CPlantBox simulations); every lateral insertion point is a main-root node.
Krs unit: cm3 hPa-1 d-1.
"""
import os, numpy as np, pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

import krs_paths as KP
_C = pd.read_csv(KP.COND).groupby('fragment')[['kr_cm_hPa_d', 'Kx_cm4_hPa_d']].median()
_A = pd.read_csv(KP.ANAT_SEG).set_index('segment')
_BR = pd.read_csv(KP.BRACKET).set_index('fragment')
RAD = {s: np.sqrt(_A.loc[s, 'em_root_mm2'] / np.pi) / 10 for s in _A.index}; RAD['MT0'] = RAD['MT']
CLASSES = ['MP', 'MM', 'MT', 'MT0', 'LP', 'LM', 'LT']


def props(scenario='v5'):
    kr = {s: _C.loc[s, 'kr_cm_hPa_d'] for s in CLASSES}; kx = {s: _C.loc[s, 'Kx_cm4_hPa_d'] for s in CLASSES}
    if scenario == 'lat_endo':
        for s in ('LP', 'LM'): kr[s] *= _BR.loc[s, 'ratio_median']
    return kr, kx


def classify(d, main, TIP=1.5, STEP=3.0):
    if main: return 'MT0' if d < TIP else 'MT' if d < TIP + STEP else 'MM' if d < TIP + 2 * STEP else 'MP'
    return 'LT' if d < TIP else 'LM' if d < TIP + STEP else 'LP'


def build(MR, lats, DX=0.1, TIP=1.5, STEP=3.0):
    """Segments (i, j, length, class). Node 0 = collar (main-root base)."""
    pos = np.unique(np.r_[np.arange(0, MR, DX), MR, [min(max(p, 0.), MR) for p, _ in lats]])
    pos = pos[np.r_[True, np.diff(pos) > 1e-9]]
    seg = []
    for k in range(len(pos) - 1):
        l = pos[k + 1] - pos[k]; d = MR - 0.5 * (pos[k] + pos[k + 1])
        seg.append((k, k + 1, l, classify(d, True, TIP, STEP)))
    nid = len(pos)
    for p, L in lats:
        if L <= 0: continue
        a = int(np.argmin(np.abs(pos - min(max(p, 0.), MR))))
        q = np.r_[np.arange(0, L, DX), L]; q = q[np.r_[True, np.diff(q) > 1e-9]]
        prev = a
        for k in range(len(q) - 1):
            l = q[k + 1] - q[k]; d = L - 0.5 * (q[k] + q[k + 1])
            seg.append((prev, nid, l, classify(d, False, TIP, STEP))); prev = nid; nid += 1
    return seg, nid


def solve(seg, N, scenario='v5'):
    kr, kx = props(scenario)
    r = []; c = []; v = []; diag = np.zeros(N); info = []
    for a, b, l, cl in seg:
        g = kx[cl] / l; R = 2 * np.pi * RAD[cl] * l * kr[cl]
        r += [a, b]; c += [b, a]; v += [-g, -g]; diag[a] += g + R / 2; diag[b] += g + R / 2; info.append((a, b, R, cl, l))
    M = coo_matrix((v + list(diag), (r + list(range(N)), c + list(range(N)))), shape=(N, N)).tocsr()
    free = np.ones(N, bool); free[0] = False
    psi = np.zeros(N); psi[0] = -1.
    if free.any():
        psi[free] = spsolve(M[free][:, free].tocsc(), M[free][:, [0]].toarray().ravel())
    up = {k: 0. for k in CLASSES}; ln = {k: 0. for k in CLASSES}
    for a, b, R, cl, l in info: up[cl] += R / 2 * (-psi[a]) + R / 2 * (-psi[b]); ln[cl] += l
    return sum(up.values()), up, ln


def krs(MR, lats, scenario='v5', TIP=1.5, STEP=3.0, DX=0.1):
    seg, N = build(MR, lats, DX, TIP, STEP)
    return solve(seg, N, scenario)
