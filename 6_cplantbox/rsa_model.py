"""CPlantBox model of a tomato (S. pimpinellifolium) seedling root system: main root + 1st-order laterals."""
import cpb_paths as CP   # paths; changes into outputs/cplantbox/
import plantbox as pb, numpy as np
_KEEP = []   # Python growth-function objects must outlive the plants that use them
import os as _os
BASE=_os.path.join(_os.path.dirname(pb.__file__),'modelparameter','structural','rootsystem','Anagallis_femina_Leitner_2010.xml')   # template; every parameter used is overwritten

def build(P, seed=1, stochastic=False, plate=False):
    """P: dict with MR: r, la, lb, ln, a ; LR: r, lmax, theta(deg), tropismN, tropismS, a"""
    p = pb.Plant(); p.readParameters(BASE, verbose=False) if 'verbose' in pb.Plant.readParameters.__doc__ else p.readParameters(BASE)
    rp = p.getOrganRandomParameter(pb.root)
    mr, lr = rp[1], rp[2]
    s = 1.0 if stochastic else 0.0
    mr.r, mr.la, mr.lb, mr.ln, mr.lmax, mr.a = P['r'], P['la'], P['lb'], P['ln'], P.get('mr_lmax',100.), 0.025
    mr.lnk = P.get('lnk', 0.)
    mr.gf = P.get('gfMR', 1)   # 1: slows towards mr_lmax; 2: constant rate
    mr.rs, mr.las, mr.lbs, mr.lns, mr.lmaxs = s*P['r']*P.get('cv_r', .1), 0, 0, s*P['ln']*0.25, 0
    mr.tropismT, mr.tropismN, mr.tropismS = 1, 2.0, 0.15 if stochastic else 0.05
    mr.theta = 0.; mr.dx = 0.1
    mr.successorST = [[2]]; mr.successorP = [[1.]]; mr.successorOT=[[pb.root]]
    lr.r, lr.lmax, lr.a = P['rL'], P['lmaxL'], 0.015
    lr.gf = P.get('gfL', 1)
    lr.rs, lr.lmaxs = s*P['rL']*0.3, s*P['lmaxL']*0.2
    lr.la, lr.lb, lr.ln, lr.lmaxs = 0.1, 0.1, 0., lr.lmaxs
    lr.theta, lr.thetas = np.radians(P['theta']), s*np.radians(10)
    lr.tropismT, lr.tropismN, lr.tropismS = 1, P['tN'], P.get('tS', .2)
    lr.ldelay = P.get('ldelay', -1); lr.dx = 0.1
    lr.successorST=[]; lr.successorP=[]; lr.successorOT=[]
    sp = p.getOrganRandomParameter(pb.seed)[0]
    sp.seedPos = pb.Vector3d(0, 0, -0.5); sp.maxB = 0; sp.firstB = 1e9; sp.delaySB=1e9; sp.firstSB=1e9
    if P.get('pot'):   # cylindrical pot: (radius cm, depth cm); seed sits 0.5 cm below the surface
        r_, d_ = P['pot']; p.setGeometry(pb.SDF_PlantContainer(r_, r_, d_, False))
    elif plate:   # thin slab: roots grow on the agar surface (x-z plane)
        p.setGeometry(pb.SDF_PlantBox(40., 0.05, 60.))
    p.setSeed(seed)
    p.initialize(False)
    if P.get('sL', 0) > 0:   # accelerating lateral elongation (growth.py); set after initialize, keep a reference
        import growth
        g = growth.AcceleratingGrowth(P['sL']); _KEEP.append(g); del _KEEP[:-5000]
        p.getOrganRandomParameter(pb.root)[2].f_gf = g
    if P.get('extend_mr', False):   # constant-rate main root that keeps branching: continue the realised spacing pattern beyond mr_lmax
        p.simulate(1e-4, False)
        prm = [o for o in p.getOrgans(pb.root) if o.getParameter('subType') == 1][0].getParam()
        v = list(prm.ln); d = max(v[-1] - v[-2], 0.) if len(v) > 1 else 0.
        while sum(v) < 200: v.append(max(v[-1] + d, 0.05))
        prm.ln = v
    return p

def measure(p):
    roots = p.getOrgans(pb.root)
    mr = [o for o in roots if o.getParameter('subType') == 1][0]
    out = dict(MR=mr.getLength(), lats=[])
    nodes = np.array([[n.x, n.y, n.z] for n in mr.getNodes()])
    cum = np.r_[0, np.cumsum(np.linalg.norm(np.diff(nodes, axis=0), axis=1))]
    for o in roots:
        if o.getParameter('subType') != 2 or o.getLength() <= 0: continue
        q = np.array([[n.x, n.y, n.z] for n in o.getNodes()])
        i = np.argmin(np.linalg.norm(nodes - q[0], axis=1))
        tip = q[-1] - q[0]; horiz = np.hypot(tip[0], tip[1])
        lc = np.r_[0, np.cumsum(np.linalg.norm(np.diff(q, axis=0), axis=1))]
        k = min(np.searchsorted(lc, min(0.3, lc[-1])), len(q)-1); ld = q[k] - q[0]
        md = nodes[min(i+2, len(nodes)-1)] - nodes[max(i-2, 0)]
        ins = np.degrees(np.arccos(np.clip(ld @ md / (np.linalg.norm(ld)*np.linalg.norm(md) + 1e-12), -1, 1)))
        out['lats'].append((cum[i], o.getLength(), np.degrees(np.arctan2(horiz, -tip[2])), ins))
    out['lats'] = np.array(out['lats']).reshape(-1, 4)
    return out

def simulate(P, T, seed=1, stochastic=False, plate=False, dt=0.1):
    p = build(P, seed, stochastic, plate)
    for _ in range(int(round(T/dt))): p.simulate(dt, False)
    return p
