"""Root water uptake from a drying pot: CPlantBox architectures (v4) + root hydraulics (v5 conductivities, conductance.py)
+ 1D Richards equation for the pot (layers 0.5 cm, van Genuchten-Mualem soils).

Pot: cylinder, radius POT_R cm, depth POT_D cm; roots are confined to it. Seed 0.5 cm below the surface.
Soil water: 1D vertical Richards equation (mixed form, implicit, Picard), no flux at the bottom (pot on a saucer).
Root uptake: the root network (radial 2 pi a l kr, axial Kx/l, gravity included) is joined to the soil of the layer each
root piece sits in through a perirhizal conductance 2 pi l K(h_layer) / ln(r_b / a), where r_b = sqrt(V_layer / (pi L_layer))
is half the mean distance between roots in that layer (Doussan / Couvreur-type coupling; soil K at the bulk layer pressure).
Shoot: potential transpiration TPOT cm3/day on a day/night sinusoid (all topologies the same); when the collar xylem pressure
would fall below PSI_WILT hPa, the collar is held at PSI_WILT and transpiration falls (stress).
Water scenarios:
  'wet'     : well watered, water added at the surface as fast as it is taken up
  'drydown' : starts at container capacity (drained pot), no water added
  'drytop'  : top DRYTOP_D cm start dry (h = -15000 hPa), the rest at container capacity, no water added
Architectures:
  'frozen'  : root system grown to 30 cm total root length inside the pot, then kept fixed
  'growing' : starts from the same 30 cm root system and keeps growing during the run
Units: cm, day, hPa (1 hPa ~ 1 cm water head)."""
import cpb_paths as CP   # paths; changes into outputs/cplantbox/
import numpy as np, json, os, sys, time, plantbox as pb, rsa_model as m, conductance as K
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from scipy.linalg import solve_banded

SOILS = {   # Carsel & Parrish (1988): theta_r, theta_s, alpha (1/cm), n, Ks (cm/d)
    'loam':            (0.078, 0.43, 0.036, 1.56, 24.96),
    'sandy loam':      (0.065, 0.41, 0.075, 1.89, 106.1),
    'silty clay loam': (0.089, 0.43, 0.010, 1.23, 1.68)}
POT_R, POT_D, DZ = float(os.environ.get("POT_R", 1.25)), 20.0, 0.5
TPOT = float(os.environ.get("TPOT", 3.0))          # daily potential transpiration, cm3 per day
PSI_WILT = -15000.
DRYTOP_D = 6.0
H_BOTTOM0 = -10.                                     # pressure head at the pot bottom after draining (container capacity)

class Soil:
    def __init__(self, name):
        self.tr, self.ts, self.a, self.n, self.Ks = SOILS[name]; self.mm = 1 - 1/self.n
    def theta(self, h):
        h = np.minimum(h, 0.); return self.tr + (self.ts - self.tr) * (1 + (self.a*np.abs(h))**self.n)**(-self.mm)
    def Se(self, h): return (self.theta(h) - self.tr)/(self.ts - self.tr)
    def K(self, h):
        se = np.clip(self.Se(h), 1e-12, 1.); return self.Ks*se**0.5*(1 - (1 - se**(1/self.mm))**self.mm)**2
    def C(self, h):
        h = np.minimum(h, -1e-6); ah = self.a*np.abs(h)
        return (self.ts - self.tr)*self.mm*self.n*self.a*ah**(self.n-1)*(1 + ah**self.n)**(-self.mm-1)

class Column:
    def __init__(self, soil, scenario):
        self.s = soil; self.N = int(round(POT_D/DZ)); self.z = -(np.arange(self.N) + 0.5)*DZ      # layer centres (cm, negative down)
        self.area = np.pi*POT_R**2
        h = H_BOTTOM0 - (self.z - self.z[-1])                                                       # hydrostatic, drained pot (pressure falls upwards)
        if scenario == 'drytop':   # dry top layer, log-smooth 1 cm transition to the wet subsoil
            w = 1/(1 + np.exp(-(self.z + DRYTOP_D)/0.25)); h = -np.exp(w*np.log(15000.) + (1 - w)*np.log(-h))
        self.h = h
    def water(self): return np.sum(self.s.theta(self.h))*DZ*self.area                                # cm3
    def step(self, dt, sink, top_in):
        """sink: cm3/day taken out of each layer (positive = uptake); top_in: cm3/day added at the surface."""
        s = self.s; N = self.N; hn = self.h.copy(); th_old = s.theta(hn); src = -sink/(self.area*DZ); qtop = top_in/self.area
        h = hn.copy()
        for it in range(60):
            Kc = s.K(h); Kh = np.sqrt(Kc[:-1]*Kc[1:]); Cc = s.C(h) + 1e-5   # geometric-mean interface K: stable across sharp wet/dry fronts
            # dtheta/dt dz = q_in_from_below - q_out_top ; q_up(i+1/2) = -Kh*((h_i - h_{i+1})/dz + 1)
            ab = np.zeros((3, N)); rhs = np.zeros(N)
            diag = Cc*DZ/dt; rhs[:] = Cc*DZ/dt*h - (s.theta(h) - th_old)*DZ/dt + src*DZ
            # flux through interface i+1/2 into layer i (from below) = -Kh*((h_i - h_{i+1})/dz + 1); into layer i+1 the negative
            for_i = Kh/DZ
            diag[:-1] += for_i; diag[1:] += for_i
            ab[0, 1:] = -for_i; ab[2, :-1] = -for_i
            rhs[:-1] += -Kh; rhs[1:] += Kh
            rhs[0] += qtop
            ab[1] = diag
            hnew = solve_banded((1, 1), ab, rhs)
            hnew = np.maximum(hnew, -1e6)
            err = np.max(np.abs(hnew - h)/np.maximum(10., np.abs(h))); h = hnew
            if err < 1e-5: break
        ok = it < 59
        if ok: self.h = h
        return ok

def plant_segments(p):
    S = K.segments(p); return S

class RootHydraulics:
    def __init__(self, S, col, scenario='v5'):
        kr, kx = K.props(scenario); self.idx = {}
        for s in S:
            for q in (s[0], s[1]): self.idx.setdefault(q, len(self.idx))
        N = len(self.idx); self.N = N
        self.a = np.array([self.idx[s[0]] for s in S]); self.b = np.array([self.idx[s[1]] for s in S])
        self.l = np.array([s[2] for s in S]); self.cl = [s[3] for s in S]
        mid = np.array([(s[4] + s[5])/2 for s in S]); self.zseg = mid[:, 2]
        znode = np.zeros(N); znode[self.a] = np.array([s[4][2] for s in S]); znode[self.b] = np.array([s[5][2] for s in S]); self.znode = znode
        self.layer = np.clip((-(self.zseg + 0.5)/DZ).astype(int), 0, col.N - 1)                   # soil z = plant z + 0.5 (seed at -0.5)
        self.rad = np.array([K.RAD[c] for c in self.cl]); self.Rroot = 2*np.pi*self.rad*self.l*np.array([kr[c] for c in self.cl])
        self.gax = np.array([kx[c] for c in self.cl])/self.l
        Lay = np.bincount(self.layer, weights=self.l, minlength=col.N); Vlay = col.area*DZ
        rb = np.sqrt(Vlay/(np.pi*np.maximum(Lay, 1e-9))); self.rb = np.maximum(rb[self.layer], self.rad*1.5)
        self.collar = self.idx[S[0][0]]; self.cls = np.array(self.cl)
    def solve(self, hsoil_layer, Ksoil_layer, zsoil_layer):
        gp = 2*np.pi*self.l*Ksoil_layer[self.layer]/np.log(self.rb/self.rad)
        G = 1/(1/self.Rroot + 1/np.maximum(gp, 1e-30))                                             # radial path: soil -> root xylem
        Hs = hsoil_layer[self.layer] + (self.zseg + 0.5)*0 + zsoil_layer[self.layer]               # total soil head of each piece
        N = self.N; a, b = self.a, self.b
        diag = np.zeros(N); np.add.at(diag, a, self.gax + G/2); np.add.at(diag, b, self.gax + G/2)
        M = coo_matrix((np.r_[-self.gax, -self.gax, diag], (np.r_[a, b, np.arange(N)], np.r_[b, a, np.arange(N)])), shape=(N, N)).tocsr()
        rhs = np.zeros(N); np.add.at(rhs, a, G/2*Hs); np.add.at(rhs, b, G/2*Hs)
        free = np.ones(N, bool); free[self.collar] = False; Mf = M[free][:, free].tocsc(); Mc = M[free][:, [self.collar]].toarray().ravel()
        zc = self.znode[self.collar] + 0.5
        def run(Hc, rhs_):
            H = np.zeros(N); H[self.collar] = Hc; H[free] = spsolve(Mf, rhs_[free] - Mc*Hc)
            u = G/2*(Hs - H[a]) + G/2*(Hs - H[b]); return u
        uA = run(zc, rhs)                     # collar pressure 0 (total head = elevation), actual soil
        uB = run(zc - 1., np.zeros(N)) - 0    # unit drop at the collar, soil heads 0 ...
        uB = uB - run(zc, np.zeros(N))        # ... minus the elevation-only part -> response to a 1 hPa drop
        return uA, uB, zc

def run(topology, soilname, scenario, arch, rep=0, days=8., lat='v5', record=0.25, verbose=False):
    P = json.load(open('params_v4.json'))[topology]['3D']; P = dict(P, pot=(POT_R, POT_D))
    p = m.build(P, 3000 + rep, stochastic=True); t = 0.
    while sum(o.getLength() for o in p.getOrgans(pb.root)) < 30:
        p.simulate(0.01, False); t += 0.01
    age0 = t; soil = Soil(soilname); col = Column(soil, scenario)
    S = plant_segments(p); RH = RootHydraulics(S, col, lat)
    W0 = col.water(); out = []; tt = 0.; dt = 1e-3; cumT = 0.; cumTpot = 0.; added = 0.; next_rec = 0.; next_grow = 0.1
    while tt < days - 1e-9:
        if arch == 'growing' and tt >= next_grow - 1e-9:
            p.simulate(0.1, False); S = plant_segments(p); RH = RootHydraulics(S, col, lat); next_grow += 0.1
        Tp = TPOT*np.pi*max(0., np.sin(2*np.pi*(tt - 0.25)))                                          # day 06:00-18:00, daily total TPOT
        Ks_l = soil.K(col.h); uA, uB, zc = RH.solve(col.h, Ks_l, col.z)
        QA, QB = uA.sum(), uB.sum()                                                                   # uptake = QA + QB * (drop below 0)
        drop = max(0., (Tp - QA)/QB) if QB > 0 else 0.
        psi_c = -drop
        if psi_c < PSI_WILT: psi_c = PSI_WILT; drop = -PSI_WILT
        u = uA + uB*drop; Ta = u.sum()
        sink = np.bincount(RH.layer, weights=u, minlength=col.N)
        top = Ta if scenario == 'wet' else 0.
        while True:
            if col.step(dt, sink, top): break
            dt /= 2
            if dt < 1e-7: raise RuntimeError('soil solver failed')
        tt += dt; cumT += Ta*dt; cumTpot += Tp*dt; added += top*dt
        dt = min(dt*1.3, 0.005, max(next_rec - tt, 1e-5) if next_rec > tt else 0.005)
        if tt >= next_rec - 1e-9:
            cls_u = {c: float(u[RH.cls == c].sum()) for c in K.CLASSES}
            out.append(dict(topology=topology, soil=soilname, scenario=scenario, arch=arch, rep=rep, t=tt, age=age0 + tt, Tpot=Tp, Tact=Ta,
                            psi_collar=psi_c, cumT=cumT, cumTpot=cumTpot, water=col.water(), water0=W0, added=added,
                            TRL=float(RH.l.sum()), Krs_soil=QB, h_profile=col.h.tolist(), uptake_profile=(sink/(col.area*DZ)).tolist(), **{f'u_{c}': v for c, v in cls_u.items()}))
            next_rec += record
            if verbose: print(f'{tt:5.2f} Tp {Tp:5.2f} Ta {Ta:5.2f} psi_c {psi_c:8.0f} W {col.water():6.1f}')
    return out, p, age0

if __name__ == '__main__':
    t0 = time.time(); o, p, a = run('C', 'loam', 'drydown', 'frozen', days=2, verbose=True); print('seconds', time.time() - t0)
