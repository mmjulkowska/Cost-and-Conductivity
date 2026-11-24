"""Accelerating lateral growth, as measured in the PimpiBig tracings: elongation rate r + s*age (cm/day), capped at vmax
(the fastest median rate measured, 2 cm/day) and at lmax."""
import cpb_paths as CP   # paths; changes into outputs/cplantbox/
import plantbox as pb, numpy as np
VMAX = 2.0
class AcceleratingGrowth(pb.GrowthFunction):
    def __init__(self, s, vmax=VMAX): pb.GrowthFunction.__init__(self); self.s = s; self.vmax = vmax
    def _astar(self, r): return max((self.vmax - r)/self.s, 0.) if self.s > 0 else 1e12
    def getLength(self, t, r, k, o):
        if t <= 0: return 0.
        a = self._astar(r)
        l = r*t + 0.5*self.s*t*t if t <= a else r*a + 0.5*self.s*a*a + self.vmax*(t - a)
        return min(k, l)
    def getAge(self, l, r, k, o):
        if l >= k: return 1e9
        a = self._astar(r); la = r*a + 0.5*self.s*a*a
        if l <= la: return (-r + np.sqrt(r*r + 2*self.s*l))/self.s if self.s > 0 else l/r
        return a + (l - la)/self.vmax
    def copy(self): return AcceleratingGrowth(self.s, self.vmax)
