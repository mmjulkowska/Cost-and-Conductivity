import sys, time, io, contextlib
from openalea.mecha import Mecha, InData
def run(f, barriers, im_scale=1000.0, kw_equal=False):
    buf=io.StringIO()
    with contextlib.redirect_stdout(buf):
        d=InData(cellset_file=f)
        d.geometry.im_scale=im_scale
        d.geometry.set_maturity_stages(barriers, [200.0]*len(barriers))
        if kw_equal:
            kw=d.hydraulic.kw[0]
            for e in d.hydraulic.kw_barrier_elems:
                for k in ('value','Casp','Sub','Sub_in','Sub_out','Lig'): e[k]=kw
            d.hydraulic.kw_phi_thick=[kw]
        m=Mecha(d); m.compute_conductivities()
    return m.root_hydraulic_properties
if __name__=='__main__':
    t=time.time()
    print(run(sys.argv[1],[int(x) for x in sys.argv[2].split(',')],float(sys.argv[3])), time.time()-t)
