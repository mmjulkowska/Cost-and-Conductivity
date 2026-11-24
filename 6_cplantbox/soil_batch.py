"""Run every combination: 3 topologies x 3 soils x 3 water scenarios x 2 architecture modes x NREP replicates (10 days each)."""
import cpb_paths as CP   # paths; changes into outputs/cplantbox/
import json, os, sys, itertools, time
from multiprocessing import Pool
OUT='soil_runs'; os.makedirs(OUT,exist_ok=True); NREP=int(os.environ.get('NREP',3)); DAYS=10.
def job(a):
    t,soil,sc,arch,rep=a; fn=f"{OUT}/{t}_{soil.replace(' ','-')}_{sc}_{arch}_r{rep}.json"
    if os.path.exists(fn): return fn
    import soil_sim as S
    t0=time.time(); o,p,age0=S.run(t,soil,sc,arch,rep=rep,days=DAYS)
    json.dump(dict(meta=dict(topology=t,soil=soil,scenario=sc,arch=arch,rep=rep,age0=age0,secs=time.time()-t0),rows=o),open(fn,'w'))
    return fn
if __name__=='__main__':
    jobs=list(itertools.product(['C','B','dT'],['loam','sandy loam','silty clay loam'],['wet','drydown','drytop'],['frozen','growing'],range(NREP)))
    jobs.sort(key=lambda j:(j[1]!='sandy loam' or j[2]!='drytop'))   # slow ones first
    with Pool(2) as pool:
        for i,fn in enumerate(pool.imap_unordered(job,jobs)): print(i+1,len(jobs),fn,flush=True)
