"""MECHA on the v5 anatomies (run from anywhere: python 2_mecha/batch.py [ids]).
Barriers: MP, MM B6; MT B1; MT0 = MT anatomies with no barrier; LP, LM B10; LT no barrier (all barrier walls = kw).
Writes <outputs>/mecha/mecha_raw.csv and conductivities_per_anatomy.csv."""
import os, sys; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from paths import *
import csv, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from concurrent.futures import ProcessPoolExecutor
from run_one import run
import pandas as pd
SPEC={'MP':('MP',6,False),'MM':('MM',6,False),'MT':('MT',1,False),'MT0':('MT',6,True),'LP':('LP',10,False),'LM':('LM',10,False),'LT':('LT',10,True)}
META={'MP':('Main proximal (3 cm from MM)','B6','lignin','lignin + suberin'),'MM':('Main middle (3 cm from MT)','B6','lignin','lignin + suberin'),
      'MT':('Main tip zone (3 cm from tip)','B1','lignin','none'),'MT0':('Main root true tip (MT geometry)','none','none','none'),
      'LP':('Lateral proximal (3 cm from LM)','B10','none','lignin + suberin'),'LM':('Lateral middle (3 cm from LT)','B10','none','lignin + suberin'),'LT':('Lateral true tip','none','none','none')}
ORDER=['MP','MM','MT','MT0','LP','LM','LT']
seg=lambda n:'MP' if n<=9 else 'MM' if n<=20 else 'MT' if n<=26 else 'LP' if n<=40 else 'LM' if n<=50 else 'LT'
def job(a):
    frag,f=a; src,b,kweq=SPEC[frag]
    try: r=run(os.path.join(ANAT,f+'_root.xml'),[b],1000.0,kw_equal=kweq)[0]; return dict(fragment=frag,anatomy=f,Kx=r['Kx'],kr=r['kr'])
    except Exception as e: return dict(fragment=frag,anatomy=f,Kx=float('nan'),kr=float('nan'))
if __name__=='__main__':
    ids=sys.argv[1:] or sorted(x[:-9] for x in os.listdir(ANAT) if x.endswith('_root.xml'))
    jobs=[(fr,f) for f in ids for fr,(src,_,_) in SPEC.items() if src==seg(int(f.split('_')[0][2:]))]
    with ProcessPoolExecutor(2) as ex: res=list(ex.map(job,jobs))
    new=pd.DataFrame(res)
    if sys.argv[1:] and os.path.exists(MECHA_RAW):
        old=pd.read_csv(MECHA_RAW); new=pd.concat([old[~old.anatomy.isin(ids)],new])
    new.to_csv(MECHA_RAW,index=False)
    t=new.copy()
    for i,col in enumerate(['description','mecha_barrier','endodermis','exodermis']): t[col]=t.fragment.map(lambda f:META[f][i])
    t['o']=t.fragment.map(ORDER.index); t=t.sort_values(['o','anatomy']).rename(columns={'Kx':'Kx_cm4_hPa_d','kr':'kr_cm_hPa_d'})
    t[['fragment','description','anatomy','mecha_barrier','endodermis','exodermis','Kx_cm4_hPa_d','kr_cm_hPa_d']].to_csv(COND,index=False)
    print(len(new),'runs;',new.kr.isna().sum(),'without kr:',sorted(new[new.kr.isna()].anatomy.unique()))
