"""Read SmartRoot RSML of PimpiBig2, days 1-5, into one table per root per day (input preparation for step 5).
Source: Dropbox Tomato_Pimpinellifolium/Pimp_RSA/SmartRooted/pimpiBIG2 (set PIMPIBIG2_RSML to that folder).
Writes timeseries_roots.csv into KRS_IN: day, file, plant label, root id, parent id, order, length (cm),
insertion distance from MR base (cm, laterals), base and tip coordinates (cm). 300 dpi -> cm."""
import xml.etree.ElementTree as ET, numpy as np, pandas as pd, glob, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import krs_paths as KP
B=os.path.join(os.path.expanduser(os.environ.get('PIMPIBIG2_RSML', '~/mnt/Pimp_RSA/SmartRooted/pimpiBIG2')), '')
DAYS={1:B+'done/d1/_set1_day1_*.rsml',2:B+'done/d2/_set1_day2_*.rsml',3:B+'done/d3/_set1_day3_*.rsml',4:B+'_set1_day4_*.rsml',5:B+'_set1_day5_*.rsml'}
PX=2.54/300
def poly(r):
    g=r.find('geometry/polyline'); return None if g is None else np.array([[float(q.get('x')),float(q.get('y'))] for q in g])*PX
rows=[]; bad=[]
for day,pat in DAYS.items():
    for f in sorted(glob.glob(pat)):
        try: t=ET.parse(f).getroot()
        except Exception as e: bad.append((day,f)); continue
        num=os.path.basename(f).split('_')[-1].split('.')[0]
        for pl in t.iter('plant'):
            for mr in pl.findall('root'):
                P=poly(mr)
                if P is None or len(P)<2: continue
                cum=np.r_[0,np.cumsum(np.hypot(*np.diff(P,axis=0).T))]
                rows.append(dict(day=day,num=num,label=mr.get('label'),id=mr.get('ID'),parent=None,order=0,length=cum[-1],ins=np.nan,x0=P[0,0],y0=P[0,1],x1=P[-1,0],y1=P[-1,1]))
                for lr in mr.findall('root'):
                    Q=poly(lr)
                    if Q is None or len(Q)<2: continue
                    i=np.argmin(np.hypot(*(P-Q[0]).T))
                    rows.append(dict(day=day,num=num,label=mr.get('label'),id=lr.get('ID'),parent=mr.get('ID'),order=1,length=np.hypot(*np.diff(Q,axis=0).T).sum(),ins=cum[i],x0=Q[0,0],y0=Q[0,1],x1=Q[-1,0],y1=Q[-1,1]))
D=pd.DataFrame(rows); D.to_csv(os.path.join(KP.KRS_IN,'timeseries_roots.csv'),index=False)
print(D.groupby(['day','order']).size().unstack(), '\nunreadable files:', len(bad))
