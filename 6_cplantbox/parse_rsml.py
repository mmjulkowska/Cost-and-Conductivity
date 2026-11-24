import cpb_paths as CP   # paths; changes into outputs/cplantbox/
import xml.etree.ElementTree as ET, numpy as np, pandas as pd, os, sys
D=CP.RSML_DIR+'/'
sel=pd.read_csv(CP.RSML_SELECTED)
PX=2.54/300  # cm per px (300 dpi)
def plen(p): return np.sum(np.hypot(*np.diff(p,axis=0).T))
rows=[];lat=[]
for _,s in sel.iterrows():
    f=CP.rsml(s.rsml)
    t=ET.parse(f).getroot()
    for pl in t.iter('plant'):
        for mr in pl.findall('root'):
            lab=mr.get('label')
            if not s.arbor_id.endswith(lab): continue
            P=np.array([[float(q.get('x')),float(q.get('y'))] for q in mr.find('geometry/polyline')])*PX
            L=plen(P); cum=np.r_[0,np.cumsum(np.hypot(*np.diff(P,axis=0).T))]
            mrang=np.degrees(np.arctan2(P[-1,0]-P[0,0],P[-1,1]-P[0,1]))
            nl=0
            for lr in mr.findall('root'):
                Q=np.array([[float(q.get('x')),float(q.get('y'))] for q in lr.find('geometry/polyline')])*PX
                if len(Q)<2: continue
                l=plen(Q); i=np.argmin(np.hypot(*(P-Q[0]).T)); ins=cum[i]
                # direction of MR at insertion and of lateral first ~0.3 cm / overall
                j=min(i+3,len(P)-1); k=max(i-3,0); mdir=P[j]-P[k]
                k2=np.searchsorted(np.r_[0,np.cumsum(np.hypot(*np.diff(Q,axis=0).T))],min(0.3,l)); ldir=Q[min(k2,len(Q)-1)]-Q[0]
                ang=np.degrees(np.arccos(np.clip(np.dot(mdir,ldir)/np.linalg.norm(mdir)/np.linalg.norm(ldir),-1,1))) if np.linalg.norm(ldir)>0 else np.nan
                tip=Q[-1]-Q[0]; tipv=np.degrees(np.arctan2(abs(tip[0]),tip[1]))  # angle from vertical down, whole lateral
                lat.append(dict(arbor=s.arbor_id,type=s.type,cond=s.condition,MR=L,ins=ins,rel=ins/L,len=l,ins_angle=ang,tip_from_vertical=tipv))
                nl+=1
            rows.append(dict(arbor=s.arbor_id,type=s.type,cond=s.condition,MR=L,nLR=nl,MR_angle=mrang))
R=pd.DataFrame(rows); Lt=pd.DataFrame(lat)
R.to_csv('rsml_plants.csv',index=False); Lt.to_csv('rsml_laterals.csv',index=False)
print(R.groupby(['type','cond']).agg(n=('MR','size'),MR=('MR','median'),nLR=('nLR','median')).round(2))
Lt['zone']=pd.cut(Lt.rel,[0,.33,.66,1.01],labels=['top','mid','bottom'])
print(Lt[Lt.cond=='control'].groupby(['type','zone']).agg(n=('len','size'),len=('len','median'),ang=('ins_angle','median'),tipv=('tip_from_vertical','median')).round(2))
g=Lt[Lt.cond=='control'].groupby('type')
print(g.agg(lbuz=('ins','min'),lmax=('ins','max'),len=('len','median'),ang=('ins_angle','median')).round(2))
