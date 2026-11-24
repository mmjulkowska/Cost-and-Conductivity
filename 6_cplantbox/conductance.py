"""Whole-root-system hydraulic conductance (Krs) of the CPlantBox architectures with the v5 segment conductivities.
Each root piece gets the kr / Kx of the sampled segment nearest to it by distance from its own tip:
  main root: 0-1.5 cm from tip = MT0 (true tip, no barrier), 1.5-4.5 = MT (B1), 4.5-7.5 = MM (B6), > 7.5 = MP (B6)
  laterals:  0-1.5 cm = LT (no barrier), 1.5-4.5 = LM (B10), > 4.5 = LP (B10)
(samples were taken at the tip and 3, 6, 9 cm from it; boundaries are the midpoints; TIP / STEP can be changed).
Values: medians per segment of Cost-and-Conductivity outputs/mecha/conductivities_per_anatomy.csv (v5, MECHA-dev 945b2be).
Radius per segment: Emily's measured cross-section area (anatomy_by_segment.csv, em_root_mm2), d = 2 sqrt(A/pi).
Scenario 'lat_endo': lateral LP/LM with an endodermal Casparian strip (B6 instead of B10; outputs/mecha/checks/bracket_laterals_kr.csv medians).
Krs (cm3 hPa-1 d-1) = water uptake / (soil - collar water potential), uniform soil, no gravity (Doussan-type network, radial
conductance 2 pi a l kr split equally to a segment's two nodes, axial Kx / l)."""
import cpb_paths as CP   # paths; changes into outputs/cplantbox/
import numpy as np, pandas as pd, os, plantbox as pb
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
C=pd.read_csv(CP.COND).groupby('fragment')[['kr_cm_hPa_d','Kx_cm4_hPa_d']].median()
A=pd.read_csv(CP.ANAT_SEG).set_index('segment')
BR=pd.read_csv(CP.BRACKET).set_index('fragment')
RAD={s:np.sqrt(A.loc[s,'em_root_mm2']/np.pi)/10 for s in A.index}; RAD['MT0']=RAD['MT']          # cm
CLASSES=['MP','MM','MT','MT0','LP','LM','LT']
COLOR={'MP':'#08306b','MM':'#2171b5','MT':'#6baed6','MT0':'#d1477a','LP':'#00441b','LM':'#41ab5d','LT':'#f39cc0'}
LABEL={'MP':'MP main proximal (B6)','MM':'MM main middle (B6)','MT':'MT main tip zone (B1)','MT0':'MT0 main true tip (no barrier)',
       'LP':'LP lateral proximal (B10)','LM':'LM lateral middle (B10)','LT':'LT lateral true tip (no barrier)'}
def props(scenario='v5'):
    kr={s:C.loc[s,'kr_cm_hPa_d'] for s in CLASSES}; kx={s:C.loc[s,'Kx_cm4_hPa_d'] for s in CLASSES}
    if scenario=='lat_endo':
        for s in ('LP','LM'): kr[s]=kr[s]*BR.loc[s,'ratio_median']
    return kr,kx
def classify(dist,is_main,TIP=1.5,STEP=3.0):
    if is_main: return 'MT0' if dist<TIP else 'MT' if dist<TIP+STEP else 'MM' if dist<TIP+2*STEP else 'MP'
    return 'LT' if dist<TIP else 'LM' if dist<TIP+STEP else 'LP'
def segments(p,TIP=1.5,STEP=3.0):
    """List of (node_i, node_j, length, class, xyz_i, xyz_j) for every root segment."""
    out=[]
    for o in p.getOrgans(pb.root):
        n=o.getNumberOfNodes()
        if n<2: continue
        X=np.array([[q.x,q.y,q.z] for q in o.getNodes()]); ids=[o.getNodeId(i) for i in range(n)]
        seg=np.linalg.norm(np.diff(X,axis=0),axis=1); cum=np.r_[0,np.cumsum(seg)]; L=cum[-1]; main=o.getParameter('subType')==1
        for k in range(n-1):
            if seg[k]<=1e-9: continue
            d=L-0.5*(cum[k]+cum[k+1]); out.append((ids[k],ids[k+1],seg[k],classify(d,main,TIP,STEP),X[k],X[k+1]))
    return out
def krs(p,scenario='v5',TIP=1.5,STEP=3.0,segs=None):
    S=segs or segments(p,TIP,STEP)
    if not S: return 0.,{c:0. for c in CLASSES},S
    kr,kx=props(scenario); idx={}
    for s in S:
        for q in (s[0],s[1]): idx.setdefault(q,len(idx))
    N=len(idx); r=[];c=[];v=[]; diag=np.zeros(N); rad=np.zeros(N); cls_rad=[]
    for i,j,l,cl,_,_ in S:
        a=idx[i]; b=idx[j]; g=kx[cl]/l; R=2*np.pi*RAD[cl]*l*kr[cl]
        r+=[a,b]; c+=[b,a]; v+=[-g,-g]; diag[a]+=g+R/2; diag[b]+=g+R/2; rad[a]+=R/2; rad[b]+=R/2; cls_rad.append((a,b,R,cl))
    # collar = the main root's first node (seed); psi_soil = 0, psi_collar = -1 -> Krs = uptake
    collar=idx[S[0][0]]
    M=coo_matrix((v+list(diag),(r+list(range(N)),c+list(range(N)))),shape=(N,N)).tocsr()
    rhs=np.zeros(N)            # sum_j g(psi_j - psi_i) + R(0 - psi_i) = 0  ->  M psi = 0 for free nodes
    free=np.ones(N,bool); free[collar]=False
    Mf=M[free][:,free]; rhs_f=-(M[free][:,[collar]].toarray().ravel())*(-1.)
    psi=np.zeros(N); psi[collar]=-1.; psi[free]=spsolve(Mf.tocsc(),rhs_f)
    up={k:0. for k in CLASSES}
    for a,b,R,cl in cls_rad: up[cl]+=R/2*(0-psi[a])+R/2*(0-psi[b])
    return sum(up.values()), up, S
if __name__=='__main__':
    print(pd.DataFrame({'kr':props()[0],'kr_lat_endo':props('lat_endo')[0],'Kx':props()[1],'radius_um':{k:RAD[k]*1e4 for k in CLASSES}}).loc[CLASSES])
