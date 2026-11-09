"""Half-and-half cross-section icons: left = Emily Miller's photo, right = GRANAR v5 anatomy (same µm scale),
barrier walls colour-coded. Needs the photos (paths.IMAGES) and the outputs of 1_granar and 2_mecha.
Writes <outputs>/figures/icons/icon_<frag>.png and icons_info.json."""
import os, sys; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from paths import *
import re, json, numpy as np, pandas as pd
from PIL import Image, ImageDraw, ImageFilter
from scipy.spatial import ConvexHull
_IS=json.load(open(ICON_SECTIONS))['sections']
C={int(i):[int(i),v['segment'],os.path.join(IMAGES,v['image']),v['centre_x_px'],v['centre_y_px'],v['radius_px']] for i,v in _IS.items()}
PICK={'MP':17,'MM':8,'MT':0,'MT0':0,'LP':58,'LM':47,'LT':30}   # section index in the image set (chosen for clarity, root size near segment median)
BAR={'MP':(1,1),'MM':(1,1),'MT':(0,1),'MT0':(0,0),'LP':(1,0),'LM':(1,0),'LT':(0,0)}
EXO,ENDO=(237,161,0),(227,73,72)
FILL={'epidermis':(233,220,198),'exodermis':(240,214,160),'cortex':(241,233,214),'endodermis':(244,205,196),'passage_cell':(244,205,196),
      'pericycle':(226,220,238),'stele':(218,213,234),'phloem':(196,217,178),'companion_cell':(196,217,178),'cambium':(196,217,178),
      'xylem':(52,86,140),'metaxylem':(52,86,140),'inter_cellular_space':(255,255,255),'aerenchyma':(255,255,255)}
SI={int(i):dict(row=v['emily_row'],seg=v['segment'],file=v['image']) for i,v in _IS.items()}
X=pd.read_excel(SHEET); X.insert(0,'MPn',range(1,len(X)+1))
rd=dict(zip(X.MPn,2*np.sqrt(pd.to_numeric(X['Root Area (pixels)'],errors='coerce')*1.03214286e-7/np.pi)*1000))
K=pd.read_csv(COND)
def parse(path):
    t=open(path).read()
    names={}
    for gid,nm in re.findall(r'<group id="(\d+)" name="(\w+)"',t.split('<cellgroups>')[1].split('</cellgroups>')[0]): names.setdefault(int(gid),[]).append(nm)
    cells_txt=t.split('<cells')[1].split('</cells>')[0]
    cells=[(int(g),[int(w) for w in re.findall(r'<wall id="\s*(\d+)"',body)]) for g,body in re.findall(r'<cell id="\s*\d+" group="(\d+)"[^>]*>(.*?)</cell>',cells_txt,re.S)]
    walls={}
    for wid,body in re.findall(r'<wall id="\s*(\d+)"[^>]*>(.*?)</wall>',t.split('</cells>')[1],re.S):
        walls[int(wid)]=[(float(x),float(y)) for x,y in re.findall(r'x="\s*([\d.eE+-]+)" y="\s*([\d.eE+-]+)"',body)]
    return names,cells,walls
GS=pd.read_csv(GRANAR_SUMMARY).set_index('id')
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu
def recentre(f,cx,cy,r):
    w=int(2.2*r); im=np.asarray(Image.open(f).convert('RGB').crop((cx-w,cy-w,cx+w,cy+w))).astype(float)
    p=ndi.gaussian_filter((im[...,0]+im[...,2])/2-im[...,1],6)+ndi.gaussian_filter(np.abs(im.mean(2)-ndi.gaussian_filter(im.mean(2),3)),6)
    m=p>threshold_otsu(p); lab,n=ndi.label(m); k=np.argmax(ndi.sum(m,lab,range(1,n+1)))+1; b=ndi.binary_fill_holes(lab==k)
    yy,xx=ndi.center_of_mass(b); return cx-w+xx, cy-w+yy
def celltype(g,names,ncell_name):
    n=names.get(g,['?']); return n[0] if len(n)==1 else ncell_name.get(g,n[0])
info={}
for frag,i in PICK.items():
    s=SI[i]; row=s['row']; seg=s['seg']
    k=K[(K.anatomy.str.startswith(f'MP{row}_'))&(K.fragment==frag)].copy()
    med=K[K.fragment==frag].kr_cm_hPa_d.median(); k['dev']=(np.log(k.kr_cm_hPa_d)-np.log(med)).abs(); anat=k.sort_values('dev').anatomy.iloc[0]
    names,cells,walls=parse(os.path.join(ANAT,anat+'_root.xml'))
    pts=np.array([p for w in walls.values() for p in w]); cen=np.median(pts,0)
    ga=GS.loc[anat,'area_all']; Rg=np.sqrt(ga/np.pi)*1000   # µm, from the anatomy's own area (root hairs excluded)
    Re=rd[row]/2
    N=900; sc=(N/2*0.94)/Rg   # each half is drawn at its own scale so both roots fill the circle
    # right half: GRANAR
    g=Image.new('RGB',(N,N),(255,255,255)); d=ImageDraw.Draw(g)
    def P(x,y): return ((x-cen[0])*1000*sc+N/2, N/2-(y-cen[1])*1000*sc)
    gname={1:'exodermis',2:'epidermis',3:'endodermis',4:'cortex',5:'stele',13:'xylem',16:'pericycle',12:'companion_cell',11:'phloem'}
    exo_walls,endo_walls=set(),set()
    for gid,ws in cells:
        tp=gname.get(gid,'cortex'); q=[P(*p) for w in ws for p in walls.get(w,[])]
        if len(q)<3: continue
        q=np.array(q); 
        try: h=ConvexHull(q); poly=[tuple(q[j]) for j in h.vertices]
        except Exception: continue
        d.polygon(poly,fill=FILL.get(tp,(240,240,240)),outline=(120,110,100))
        if tp=='exodermis': exo_walls.update(ws)
        if tp=='endodermis': endo_walls.update(ws)
    exo,endo=BAR[frag]
    for ws,on,col in ((exo_walls,exo,EXO),(endo_walls,endo,ENDO)):
        if not on: continue
        for w in ws:
            q=[P(*p) for p in walls.get(w,[])]
            if len(q)>=2: d.line(q,fill=col,width=4)
    # left half: photo at the same µm scale
    c=C[i]; f=c[2]; cx,cy=recentre(f,c[3],c[4],Re/0.32127); um_px=0.32127
    half=Re/um_px/0.94; ph=Image.open(f).convert('RGB').crop((cx-half,cy-half,cx+half,cy+half)).resize((N,N),Image.LANCZOS)
    # compose
    out=Image.new('RGBA',(N,N),(0,0,0,0)); m=Image.new('L',(N,N),0); ImageDraw.Draw(m).ellipse([2,2,N-3,N-3],fill=255)
    left=Image.new('L',(N,N),0); ImageDraw.Draw(left).rectangle([0,0,N//2,N],fill=255)
    comp=Image.composite(ph,g,left); out.paste(comp,(0,0),m)
    dd=ImageDraw.Draw(out); dd.line([(N//2,4),(N//2,N-4)],fill=(255,255,255,255),width=5); dd.ellipse([2,2,N-3,N-3],outline=(90,90,90,255),width=4)
    out.save(os.path.join(ICONS,f'icon_{frag}.png'))
    info[frag]=dict(section_image=s['file'],emily_row=row,anatomy=anat,emily_root_d_um=round(2*Re),granar_root_d_um=round(2*Rg))
json.dump(info,open(os.path.join(ICONS,'icons_info.json'),'w'),indent=1); print(json.dumps(info,indent=1))
