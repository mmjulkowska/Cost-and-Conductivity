"""Pool the SmartRoot root exports of PimpiBig1-6 with the panel types (input preparation for step 5).
Source: Tomato_Pareto/ZENODO_deposit/data/06_smartroot_root_exports (set SMARTROOT_EXPORTS to that folder).
Reads 01_panel_root_geometry.csv from KRS_IN; writes pooled_plants_by_day.csv (one row per plant and tracing day)
and pooled_laterals.csv (one row per lateral >= 0.1 cm and day) into KRS_IN."""
import pandas as pd, numpy as np, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import krs_paths as KP
E=os.path.join(os.path.expanduser(os.environ.get('SMARTROOT_EXPORTS', '~/mnt/Tomato_Pareto/ZENODO_deposit/data/06_smartroot_root_exports')), '')
FILES={'PimpiBig1':'pimpi_Big1_D4_D5.csv','PimpiBig2':'pimpi_Big2_D4_D5.csv','PimpiBig3':'pimpi_Big3_D4_D5.csv','PimpiBig4':'pimpi_Big4_D4_D5.csv','PimpiBig5':'pimpi_Big5_all_days.csv','PimpiBig6':'PimpBig6_all_days.csv'}
p=pd.read_csv(KP.inp('01_panel_root_geometry.csv'))
p['label']=p.arbor_id.str.extract(r'_\d{8}_\d{3}_(.+)$')[0].str.replace('-','_')
p=p[['experiment','label','architectural_type','condition','imaging_day']].drop_duplicates(['experiment','label'])
P=[];L=[]
for ex,f in FILES.items():
    d=pd.read_csv(E+f); d['day']=d.image.str.extract(r'day(\d)')[0].astype(int); d['date']=d.image.str.extract(r'day\d_(\d{8})')[0]
    m=d[d.root_order==0].copy(); m['label']=m.root_name.str.replace('-','_')
    m=m.sort_values('length').drop_duplicates(['day','label'],keep='last')            # duplicate tracings of the same plant on a day: keep the longest main root
    lat=d[d.root_order==1].merge(m[['root','label','day','length']].rename(columns={'root':'parent','length':'MR'}),on=['parent','day'])
    lat=lat[lat.length>0.1]
    g=lat.groupby(['label','day']); agg=pd.DataFrame({'nLR':g.size(),'totLR':g.length.sum()}).reset_index()
    pl=m[['label','day','date','length']].rename(columns={'length':'MR'}).merge(agg,on=['label','day'],how='left').fillna({'nLR':0,'totLR':0})
    pl['experiment']=ex; lat['experiment']=ex
    P.append(pl); L.append(lat[['experiment','label','day','root','length','insertion_position','insertion_angle','MR']])
P=pd.concat(P).merge(p,on=['experiment','label'],how='inner'); L=pd.concat(L).merge(p,on=['experiment','label'],how='inner')
P.to_csv(os.path.join(KP.KRS_IN,'pooled_plants_by_day.csv'),index=False); L.to_csv(os.path.join(KP.KRS_IN,'pooled_laterals.csv'),index=False)
c=P[P.condition=='control']
print(c.groupby(['experiment','day']).date.first().unstack())
print(c.groupby(['architectural_type','experiment']).label.nunique().unstack())
print('median MR by type/exp/day'); print(c.pivot_table(index=['architectural_type','experiment'],columns='day',values='MR',aggfunc='median').round(1))
