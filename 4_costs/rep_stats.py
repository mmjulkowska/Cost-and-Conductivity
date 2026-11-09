"""Cost vs conductivity correlations per replicate (Spearman, mixed model with section as random effect, repeated-measures r). Outputs in <outputs>/costs/."""
import os, sys; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from paths import *
os.chdir(COSTS)  # cost tables live in <outputs>/costs/
import pandas as pd, numpy as np, warnings; warnings.filterwarnings('ignore')
from scipy import stats
import statsmodels.formula.api as smf
P=pd.read_csv('replicate_level_costs_conductivity.csv')
COST=[('C_ng_mm','Building C'),('maintC_ng_mm_d','Maintenance C'),('N_ng_mm','Total N')]
out=[]
for cv,cn in COST:
    for kv in ('kr','Kx'):
        d=P[[cv,kv,'section','segment']].dropna().copy(); d['x']=np.log10(d[cv]); d['y']=np.log10(d[kv])
        rho,ps=stats.spearmanr(d[cv],d[kv]); r,pr=stats.pearsonr(d.x,d.y)
        # (2) mixed model: random intercept per cross-section
        m1=smf.mixedlm('y ~ x',d,groups=d.section).fit(reml=True)
        # (2b) + segment fixed effect (removes segment/barrier differences)
        m2=smf.mixedlm('y ~ x + C(segment)',d,groups=d.section).fit(reml=True)
        # (3) repeated-measures correlation within cross-section (Bakdash & Marusich 2017): ANCOVA y ~ section + x
        a=smf.ols('y ~ C(section) + x',d).fit(); t=a.tvalues['x']; dfr=a.df_resid; rrm=np.sign(t)*np.sqrt(t**2/(t**2+dfr))
        out.append(dict(cost=cn,conductivity=kv,n_replicates=len(d),n_sections=d.section.nunique(),
            spearman_rho=rho,spearman_p=ps,pearson_r_log=r,pearson_p=pr,
            mixed_slope=m1.params['x'],mixed_se=m1.bse['x'],mixed_p=m1.pvalues['x'],
            mixed_seg_slope=m2.params['x'],mixed_seg_se=m2.bse['x'],mixed_seg_p=m2.pvalues['x'],
            rm_corr=rrm,rm_df=int(dfr),rm_p=a.pvalues['x'],rm_slope=a.params['x']))
R=pd.DataFrame(out); R.to_csv('correlations_replicate_level.csv',index=False)
pd.set_option('display.width',250)
print(R[['cost','conductivity','n_replicates','spearman_rho','spearman_p','mixed_slope','mixed_p','mixed_seg_slope','mixed_seg_p','rm_corr','rm_df','rm_p']].round(3).to_string(index=False))
