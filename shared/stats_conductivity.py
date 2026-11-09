"""One-way ANOVA + Tukey HSD on log10 values, with a compact letter display (helpers for 2_mecha/analyse.py and 4_costs/rep_costs.py)."""
import pandas as pd, numpy as np, sys, json
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.oneway import anova_oneway
ORDER=['MP','MM','MT','MT0','LP','LM','LT']
def cld(groups, sig):  # sig[(a,b)] True if significantly different; insert-and-absorb letters
    letters=[set(groups)]
    for (a,b),s in sig.items():
        if not s: continue
        new=[]
        for L in letters:
            if a in L and b in L: new+= [L-{a}, L-{b}]
            else: new.append(L)
        # absorb
        new=[L for L in new if L]; keep=[]
        for L in new:
            if not any(L<M for M in new if M is not L) and L not in keep: keep.append(L)
        letters=keep
    letters=sorted(letters,key=lambda L:min(groups.index(g) for g in L))
    out={g:'' for g in groups}
    for i,L in enumerate(letters):
        for g in L: out[g]+='abcdefghij'[i]
    return out
def analyse(df,var,groups,label):
    d=df[df.fragment.isin(groups)].dropna(subset=[var]).copy(); d['y']=np.log10(d[var])
    samples=[d.y[d.fragment==g].values for g in groups]
    F,p=stats.f_oneway(*samples); lev=stats.levene(*samples,center='median')
    w=anova_oneway(samples,use_var='unequal')
    tk=pairwise_tukeyhsd(d.y,d.fragment,alpha=0.05)
    t=pd.DataFrame(tk._results_table.data[1:],columns=tk._results_table.data[0])
    t['fold_change(group2/group1)']=10**t.meandiff.astype(float)
    sig={(r.group1,r.group2):bool(r.reject) for r in t.itertuples()}
    letters=cld(groups,sig)
    summ=d.groupby('fragment')[var].agg(n='count',mean='mean',sd='std',median='median').reindex(groups)
    summ['geo_mean']=d.groupby('fragment').y.mean().reindex(groups).pipe(lambda s:10**s)
    summ['tukey_letter']=[letters[g] for g in groups]
    dfb=len(groups)-1; dfw=len(d)-len(groups)
    res=dict(dataset=label,variable=var,groups=groups,F=float(F),df_between=dfb,df_within=dfw,p=float(p),
             levene_W=float(lev.statistic),levene_p=float(lev.pvalue),welch_F=float(w.statistic),welch_df=[float(x) for x in w.df],welch_p=float(w.pvalue))
    return summ,t,res
