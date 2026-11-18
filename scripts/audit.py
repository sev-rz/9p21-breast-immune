import pandas as pd, numpy as np
from scipy.stats import spearmanr
D='data/raw/'
l=pd.read_csv(D+'leuk.tsv',sep='\t',header=None,names=['ct','bc','LF'])
l=l[l.ct=='BRCA'].copy(); l['s']=l.bc.str[:15]; l['p']=l.bc.str[:12]
print('BRCA LF rows',len(l),'unique samples',l.s.nunique(),'tumor(-01)',l.s.str.endswith('-01').sum())
l=l[l.s.str.endswith('-01')].drop_duplicates('s')
print(l.LF.describe().round(3).to_dict())
m=pd.read_csv(D+'mastercalls.txt',sep='\t'); m['s']=m['array'].str[:15]
sg=pd.read_csv(D+'seg_scores.tsv',sep='\t'); sg['s']=sg.Sample.str[:15]
ab=pd.read_csv(D+'absolute_scores.tsv',sep='\t',index_col=0); ab['s']=ab.index.str[:15]
st=pd.read_csv(D+'subtype.tsv',sep='\t'); st=st[st['cancer.type']=='BRCA'].copy(); st['p']=st['pan.samplesID']
d=l.merge(m[['s','purity','ploidy','Genome doublings']].drop_duplicates('s'),on='s',how='left').merge(sg[['s','frac_altered','n_segs']].drop_duplicates('s'),on='s',how='left').merge(ab[['s','AS','LOH_frac_altered']].drop_duplicates('s'),on='s',how='left').merge(st[['p','Subtype_mRNA']],on='p',how='left')
print(len(d)); print(d.isna().sum().to_dict())
print(d.Subtype_mRNA.value_counts(dropna=False).to_dict())
for c in ['purity','frac_altered','AS','n_segs','ploidy']:
    x=d[['LF',c]].dropna(); print(c,'vs LF rho',round(spearmanr(x.LF,x[c])[0],3),len(x))
x=d[['purity','frac_altered']].dropna(); print('purity vs FGA rho',round(spearmanr(x.purity,x.frac_altered)[0],3))
print(d.groupby('Subtype_mRNA')[['LF','frac_altered','purity']].median().round(3))
import statsmodels.formula.api as smf
x=d.dropna(subset=['LF','frac_altered','purity'])
for f in ['LF~frac_altered','LF~frac_altered+purity']:
    r=smf.ols(f,x).fit(); print(f,'FGA coef',round(r.params['frac_altered'],3),'p',r.pvalues['frac_altered'],'R2',round(r.rsquared,3))
