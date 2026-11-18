import pandas as pd, numpy as np, statsmodels.formula.api as smf, warnings; warnings.filterwarnings('ignore')
D='data/raw/'
l=pd.read_csv(D+'leuk.tsv',sep='\t',header=None,names=['ct','bc','LF']); l=l[l.ct=='BRCA']; l['s']=l.bc.str[:15]; l=l[l.s.str.endswith('-01')].drop_duplicates('s').set_index('s')
m=pd.read_csv(D+'mastercalls.txt',sep='\t'); m['s']=m['array'].str[:15]; m=m.drop_duplicates('s').set_index('s')
g=pd.read_csv(D+'cbp_gistic.csv',index_col=0); cs=pd.read_csv(D+'cbp_clinical.csv',index_col=0); pt=pd.read_csv(D+'cbp_patient.csv',index_col=0)
d=l[['LF']].join(g,how='inner').join(cs[['FRACTION_GENOME_ALTERED','ANEUPLOIDY_SCORE','MUTATION_COUNT']]).join(m[['purity','ploidy']])
d['sub']=pd.Series(d.index.str[:12],index=d.index).map(pt['SUBTYPE']).str.replace('BRCA_','')
d=d.rename(columns={'FRACTION_GENOME_ALTERED':'FGA','ANEUPLOIDY_SCORE':'AS'}); d['logit']=np.log(d.LF/(1-d.LF))
d['sub']=d['sub'].fillna('NA'); print(len(d)); 
d9=d[d['sub'].isin(['LumA','LumB','Basal','Her2'])].dropna(subset=['FGA','purity','AS']).copy()
print('n model',len(d9))
rows=[]
for gene in ['CDKN2A','MTAP','IFNA1','IFNE','B2M','JAK2','STAT1','TAP1','HLA-A','PTEN','RB1','IFNGR1','JAK1']:
    d9['x']=(d9[gene]<=-1).astype(int)
    r={'gene':gene,'n_loss':int(d9.x.sum())}
    for name,f in [('M0','logit~x'),('M1','logit~x+C(sub)'),('M2','logit~x+C(sub)+FGA'),('M3','logit~x+C(sub)+FGA+purity')]:
        fit=smf.ols(f,d9).fit(); r[name]=round(fit.params['x'],3); r[name+'p']=fit.pvalues['x']
    rows.append(r)
o=pd.DataFrame(rows); pd.set_option('display.width',250); print(o.round(4).to_string())
# per-subtype for CDKN2A
for s in ['LumA','LumB','Basal','Her2']:
    z=d9[d9['sub']==s].copy(); z['x']=(z['CDKN2A']<=-1).astype(int)
    fit=smf.ols('logit~x+FGA+purity',z).fit(); print(s,len(z),int(z.x.sum()),'coef',round(fit.params['x'],3),'p',round(fit.pvalues['x'],4))
print('corr FGA-x CDKN2A', d9[['FGA']].assign(x=(d9.CDKN2A<=-1)).groupby('x').FGA.median().round(3).to_dict())
