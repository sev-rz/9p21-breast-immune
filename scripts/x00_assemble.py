import sys; sys.path.insert(0,'scripts')
from common import *
d = load_master(); X, g, keys = load_gistic(); tf, per, n, ntot, tot = burden(X, g, keys)
gi = {s:i for i,s in enumerate(g['Gene Symbol'])}
genes = ['CDKN2A','CDKN2B','MTAP','IFNB1','IFNE','IFNW1','IFNA1','IFNA2','IFNA13','B2M','HLA-A','HLA-B','HLA-C','TAP1','TAP2','JAK1','JAK2','STAT1','IRF1','IFNGR1','IFNGR2','PTEN','RB1','MYC','CD274','PDCD1LG2','TP53']
G = pd.DataFrame({'key':keys, 'burden_all':tf})
for s in genes:
    if s in gi: G['g_'+s] = X[gi[s]]
import json
G['alt_total'] = tot
arms = [a for a in pd.unique(g.arm.values) if isinstance(a, str) and a != 'nan']
for a in arms:
    ix = np.where(g.arm.values==a)[0]
    G['af_'+a] = (np.abs(X[ix])>=1).mean(0); G['al_'+a] = (X[ix]<=-1).mean(0); G['ag_'+a] = (X[ix]>=1).mean(0)
json.dump({'n':{a:int(v) for a,v in n.items()}, 'ntot':int(ntot)}, open(P+'arm_meta.json','w'))
# leave-9p-out burden
G['burden_ex9p'] = (tot - per['9p'])/(ntot - n['9p'])
m = d.merge(G, on='key', how='inner')
cs = pd.read_csv(D+'cbp_clinical.csv', index_col=0); pt = pd.read_csv(D+'cbp_patient.csv', index_col=0)
m['sub'] = m.pat.map(pt['SUBTYPE']).str.replace('BRCA_','')
m.to_pickle(P+'panc.pkl'); print(m.shape); print(m.groupby('ct').size().to_string())
print('missing purity', m.purity.isna().sum(), 'of', len(m))
print(m[m.ct=='BRCA']['sub'].value_counts(dropna=False).to_dict())
