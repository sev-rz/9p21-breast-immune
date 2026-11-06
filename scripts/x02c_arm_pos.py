"""X6b: focal CDKN2A loss vs whole-9p-arm loss in positive-control and other cancers (pooled, cancer FE)."""
import sys, warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'scripts')
from common import *
import statsmodels.formula.api as smf
m = pd.read_pickle(P+'panc.pkl').dropna(subset=['purity','logitLF','burden_ex9p']); POS = ['LUAD','LUSC','HNSC','SKCM']
m['gene'] = (m.g_CDKN2A <= -1).astype(int); m['arm'] = (m.al_9p >= 0.8).astype(int); rows = []
for lab, z in [('pooled positive controls', m[m.ct.isin(POS)]), ('pooled other cancers (no BRCA)', m[~m.ct.isin(POS) & (m.ct!='BRCA')])] + [(c, m[m.ct==c]) for c in POS]:
    z = z.copy(); z['focal'] = ((z.gene==1)&(z.arm==0)).astype(int)
    for mod, f in [('gene loss alone','gene'),('arm loss alone','arm'),('gene | arm','gene + arm'),('arm | gene','arm + gene')]:
        fm = smf.ols(f'logitLF ~ {f} + C(ct) + burden_ex9p + purity', z).fit(cov_type='HC3'); key = f.split(' ')[0] if mod in ('gene loss alone','arm loss alone','gene | arm') else 'arm'
        k = 'gene' if mod in ('gene loss alone','gene | arm') else 'arm'; rows.append(dict(set=lab, model=mod, n=int(fm.nobs), n_gene=int(z.gene.sum()), n_arm=int(z.arm.sum()), n_focal=int(z.focal.sum()), beta=fm.params[k], lo=fm.conf_int().loc[k,0], hi=fm.conf_int().loc[k,1], p=fm.pvalues[k]))
r = pd.DataFrame(rows); r.to_csv('results/x6b_arm_pan.csv', index=False); print(r.round(3).to_string())
