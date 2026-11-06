"""X4b: IFN vs CDKN2A/MTAP classes within a-priori positive-control cancers, per cancer and pooled."""
import sys, warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'scripts')
from common import *
import statsmodels.formula.api as smf
m = pd.read_pickle(P+'panc.pkl').dropna(subset=['purity','logitLF','burden_ex9p']); POS = ['LUAD','LUSC','HNSC','SKCM']; rows = []
ifn = ['g_IFNB1','g_IFNE','g_IFNW1','g_IFNA1','g_IFNA2','g_IFNA13']
for thr in [-1,-2]:
    il = (m[ifn] <= thr).mean(axis=1) >= 0.5; cd = m.g_CDKN2A <= thr; mt = m.g_MTAP <= thr
    m['cls'] = np.select([il&cd&mt, (~il)&cd, ~(il|cd)], ['C1','C2','C4'], 'other')
    for lab, z in [('pooled positive controls', m[m.ct.isin(POS)]), ('pooled all other cancers', m[~m.ct.isin(POS) & (m.ct!='BRCA')])] + [(c, m[m.ct==c]) for c in POS]:
        z = z[z.cls != 'other']
        if (z.cls=='C1').sum() < 10 or (z.cls=='C2').sum() < 10: rows.append(dict(thr=thr, set=lab, note='too few', nC1=int((z.cls=='C1').sum()), nC2=int((z.cls=='C2').sum()))); continue
        f = smf.ols('logitLF ~ C(cls, Treatment("C4")) + C(ct) + burden_ex9p + purity', z).fit(cov_type='HC3')
        for c in ['C1','C2']:
            k = f'C(cls, Treatment("C4"))[T.{c}]'; rows.append(dict(thr=thr, set=lab, cls=c, n_cls=int((z.cls==c).sum()), beta=f.params[k], lo=f.conf_int().loc[k,0], hi=f.conf_int().loc[k,1], p=f.pvalues[k]))
r = pd.DataFrame(rows); r.to_csv('results/x4b_classes_pos.csv', index=False); print(r.round(3).to_string())
