"""X14: robustness of the primary BRCA contrast (9p21.3/CDKN2A loss on logit LF, M3)."""
import sys, warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'scripts')
from common import *
import statsmodels.formula.api as smf
from scipy.stats import norm, rankdata
m = pd.read_pickle(P+'panc.pkl'); base = m[m.ct=='BRCA'].copy(); rows = []
def run(lab, d, x='g_CDKN2A', thr=-1, cov=('C(sub)','burden_ex9p','purity'), y='logitLF', cluster=None):
    d = d.copy(); d['x'] = (d[x] <= thr).astype(int)
    if cluster is None: f = fit_ols(d, y, 'x', list(cov)); rows.append(dict(spec=lab, **f))
    else:
        r = smf.ols(f'{y} ~ x + ' + ' + '.join(cov), d).fit(cov_type='cluster', cov_kwds={'groups': pd.factorize(d[cluster])[0]}); ci = r.conf_int().loc['x']; ci90 = r.conf_int(0.10).loc['x']
        rows.append(dict(spec=lab, beta=r.params['x'], se=r.bse['x'], lo95=ci[0], hi95=ci[1], lo90=ci90[0], hi90=ci90[1], p=r.pvalues['x'], n=int(r.nobs)))
four = base[base['sub'].isin(['LumA','LumB','Basal','Her2'])].dropna(subset=['purity','logitLF'])
run('Primary (shallow loss, 4 PAM50, M3)', four)
run('Deep deletion only (<= -2)', four, thr=-2)
run('+ AS, ploidy, WGD', four.dropna(subset=['AS','ploidy','WGD']), cov=('C(sub)','burden_ex9p','purity','AS','ploidy','WGD'))
run('Normal-like included', base[base['sub'].isin(['LumA','LumB','Basal','Her2','Normal'])].dropna(subset=['purity','logitLF']))
run('Missing PAM50 included as own level', base.assign(sub=base['sub'].fillna('NA')).dropna(subset=['purity','logitLF']))
run('Exclude LF > 0.6', four[four.LF <= 0.6])
q = four.copy(); q['qLF'] = norm.ppf((rankdata(q.logitLF)-0.5)/len(q)); run('Rank-normalised LF', q, y='qLF')
run('High-purity half only (purity >= median)', four[four.purity >= four.purity.median()])
run('Low-purity half only', four[four.purity < four.purity.median()])
four2 = four.copy(); four2['tss'] = four2.key.str[5:7]; run('Cluster-robust by tissue source site', four2, cluster='tss')
run('Burden = ABSOLUTE FGA instead of gene-count', four.dropna(subset=['FGA']), cov=('C(sub)','FGA','purity'))
run('No purity adjustment (subtype + burden)', four, cov=('C(sub)','burden_ex9p'))
run('Purity only (no burden)', four, cov=('C(sub)','purity'))
run('Outcome = LF/(1-purity), logit', four.assign(lfs=LOGIT((four.LF/(1-four.purity)).clip(0.01,0.99))), y='lfs')
run('IFN-cluster loss as exposure', four.assign(g_ifnc=four[['g_IFNB1','g_IFNE','g_IFNW1','g_IFNA1','g_IFNA2','g_IFNA13']].min(axis=1)), x='g_ifnc')
r = pd.DataFrame(rows); r['equiv_0.10'] = (r.lo90 > -0.10) & (r.hi90 < 0.10); r.to_csv('results/x14_robust.csv', index=False); print(r[['spec','n','beta','lo95','hi95','lo90','hi90','p','equiv_0.10']].round(3).to_string())
