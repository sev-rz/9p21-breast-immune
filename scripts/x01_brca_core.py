"""BRCA core experiments: X1a (variance decomposition), X2 (attenuation ladder + rival comparison), X3 (dose/gain), X6 (focal vs arm), X7 (subtype), X9 (AIPW, E-value, equivalence)."""
import sys, itertools, warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'scripts')
from common import *
import statsmodels.formula.api as smf, statsmodels.api as sm
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold
R = 'results/'; rng = np.random.default_rng(20261016)
m = pd.read_pickle(P+'panc.pkl'); b = m[(m.ct=='BRCA') & m['sub'].isin(['LumA','LumB','Basal','Her2'])].copy()
b = b.dropna(subset=['purity','AS','ploidy','WGD','logitLF']).reset_index(drop=True)
print('BRCA modelled n =', len(b), b['sub'].value_counts().to_dict())
ARM = {'CDKN2A':'9p','MTAP':'9p','IFNA1':'9p','JAK2':'9p','B2M':'15q','STAT1':'2q','JAK1':'1p','HLA-A':'6p','TAP1':'6p','IFNGR1':'6q','IRF1':'5q','PTEN':'10q','RB1':'13q','MYC':'8q'}
loss = lambda s: (b['g_'+s] <= -1).astype(int)
b['ifncl'] = ((b[['g_IFNB1','g_IFNE','g_IFNW1','g_IFNA1','g_IFNA2','g_IFNA13']] <= -1).any(axis=1)).astype(int)
UNITS = {'9p21.3 (CDKN2A)': (loss('CDKN2A'), ['9p']), 'MTAP': (loss('MTAP'), ['9p']), 'IFN cluster (any)': (b.ifncl, ['9p']),
 'Antigen presentation (B2M/HLA-A/B/C/TAP1/2)': ((b[['g_B2M','g_HLA-A','g_HLA-B','g_HLA-C','g_TAP1','g_TAP2']]<=-1).any(axis=1).astype(int), ['15q','6p']),
 'IFN signalling (JAK1/2,STAT1,IRF1,IFNGR1/2)': ((b[['g_JAK1','g_JAK2','g_STAT1','g_IRF1','g_IFNGR1','g_IFNGR2']]<=-1).any(axis=1).astype(int), ['1p','9p','2q','5q','6q','21q']),
 'B2M': (loss('B2M'), ['15q']), 'HLA-A': (loss('HLA-A'), ['6p']), 'TAP1': (loss('TAP1'), ['6p']), 'STAT1': (loss('STAT1'), ['2q']), 'JAK1': (loss('JAK1'), ['1p']), 'JAK2': (loss('JAK2'), ['9p']), 'IFNGR1': (loss('IFNGR1'), ['6q']),
 'PTEN (control)': (loss('PTEN'), ['10q']), 'RB1 (control)': (loss('RB1'), ['13q']), 'MYC (control)': (loss('MYC'), ['8q'])}
# ---------- X1a: Shapley decomposition of R2 for logit(LF)
groups = {'Subtype':['C(sub)'], 'Purity':['purity'], 'Burden (FGA, 9p-out)':['burden_ex9p'], 'AS+ploidy+WGD':['AS','ploidy','WGD']}
def r2(gs):
    if not gs: return 0.0
    f = 'logitLF ~ ' + ' + '.join(sum([groups[g] for g in gs], [])); return smf.ols(f, b).fit().rsquared
names = list(groups); K = len(names); sh = {}
from math import factorial
for g in names:
    others = [x for x in names if x != g]; v = 0
    for k in range(K):
        for S in itertools.combinations(others, k):
            w = factorial(k)*factorial(K-k-1)/factorial(K); v += w*(r2(list(S)+[g]) - r2(list(S)))
    sh[g] = v
pd.Series(sh).rename('shapley_R2').to_frame().assign(total_R2=r2(names)).to_csv(R+'x1a_shapley.csv'); print('X1a', {k: round(v,3) for k,v in sh.items()}, 'total', round(r2(names),3))
b[['logitLF','purity','burden_ex9p','AS','ploidy','FGA','LF']].corr(method='spearman').round(3).to_csv(R+'x1a_corr.csv')
# ---------- X2: attenuation ladder
rows = []
for name,(x,arms) in UNITS.items():
    b['x'] = x.values; bex = burden_ex(b, arms); b['bex'] = bex
    if b.x.sum() < 20: continue
    spec = {'M0 marginal': [], 'M1 +subtype': ['C(sub)'], 'M2 +burden(arm-out)': ['C(sub)','bex'], 'M3 +purity': ['C(sub)','bex','purity'], 'M4 +AS,ploidy,WGD': ['C(sub)','bex','purity','AS','ploidy','WGD']}
    fits = {k: fit_ols(b, 'logitLF', 'x', c) for k,c in spec.items()}
    cls = classify(fits['M1 +subtype'], fits['M3 +purity'])
    for k,v in fits.items(): rows.append(dict(unit=name, model=k, n_loss=int(b.x.sum()), **v, cls=cls if k=='M3 +purity' else ''))
x2 = pd.DataFrame(rows); x2.to_csv(R+'x2_ladder.csv', index=False)
print(x2[x2.model.isin(['M0 marginal','M1 +subtype','M3 +purity'])][['unit','model','n_loss','beta','lo95','hi95','p','cls']].round(3).to_string())
# ---------- X2b rival comparison (H2): R1/R2 identical in BRCA?
b['cdkn2a_loss'] = loss('CDKN2A').values; b['arm9p_loss'] = (b['al_9p'] >= 0.8).astype(int); b['ifn_loss'] = b.ifncl
print('IFN-cluster loss w/o CDKN2A loss:', int(((b.ifn_loss==1)&(b.cdkn2a_loss==0)).sum()), '| CDKN2A loss w/o IFN loss:', int(((b.ifn_loss==0)&(b.cdkn2a_loss==1)).sum()), '| both', int(((b.ifn_loss==1)&(b.cdkn2a_loss==1)).sum()), '| 9p arm loss n', int(b.arm9p_loss.sum()), '| arm loss & CDKN2A loss', int(((b.arm9p_loss==1)&(b.cdkn2a_loss==1)).sum()))
base = 'logitLF ~ C(sub) + purity + burden_ex9p + AS + ploidy + WGD'
mods = {'R4 base (subtype+purity+burden)': base, 'R2/R1 + CDKN2A loss': base+' + cdkn2a_loss', 'R3 + 9p arm loss': base+' + arm9p_loss', 'R3+R2 both': base+' + arm9p_loss + cdkn2a_loss', 'R1 + IFN-cluster loss': base+' + ifn_loss'}
cvres = []
for k,f in mods.items():
    r = smf.ols(f, b).fit(); sc = []
    for rep in range(5):
        for tr,te in KFold(10, shuffle=True, random_state=rep).split(b):
            rr = smf.ols(f, b.iloc[tr]).fit(); pr = rr.predict(b.iloc[te]); y = b.logitLF.iloc[te]
            sc.append(1 - ((y-pr)**2).sum()/((y-b.logitLF.iloc[tr].mean())**2).sum())
    cvres.append(dict(model=k, k=int(r.df_model), AIC=r.aic, BIC=r.bic, R2=r.rsquared, cvR2=np.mean(sc), cvR2_sd=np.std(sc)))
pd.DataFrame(cvres).to_csv(R+'x2b_rivals.csv', index=False); print(pd.DataFrame(cvres).round(4).to_string())
# ---------- X3: dose response + gain symmetry (CDKN2A) and loss vs gain per unit
b['call'] = b['g_CDKN2A'].clip(upper=1).astype(int).astype(str)
r = smf.ols('logitLF ~ C(call, Treatment("0")) + C(sub) + burden_ex9p + purity', b).fit(cov_type='HC3')
dr = pd.DataFrame({'coef':r.params.filter(like='call'), 'lo':r.conf_int().loc[r.params.filter(like='call').index,0], 'hi':r.conf_int().loc[r.params.filter(like='call').index,1]}); dr['n'] = [int((b.call==c[-3:-2] if False else b.call==c.split('[T.')[1][:-1]).sum()) for c in dr.index]
dr.to_csv(R+'x3_dose_cdkn2a.csv'); print(dr.round(3))
rows = []
for name,(x,arms) in UNITS.items():
    pass
gains = {'9p21.3 (CDKN2A)': ('CDKN2A',['9p']), 'B2M':('B2M',['15q']), 'STAT1':('STAT1',['2q']), 'JAK1':('JAK1',['1p']), 'HLA-A':('HLA-A',['6p']), 'PTEN (control)':('PTEN',['10q']), 'RB1 (control)':('RB1',['13q']), 'MYC (control)':('MYC',['8q']), 'TAP1':('TAP1',['6p']), 'JAK2':('JAK2',['9p'])}
for name,(gn,arms) in gains.items():
    b['bex'] = burden_ex(b, arms); b['xl'] = (b['g_'+gn] <= -1).astype(int); b['xg'] = (b['g_'+gn] >= 1).astype(int)
    for lab,xx in [('loss','xl'),('gain','xg')]:
        if b[xx].sum() >= 20 and b[xx].sum() <= len(b)-20:
            rows.append(dict(unit=name, kind=lab, n_events=int(b[xx].sum()), **fit_ols(b,'logitLF',xx,['C(sub)','bex','purity'])))
pd.DataFrame(rows).to_csv(R+'x3_gain_loss.csv', index=False); print(pd.DataFrame(rows)[['unit','kind','n_events','beta','lo95','hi95','p']].round(3).to_string())
# ---------- X6: focal vs arm
rows = []
for name,gn,arm in [('9p21.3 / 9p','CDKN2A','9p'),('B2M / 15q','B2M','15q'),('STAT1 / 2q','STAT1','2q'),('JAK1 / 1p','JAK1','1p')]:
    b['bex'] = burden_ex(b,[arm]); b['gl'] = (b['g_'+gn] <= -1).astype(int); b['al'] = (b['al_'+arm] >= 0.8).astype(int)
    b['focal'] = ((b.gl==1)&(b.al==0)).astype(int)
    cov = ['C(sub)','bex','purity']
    for lab,(x,extra) in {'gene loss alone':('gl',[]), 'arm loss alone':('al',[]), 'gene loss | arm loss':('gl',['al']), 'arm loss | gene loss':('al',['gl'])}.items():
        rows.append(dict(pair=name, model=lab, n_gene=int(b.gl.sum()), n_arm=int(b.al.sum()), n_focal=int(b.focal.sum()), **fit_ols(b,'logitLF',x,cov+extra)))
pd.DataFrame(rows).to_csv(R+'x6_focal_arm.csv', index=False); print(pd.DataFrame(rows)[['pair','model','n_gene','n_arm','n_focal','beta','lo95','hi95','p']].round(3).to_string())
# ---------- X7: subtype heterogeneity
rows = []
for name in ['9p21.3 (CDKN2A)','Antigen presentation (B2M/HLA-A/B/C/TAP1/2)','IFN signalling (JAK1/2,STAT1,IRF1,IFNGR1/2)','STAT1','B2M']:
    x,arms = UNITS[name]; b['x'] = x.values; b['bex'] = burden_ex(b, arms)
    for s in ['LumA','LumB','Basal','Her2']:
        z = b[b['sub']==s]; rows.append(dict(unit=name, subtype=s, n_loss=int(z.x.sum()), **fit_ols(z,'logitLF','x',['bex','purity'])))
    for lab,sel in [('Luminal (A+B)',['LumA','LumB']),('Basal',['Basal'])]:
        z = b[b['sub'].isin(sel)]; rows.append(dict(unit=name, subtype=lab, n_loss=int(z.x.sum()), **fit_ols(z,'logitLF','x',['C(sub)','bex','purity'])))
    full = smf.ols('logitLF ~ x*C(sub) + bex + purity', b).fit(); red = smf.ols('logitLF ~ x + C(sub) + bex + purity', b).fit()
    F = full.compare_f_test(red); lb = b[b['sub'].isin(['LumA','LumB','Basal'])].copy(); lb['lum'] = lb['sub'].isin(['LumA','LumB']).astype(int)
    f2 = smf.ols('logitLF ~ x*lum + bex + purity + C(sub)', lb).fit(cov_type='HC3')
    rows.append(dict(unit=name, subtype='interaction x*subtype (4-way F)', p=F[1], beta=np.nan)); rows.append(dict(unit=name, subtype='interaction luminal vs basal (HC3)', beta=f2.params['x:lum'], lo95=f2.conf_int().loc['x:lum',0], hi95=f2.conf_int().loc['x:lum',1], p=f2.pvalues['x:lum']))
x7 = pd.DataFrame(rows); x7.to_csv(R+'x7_subtype.csv', index=False); print(x7[['unit','subtype','n_loss','beta','lo95','hi95','p']].round(3).to_string())
# ---------- X9: AIPW + matching + equivalence + E-value for 9p21.3 loss
b['Tr'] = b.cdkn2a_loss; cov = pd.get_dummies(b['sub'], drop_first=True).astype(float).join(b[['purity','burden_ex9p','AS','ploidy','WGD']])
Z = ((cov - cov.mean())/cov.std()).to_numpy(); Y = b.logitLF.to_numpy(); T = b['Tr'].to_numpy()
def aipw(idx):
    z, y, t = Z[idx], Y[idx], T[idx]
    e = LogisticRegression(C=1e3, max_iter=1000).fit(z, t).predict_proba(z)[:,1].clip(0.02, 0.98)
    m1 = sm.OLS(y[t==1], sm.add_constant(z[t==1], has_constant='add')).fit(); m0 = sm.OLS(y[t==0], sm.add_constant(z[t==0], has_constant='add')).fit()
    zc = sm.add_constant(z, has_constant='add'); mu1 = m1.predict(zc); mu0 = m0.predict(zc)
    return np.mean(mu1 - mu0 + t*(y-mu1)/e - (1-t)*(y-mu0)/(1-e))
est = aipw(np.arange(len(b))); boot = np.array([aipw(rng.integers(0, len(b), len(b))) for _ in range(1000)])
lo, hi = np.percentile(boot, [2.5, 97.5]); lo90, hi90 = np.percentile(boot, [5, 95])
# nearest-neighbour matching on propensity (ATT, caliper 0.2 SD of logit PS)
ps = LogisticRegression(C=1e3, max_iter=1000).fit(Z, T).predict_proba(Z)[:,1]; lp = np.log(ps/(1-ps)); cal = 0.2*lp.std()
ti = np.where(T==1)[0]; ci = np.where(T==0)[0]; d = []
for i in ti:
    j = ci[np.argmin(np.abs(lp[ci]-lp[i]))]
    if abs(lp[j]-lp[i]) <= cal: d.append(Y[i]-Y[j])
sdres = smf.ols(base, b).fit().resid.std()
def evalue(beta, lo, hi):
    f = lambda d_: (lambda rr: rr + np.sqrt(rr*(rr-1)))(np.exp(0.91*abs(d_)/sdres) if abs(d_)>0 else 1.0)
    near = lo if beta > 0 else hi
    return f(beta), (1.0 if (lo <= 0 <= hi) else f(near))
ev = evalue(est, lo, hi)
res = dict(aipw=est, lo95=lo, hi95=hi, lo90=lo90, hi90=hi90, equiv_within_0p10=bool(lo90 > -0.10 and hi90 < 0.10), matched_att=float(np.mean(d)), matched_n=len(d), matched_se=float(np.std(d, ddof=1)/np.sqrt(len(d))), n_treated=int(T.sum()), n=len(b), evalue_est=ev[0], evalue_ci=ev[1], resid_sd=sdres)
pd.Series(res).to_csv(R+'x9_aipw_brca.csv'); print('X9', {k: (round(v,3) if isinstance(v, float) else v) for k,v in res.items()})
# covariate balance
bal = pd.DataFrame({'smd_raw': (cov[T==1].mean()-cov[T==0].mean())/cov.std()}); bal.to_csv(R+'x9_balance.csv'); print(bal.round(2).T.to_string())
