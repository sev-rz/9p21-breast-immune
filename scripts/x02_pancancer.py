"""Pan-cancer: X1b purity audit, X8 positive control / per-cancer 9p21.3 effect + meta-analysis, X4 breakpoint classes."""
import sys, warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'scripts')
from common import *
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
R = 'results/'
m = pd.read_pickle(P+'panc.pkl'); m = m.dropna(subset=['purity','logitLF','FGA','AS','burden_ex9p']).copy()
m['zFGA'] = m.groupby('ct').FGA.transform(lambda s: (s-s.mean())/s.std()); m['zAS'] = m.groupby('ct').AS.transform(lambda s: (s-s.mean())/s.std())
POS = ['LUAD','LUSC','HNSC','SKCM']
# ---------- X1b purity audit
rows = []
for ct, z in m.groupby('ct'):
    if len(z) < 100: continue
    for expo in ['zFGA','zAS']:
        a = fit_ols(z,'logitLF',expo,[]); c = fit_ols(z,'logitLF',expo,['purity'])
        rows.append(dict(ct=ct, exposure=expo, n=len(z), beta_nopurity=a['beta'], p_nopurity=a['p'], beta_purity=c['beta'], p_purity=c['p'], lo_purity=c['lo95'], hi_purity=c['hi95'], rho_LF_purity=z[['LF','purity']].corr(method='spearman').iloc[0,1]))
x1b = pd.DataFrame(rows)
for e in ['zFGA','zAS']:
    i = x1b.exposure==e; x1b.loc[i,'q_nopurity'] = multipletests(x1b.loc[i,'p_nopurity'], method='fdr_bh')[1]; x1b.loc[i,'q_purity'] = multipletests(x1b.loc[i,'p_purity'], method='fdr_bh')[1]
x1b['shrink'] = x1b.beta_purity/x1b.beta_nopurity; x1b.to_csv(R+'x1b_purity_audit.csv', index=False)
for e in ['zFGA','zAS']:
    s = x1b[x1b.exposure==e]; neg = s[(s.beta_nopurity<0)&(s.q_nopurity<0.05)]
    print(f'X1b {e}: cancers {len(s)}; negative & FDR<0.05 without purity: {len(neg)}; of these still FDR<0.05 negative with purity: {int(((neg.beta_purity<0)&(neg.q_purity<0.05)).sum())}; median LF-purity rho {s.rho_LF_purity.median():.2f}')
# ---------- X8: per-cancer 9p21.3 loss effect
rows = []
for ct, z in m.groupby('ct'):
    z = z.copy(); z['x'] = (z.g_CDKN2A <= -1).astype(int); z['xd'] = (z.g_CDKN2A <= -2).astype(int)
    if ct == 'BRCA': z = z[z['sub'].isin(['LumA','LumB','Basal','Her2'])]; cov1 = ['C(sub)']; cov3 = ['C(sub)','burden_ex9p','purity']
    else: cov1 = []; cov3 = ['burden_ex9p','purity']
    nl = int(z.x.sum())
    if nl < 20 or nl > len(z)-20: rows.append(dict(ct=ct, n=len(z), n_loss=nl, note='too few')); continue
    f0 = fit_ols(z,'logitLF','x',[]); f1 = fit_ols(z,'logitLF','x',cov1) if cov1 else f0; f3 = fit_ols(z,'logitLF','x',cov3)
    row = dict(ct=ct, n=len(z), n_loss=nl, loss_freq=nl/len(z), b0=f0['beta'], b1=f1['beta'], b3=f3['beta'], se3=f3['se'], lo3=f3['lo95'], hi3=f3['hi95'], lo90=f3['lo90'], hi90=f3['hi90'], p0=f0['p'], p3=f3['p'], cls=classify(f1, f3), positive_control=ct in POS)
    if z.xd.sum() >= 15: fd = fit_ols(z,'logitLF','xd',cov3); row.update(n_deep=int(z.xd.sum()), b_deep=fd['beta'], lo_deep=fd['lo95'], hi_deep=fd['hi95'])
    rows.append(row)
x8 = pd.DataFrame(rows); ok = x8.b3.notna(); x8.loc[ok,'q3'] = multipletests(x8.loc[ok,'p3'], method='fdr_bh')[1]
x8.to_csv(R+'x8_percancer_9p.csv', index=False)
print(x8[['ct','n','n_loss','b0','b3','lo3','hi3','p3','q3','cls','positive_control']].round(3).to_string())
def dl(b_, se_):
    w = 1/se_**2; mu = (w*b_).sum()/w.sum(); Q = (w*(b_-mu)**2).sum(); k = len(b_); tau2 = max(0, (Q-(k-1))/(w.sum()-(w**2).sum()/w.sum())); w2 = 1/(se_**2+tau2); mu2 = (w2*b_).sum()/w2.sum(); se2 = (1/w2.sum())**.5
    return dict(k=k, beta=mu2, lo=mu2-1.96*se2, hi=mu2+1.96*se2, tau2=tau2, I2=max(0,(Q-(k-1))/Q) if Q>0 else 0)
v = x8[ok]; meta = {'all cancers': dl(v.b3.values, v.se3.values), 'positive controls': dl(v[v.positive_control].b3.values, v[v.positive_control].se3.values), 'all except BRCA': dl(v[v.ct!='BRCA'].b3.values, v[v.ct!='BRCA'].se3.values), 'all except positive controls & BRCA': dl(v[(~v.positive_control)&(v.ct!='BRCA')].b3.values, v[(~v.positive_control)&(v.ct!='BRCA')].se3.values)}
pd.DataFrame(meta).T.to_csv(R+'x8_meta.csv'); print(pd.DataFrame(meta).T.round(3))
# ---------- X4: breakpoint classes (gene-level calls; rows in genomic order IFN cluster -> MTAP -> CDKN2A/B)
ifn = [c for c in ['g_IFNB1','g_IFNE','g_IFNW1','g_IFNA1','g_IFNA2','g_IFNA13'] if c in m.columns]
rows = []; out = {}
for thr in [-1, -2]:
    ifn_lost = (m[ifn] <= thr).mean(axis=1) >= 0.5; cdk = m.g_CDKN2A <= thr; mt = m.g_MTAP <= thr
    cls = np.select([ifn_lost & cdk & mt, (~ifn_lost) & cdk, ifn_lost & (~cdk), ~(ifn_lost|cdk)], ['C1 IFN+MTAP+CDKN2A lost', 'C2 CDKN2A lost, IFN intact', 'C3 IFN lost, CDKN2A intact', 'C4 neither'], 'C5 other')
    out[thr] = cls; m[f'cls{thr}'] = cls
    tab = pd.crosstab(m.ct, m[f'cls{thr}']); tab.to_csv(R+f'x4_class_counts_thr{abs(thr)}.csv'); tot = tab.sum(); print(f'X4 thr {thr}:', tot.to_dict())
    for k in ['C2 CDKN2A lost, IFN intact','C3 IFN lost, CDKN2A intact']:
        print('   cancers with >=30 in', k, ':', int((tab[k]>=30).sum()) if k in tab else 0)
m['arm9p_loss'] = (m['al_9p'] >= 0.8).astype(int)
rows = []
for thr in [-1, -2]:
    z = m[~m.ct.isin(['BRCA'])|m['sub'].isin(['LumA','LumB','Basal','Her2'])].copy(); z['cls'] = z[f'cls{thr}']
    z = z[z.cls != 'C5 other']
    import statsmodels.formula.api as smf
    for refcls in ['C4 neither']:
        f = smf.ols('logitLF ~ C(cls, Treatment("C4 neither")) + C(ct) + burden_ex9p + purity', z).fit(cov_type='HC3')
        for k in f.params.index:
            if 'cls' in k: rows.append(dict(thr=thr, cls=k.split('[T.')[1][:-1], beta=f.params[k], lo=f.conf_int().loc[k,0], hi=f.conf_int().loc[k,1], p=f.pvalues[k], n=int((z.cls==k.split('[T.')[1][:-1]).sum())))
    # restricted to no 9p arm loss
    zz = z[z.arm9p_loss==0]; f = smf.ols('logitLF ~ C(cls, Treatment("C4 neither")) + C(ct) + burden_ex9p + purity', zz).fit(cov_type='HC3')
    for k in f.params.index:
        if 'cls' in k: rows.append(dict(thr=thr, cls=k.split('[T.')[1][:-1]+' | no 9p arm loss', beta=f.params[k], lo=f.conf_int().loc[k,0], hi=f.conf_int().loc[k,1], p=f.pvalues[k], n=int((zz.cls==k.split('[T.')[1][:-1]).sum())))
x4 = pd.DataFrame(rows); x4.to_csv(R+'x4_class_effects.csv', index=False); print(x4.round(3).to_string())
