"""X5: genome-wide cytoband scan in BRCA with stratified permutation null (Frisch-Waugh-Lovell), loss and gain; loss-vs-gain symmetry."""
import sys, warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'scripts')
from common import *
from statsmodels.stats.multitest import multipletests
R = 'results/'; rng = np.random.default_rng(7); P = 1000
X, g, keys = load_gistic(); m = pd.read_pickle(P+'panc.pkl')
b = m[(m.ct=='BRCA') & m['sub'].isin(['LumA','LumB','Basal','Her2'])].dropna(subset=['purity','logitLF']).reset_index(drop=True)
kidx = {k:i for i,k in enumerate(keys)}; cols = np.array([kidx[k] for k in b.key]); Xb = X[:, cols]
g = g.reset_index(drop=True); ok = g.Cytoband.notna().values & g.arm.notna().values
bands = []; B = []; order = []
for cb, ix in g[ok].groupby('Cytoband', sort=False).groups.items():
    ix = np.array(list(ix)); med = np.median(Xb[ix], axis=0); B.append(med); bands.append((cb, g.arm[ix[0]], ix[0]))
B = np.array(B).T  # n x bands (median call)
meta = pd.DataFrame(bands, columns=['band','arm','pos']); meta['nloss'] = (B<=-1).sum(0); meta['ngain'] = (B>=1).sum(0)
n = len(b); y = b.logitLF.to_numpy(); sub = pd.get_dummies(b['sub'], drop_first=True).astype(float).to_numpy()
strata = pd.factorize(b['sub'].astype(str) + '_' + pd.qcut(b.burden_all, 10, labels=False, duplicates='drop').astype(str) + '_' + pd.qcut(b.purity, 3, labels=False).astype(str))[0]
sidx = [np.where(strata==s)[0] for s in np.unique(strata)]
def permute():
    p = np.arange(n)
    for ix in sidx: p[ix] = rng.permutation(ix)
    return p
perms = [permute() for _ in range(P)]
res = []; maxnull = {'loss': np.zeros(P), 'gain': np.zeros(P)}; allnull = {}
for arm, ms in meta.groupby('arm'):
    bex = burden_ex(b, [arm]).to_numpy()
    out = {}
    for kind, ind, thr in [('m0', None, None)]: pass
    for adj, Zc in [('M1', np.column_stack([np.ones(n), sub])), ('M3', np.column_stack([np.ones(n), sub, bex, b.purity.to_numpy()]))]:
        Q, _ = np.linalg.qr(Zc); yr = y - Q@(Q.T@y)
        for kind, f in [('loss', lambda M: (M<=-1)), ('gain', lambda M: (M>=1))]:
            cand = ms.index[(f(B[:, ms.index]).sum(0) >= 50) & (f(B[:, ms.index]).sum(0) <= n-50)]
            if len(cand) == 0: continue
            Xm = f(B[:, cand]).astype(float); Xr = Xm - Q@(Q.T@Xm); ss = (Xr**2).sum(0)
            beta = (Xr*yr[:,None]).sum(0)/ss; rss = ((yr[:,None]-Xr*beta)**2).sum(0); se = np.sqrt(rss/(n-Zc.shape[1]-1)/ss); t = beta/se
            for j, c in enumerate(cand): out[(c, kind, adj)] = (beta[j], se[j], t[j])
            if adj == 'M3':
                nullt = np.zeros((P, len(cand)))
                for pi, p in enumerate(perms):
                    Xp = Xm[p]; Xpr = Xp - Q@(Q.T@Xp); ssp = (Xpr**2).sum(0); bp = (Xpr*yr[:,None]).sum(0)/ssp; rp = ((yr[:,None]-Xpr*bp)**2).sum(0); nullt[pi] = bp/np.sqrt(rp/(n-Zc.shape[1]-1)/ssp)
                maxnull[kind] = np.maximum(maxnull[kind], np.abs(nullt).max(1))
                for j, c in enumerate(cand):
                    res.append(dict(idx=c, kind=kind, perm_p=(1+(np.abs(nullt[:,j]) >= abs(t[j])).sum())/(P+1)))
rows = []
for (c, kind, adj), (bt, se, t) in out.items() if False else []: pass
# assemble from res + re-fit table
tab = pd.DataFrame(res); tab['band'] = meta.loc[tab.idx, 'band'].values; tab['arm'] = meta.loc[tab.idx, 'arm'].values; tab['pos'] = meta.loc[tab.idx, 'pos'].values
# recompute adjusted + subtype-only estimates for each (band, kind)
est = []
for arm, ms in meta.groupby('arm'):
    bex = burden_ex(b, [arm]).to_numpy()
    for adj, Zc in [('M0', np.ones((n,1))), ('M1', np.column_stack([np.ones(n), sub])), ('M3', np.column_stack([np.ones(n), sub, bex, b.purity.to_numpy()]))]:
        Q, _ = np.linalg.qr(Zc); yr = y - Q@(Q.T@y)
        for kind, f in [('loss', lambda M: (M<=-1)), ('gain', lambda M: (M>=1))]:
            cand = ms.index[(f(B[:, ms.index]).sum(0) >= 50) & (f(B[:, ms.index]).sum(0) <= n-50)]
            if len(cand)==0: continue
            Xm = f(B[:, cand]).astype(float); Xr = Xm - Q@(Q.T@Xm); ss = (Xr**2).sum(0); beta = (Xr*yr[:,None]).sum(0)/ss; rss = ((yr[:,None]-Xr*beta)**2).sum(0); se = np.sqrt(rss/(n-Zc.shape[1]-1)/ss)
            for j, c in enumerate(cand): est.append(dict(idx=c, kind=kind, model=adj, beta=beta[j], se=se[j], t=beta[j]/se[j], n_events=int(Xm[:,j].sum())))
est = pd.DataFrame(est).pivot_table(index=['idx','kind'], columns='model', values=['beta','t','se','n_events'], aggfunc='first'); est.columns = [f'{a}_{b_}' for a,b_ in est.columns]; est = est.reset_index()
tab = tab.merge(est, on=['idx','kind'])
for kind in ['loss','gain']:
    i = tab.kind==kind; tab.loc[i,'q_bh'] = multipletests(tab.loc[i,'perm_p'], method='fdr_bh')[1]
    tab.loc[i,'fwer_p'] = [(1+(maxnull[kind] >= abs(tv)).sum())/(P+1) for tv in tab.loc[i,'t_M3']]
tab.to_csv(R+'x5_scan.csv', index=False)
for kind in ['loss','gain']:
    s = tab[tab.kind==kind]; print(kind, 'bands tested', len(s), '| perm_p<0.05:', int((s.perm_p<0.05).sum()), '| FDR<0.1:', int((s.q_bh<0.1).sum()), '| FWER<0.05:', int((s.fwer_p<0.05).sum()))
    print(s.sort_values('perm_p').head(8)[['band','arm','n_events_M3','beta_M1','beta_M3','t_M3','perm_p','q_bh','fwer_p']].round(3).to_string())
s = tab[(tab.kind=='loss')&(tab.band.str.startswith('9p21'))]; print('9p21 bands:'); print(s[['band','n_events_M3','beta_M1','beta_M3','t_M3','perm_p','q_bh','fwer_p']].round(3).to_string())
for band in ['9p21.3']:
    s = tab[(tab.kind=='loss')]; r = (s.perm_p <= s[s.band==band].perm_p.iloc[0]).sum(); print(band, 'rank by perm_p among', len(s), '=', r)
# symmetry: loss vs gain adjusted betas for the same band
wide = tab.pivot_table(index='band', columns='kind', values=['beta_M3','beta_M1','beta_M0'], aggfunc='first').dropna()
wide.columns = [f'{a}_{b_}' for a,b_ in wide.columns]; wide.to_csv(R+'x5_symmetry.csv')
from scipy.stats import spearmanr, pearsonr
for mdl in ['M0','M1','M3']:
    r = spearmanr(wide[f'beta_{mdl}_loss'], wide[f'beta_{mdl}_gain']); print(f'symmetry {mdl}: bands {len(wide)}, spearman loss-vs-gain beta = {r[0]:.3f} (p={r[1]:.2g})')
print('null SD of t: ', 'overall mean |t| observed loss M3:', tab[tab.kind=="loss"].t_M3.abs().mean().round(3))
