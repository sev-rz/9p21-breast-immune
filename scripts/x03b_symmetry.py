"""Joint loss+gain model per band (reference = neutral), to avoid reference-category contamination in the loss-vs-gain comparison."""
import sys, warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'scripts')
from common import *
from scipy.stats import spearmanr
X, g, keys = load_gistic(); m = pd.read_pickle(P+'panc.pkl')
b = m[(m.ct=='BRCA') & m['sub'].isin(['LumA','LumB','Basal','Her2'])].dropna(subset=['purity','logitLF']).reset_index(drop=True)
kidx = {k:i for i,k in enumerate(keys)}; Xb = X[:, [kidx[k] for k in b.key]]; g = g.reset_index(drop=True); ok = g.Cytoband.notna().values & g.arm.notna().values
n = len(b); y = b.logitLF.to_numpy(); sub = pd.get_dummies(b['sub'], drop_first=True).astype(float).to_numpy(); rows = []
for cb, ix in g[ok].groupby('Cytoband', sort=False).groups.items():
    ix = np.array(list(ix)); med = np.median(Xb[ix], axis=0); L = (med<=-1).astype(float); G = (med>=1).astype(float)
    if L.sum() < 40 or G.sum() < 40: continue
    arm = g.arm[ix[0]]; Z = np.column_stack([np.ones(n), sub, burden_ex(b,[arm]).to_numpy(), b.purity.to_numpy(), L, G])
    beta, *_ = np.linalg.lstsq(Z, y, rcond=None); res = y - Z@beta; s2 = res@res/(n-Z.shape[1]); cov = s2*np.linalg.inv(Z.T@Z)
    rows.append(dict(band=cb, arm=arm, nL=int(L.sum()), nG=int(G.sum()), bL=beta[-2], seL=np.sqrt(cov[-2,-2]), bG=beta[-1], seG=np.sqrt(cov[-1,-1])))
d = pd.DataFrame(rows); d.to_csv('results/x5b_joint_loss_gain.csv', index=False)
r = spearmanr(d.bL, d.bG); print('bands', len(d), 'spearman(loss,gain) joint model =', round(r[0],3), 'p=', r[1])
print('same sign (both neg or both pos):', int((np.sign(d.bL)==np.sign(d.bG)).sum()), 'opposite:', int((np.sign(d.bL)!=np.sign(d.bG)).sum()))
print('mean beta loss', d.bL.mean().round(3), 'mean beta gain', d.bG.mean().round(3))
b9 = d[d.band=='9p21.3']; print(b9.round(3).to_string())
