"""Shared loaders. Raw downloads are read from data/raw/, files built by the scripts go to data/processed/."""
import os, re, numpy as np, pandas as pd
D = os.environ.get('DATA_DIR', 'data/raw/'); P = os.environ.get('PROC_DIR', 'data/processed/'); BIG = os.environ.get('BIG_DIR', P)
LOGIT = lambda p: np.log(p/(1-p))

def load_master():
    """One row per primary tumour (sample type 01): LF, purity, ploidy, WGD, global burden, per-arm burden."""
    f = P + 'master.pkl'
    if os.path.exists(f): return pd.read_pickle(f)
    l = pd.read_csv(D+'leuk.tsv', sep='\t', header=None, names=['ct','bc','LF'])
    l['key'] = l.bc.str[:15]; l = l[l.key.str.endswith('-01') | ((l.ct=='SKCM') & l.key.str.endswith('-06'))].drop_duplicates('key')  # SKCM primaries are rare; metastatic (-06) allowed for SKCM only (deviation from prereg, documented)
    m = pd.read_csv(D+'mastercalls.txt', sep='\t'); m['key'] = m['array'].str[:15]
    m = m.drop_duplicates('key')[['key','purity','ploidy','Genome doublings']].rename(columns={'Genome doublings':'WGD'})
    sg = pd.read_csv(D+'seg_scores.tsv', sep='\t'); sg['key'] = sg.Sample.str[:15]; sg = sg.drop_duplicates('key')[['key','frac_altered','n_segs']].rename(columns={'frac_altered':'FGA'})
    ab = pd.read_csv(D+'absolute_scores.tsv', sep='\t', index_col=0); ab['key'] = ab.index.str[:15]; ab = ab.drop_duplicates('key')[['key','AS','LOH_frac_altered']]
    d = l[['key','ct','LF']].merge(m, on='key', how='left').merge(sg, on='key', how='left').merge(ab, on='key', how='left')
    d['logitLF'] = LOGIT(d.LF.clip(1e-3, 1-1e-3)); d['pat'] = d.key.str[:12]
    d.to_pickle(f); return d

def load_gistic():
    """int8 genes x samples thresholded calls (rows in genomic order), gene meta with arm, sample keys (first occurrence kept)."""
    X = np.load(BIG+'thr.npy', mmap_mode='r'); g = pd.read_csv(BIG+'thr_genes.csv')
    s = pd.read_csv(BIG+'thr_samples.csv').barcode.values; keys = np.array([b[:15] for b in s])
    _, first = np.unique(keys, return_index=True); first.sort()
    keep = first[[k.endswith('-01') or k.endswith('-06') for k in keys[first]]]
    g['arm'] = g.Cytoband.str.extract(r'^(\d+|X|Y)([pq])')[0].astype(str) + g.Cytoband.str.extract(r'^(?:\d+|X|Y)([pq])')[0].astype(str)
    g.loc[g.Cytoband.isna(), 'arm'] = np.nan
    return np.asarray(X[:, keep]), g, keys[keep]

def burden(X, g, keys):
    """Gene-count burden: global fraction altered, and per-arm altered counts / gene counts for leave-arm-out."""
    alt = (np.abs(X) >= 1)
    arms = g.arm.values; ok = pd.notna(arms)
    tot = alt[ok].sum(0).astype(float); ntot = ok.sum()
    per = {}; n = {}
    for a in pd.unique(arms[ok]):
        idx = np.where(arms == a)[0]; per[a] = alt[idx].sum(0).astype(float); n[a] = len(idx)
    return tot/ntot, per, n, ntot, tot

def fit_ols(df, y, x, covs, robust='HC3'):
    import statsmodels.formula.api as smf
    f = f'{y} ~ {x}' + ''.join(f' + {c}' for c in covs)
    r = smf.ols(f, df).fit(cov_type=robust)
    b = r.params[x]; se = r.bse[x]; ci = r.conf_int(0.05).loc[x].values; ci90 = r.conf_int(0.10).loc[x].values
    return dict(beta=b, se=se, lo95=ci[0], hi95=ci[1], lo90=ci90[0], hi90=ci90[1], p=r.pvalues[x], n=int(r.nobs), r2=r.rsquared, aic=r.aic)

def classify(m1, mf, sesoi=0.10):
    """Pre-registered outcome class from subtype-only (m1) and fully adjusted (mf) fits."""
    eq = (mf['lo90'] > -sesoi) and (mf['hi90'] < sesoi)
    base_sig = (m1['lo95'] > 0) or (m1['hi95'] < 0)
    adj_sig = (mf['lo95'] > 0) or (mf['hi95'] < 0)
    if base_sig and adj_sig and np.sign(mf['beta']) != np.sign(m1['beta']): return 'Suppressed (sign flip)'
    if adj_sig and base_sig and abs(mf['beta']) >= 0.5*abs(m1['beta']): return 'Locus-specific'
    if adj_sig and base_sig: return 'Attenuated, still significant'
    if base_sig and not adj_sig and abs(mf['beta']) < 0.5*abs(m1['beta']): return 'Burden/purity-driven'
    if eq: return 'Equivalent to null'
    if not adj_sig and not base_sig: return 'Underpowered (CI beyond SESOI)'
    if adj_sig and not base_sig: return 'Emerges only after adjustment'
    return 'Mixed/other'

import json
def burden_ex(df, arms):
    """Gene-count burden excluding the given arms (leave-arm-out)."""
    meta = json.load(open(P+'arm_meta.json')); arms = [a for a in dict.fromkeys(arms)]
    alt = df['alt_total'].copy(); nn = meta['ntot']
    for a in arms: alt = alt - df['af_'+a]*meta['n'][a]; nn -= meta['n'][a]
    return alt/nn
