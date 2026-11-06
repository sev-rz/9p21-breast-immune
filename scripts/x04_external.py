"""X10 (RNA outcomes), H4 mechanism, X11 (CCLE tumour-intrinsic), X12 (METABRIC), via cBioPortal API."""
import sys, json, time, warnings, urllib.request; warnings.filterwarnings('ignore'); sys.path.insert(0,'scripts')
from common import *
import statsmodels.formula.api as smf
R = 'results/'; B = 'https://www.cbioportal.org/api'
def req(path, body=None, tries=4):
    for i in range(tries):
        try:
            r = urllib.request.Request(B+path, json.dumps(body).encode() if body is not None else None, {'Content-Type':'application/json','Accept':'application/json'})
            return json.load(urllib.request.urlopen(r, timeout=300))
        except Exception as e:
            if i == tries-1: raise
            time.sleep(3)
def entrez(symbols):
    out = req('/genes/fetch?geneIdType=HUGO_GENE_SYMBOL', symbols); return {g['hugoGeneSymbol']: g['entrezGeneId'] for g in out}
def mol(profile, samplelist, ent, chunk=150):
    rows = []; ids = list(ent.values())
    for i in range(0, len(ids), chunk):
        rows += req(f'/molecular-profiles/{profile}/molecular-data/fetch?projection=SUMMARY', {'sampleListId': samplelist, 'entrezGeneIds': ids[i:i+chunk]})
    d = pd.DataFrame(rows); inv = {v:k for k,v in ent.items()}; d['g'] = d.entrezGeneId.map(inv)
    return d.pivot_table(index='sampleId', columns='g', values='value', aggfunc='first')
CYT = ['GZMA','PRF1']; CD8 = ['CD8A','CD8B']; AYERS = ['CCL5','CD27','CD274','CD276','CD8A','CMKLR1','CXCL9','CXCR6','HLA-DQA1','HLA-DRB1','HLA-E','IDO1','LAG3','NKG7','PDCD1LG2','PSMB10','STAT1','TIGIT']
ISG = ['IFI6','IFI27','IFIT1','IFIT3','ISG15','MX1','OAS1','OAS3','IRF7','IFI44','RSAD2','IFI35']; MHC = ['HLA-A','HLA-B','HLA-C','B2M','TAP1','TAP2','PSMB9']
CIS = ['CDKN2A','MTAP','IFNE','B2M','STAT1']; EXTRA = ['PTPRC','CXCL9','CXCL10']
allg = sorted(set(CYT+CD8+AYERS+ISG+MHC+CIS+EXTRA)); ent = entrez(allg); print('genes resolved', len(ent), 'of', len(allg), [g for g in allg if g not in ent])
z = lambda s: (s-s.mean())/s.std()
def scores(E, log=False):
    L = np.log2(E.astype(float)+1) if log else E.astype(float)
    S = pd.DataFrame({'CYT': L[CYT].mean(1), 'CD8': L[CD8].mean(1), 'IFNG_sig': L[[g for g in AYERS if g in L]].mean(1), 'ISG': L[[g for g in ISG if g in L]].mean(1), 'MHC_I': L[[g for g in MHC if g in L]].mean(1), 'CXCL9_10': L[['CXCL9','CXCL10']].mean(1), 'PTPRC': L['PTPRC']})
    for g in ['CDKN2A','MTAP','B2M','STAT1']: S['cis_'+g] = L[g]
    return S.apply(z)
# ================= TCGA-BRCA RNA (X10 + H4)
m = pd.read_pickle(P+'panc.pkl'); b = m[(m.ct=='BRCA') & m['sub'].isin(['LumA','LumB','Basal','Her2'])].dropna(subset=['purity','logitLF']).copy()
E = mol('brca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna', 'brca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna', ent); S = scores(E, log=True); S.index = [i[:15] for i in S.index]
b = b.merge(S, left_on='key', right_index=True, how='inner'); print('TCGA-BRCA with RNA:', len(b))
b['bex'] = burden_ex(b, ['9p']); b['x'] = (b.g_CDKN2A <= -1).astype(int); b['lf_share'] = LOGIT((b.LF/(1-b.purity)).clip(0.01, 0.99)); b['z_logitLF'] = z(b.logitLF)
rows = []
for out in ['z_logitLF','lf_share','CYT','CD8','IFNG_sig','PTPRC','ISG','MHC_I','CXCL9_10','cis_CDKN2A','cis_MTAP','cis_B2M','cis_STAT1']:
    f1 = fit_ols(b, out, 'x', ['C(sub)']); f3 = fit_ols(b, out, 'x', ['C(sub)','bex','purity']); rows.append(dict(outcome=out, M1_beta=f1['beta'], M1_lo=f1['lo95'], M1_hi=f1['hi95'], M3_beta=f3['beta'], M3_lo=f3['lo95'], M3_hi=f3['hi95'], M3_p=f3['p'], n=f3['n']))
x10 = pd.DataFrame(rows); x10.to_csv(R+'x10_rna_outcomes.csv', index=False); print(x10.round(3).to_string())
# other loci for RNA outcomes: STAT1/B2M loss on own expression
for gn, arm, ce in [('STAT1','2q','cis_STAT1'),('B2M','15q','cis_B2M'),('MTAP','9p','cis_MTAP')]:
    b['bx'] = burden_ex(b, [arm]); b['xx'] = (b['g_'+gn] <= -1).astype(int); f = fit_ols(b, ce, 'xx', ['C(sub)','bx','purity']); print('cis', gn, round(f['beta'],3), round(f['lo95'],3), round(f['hi95'],3))
# dose of copy number on own expression
b['callC'] = b.g_CDKN2A.clip(-2,1).astype(int).astype(str); r = smf.ols('cis_CDKN2A ~ C(callC) + C(sub) + bex + purity', b).fit(cov_type='HC3'); print('CDKN2A expression by call:\n', r.params.filter(like='callC').round(3))
# ================= CCLE (X11)
cl = pd.DataFrame(req('/studies/ccle_broad_2019/clinical-data?clinicalDataType=SAMPLE&projection=SUMMARY&pageSize=100000')); cw = cl.pivot_table(index='sampleId', columns='clinicalAttributeId', values='value', aggfunc='first')
ge = sorted(set(ISG+MHC+CIS+['CDKN2B','IFNB1','IFNA1','CYT'])-{'CYT'}); ent2 = entrez(ge)
Ec = mol('ccle_broad_2019_rna_seq_mrna', 'ccle_broad_2019_all', ent2); Cc = mol('ccle_broad_2019_cna', 'ccle_broad_2019_all', {g:ent2[g] for g in ['CDKN2A','CDKN2B','MTAP','IFNB1','IFNA1']})
Sc = pd.DataFrame({'ISG': Ec[[g for g in ISG if g in Ec]].mean(1), 'MHC_I': Ec[[g for g in MHC if g in Ec]].mean(1), 'cis_CDKN2A': Ec['CDKN2A'], 'cis_MTAP': Ec['MTAP']}).astype(float)
cc = cw.join(Sc, how='inner').join(Cc.add_prefix('cna_'), how='inner'); cc['FGA'] = pd.to_numeric(cc['FRACTION_GENOME_ALTERED'], errors='coerce'); cc['PUR'] = pd.to_numeric(cc.get('PURITY'), errors='coerce')
cc['x'] = (pd.to_numeric(cc.cna_CDKN2A, errors='coerce') <= -1).astype(int); cc['xd'] = (pd.to_numeric(cc.cna_CDKN2A, errors='coerce') <= -2).astype(int); cc['ifn'] = (pd.to_numeric(cc.cna_IFNB1, errors='coerce') <= -1).astype(int)
cc = cc.dropna(subset=['FGA','ISG','MHC_I']); print('CCLE lines:', len(cc), '| breast:', int((cc.CANCER_TYPE=='Breast Cancer').sum()))
rows = []
for lab, sel in [('Breast lines', cc.CANCER_TYPE=='Breast Cancer'), ('All lineages (lineage FE)', cc.CANCER_TYPE.notna())]:
    z_ = cc[sel].copy(); z_['lin'] = z_.CANCER_TYPE
    for out in ['ISG','MHC_I','cis_CDKN2A','cis_MTAP']:
        z_['y'] = (z_[out]-z_[out].mean())/z_[out].std()
        for xx in ['x','xd']:
            if z_[xx].sum() < 8: continue
            cov = ['FGA'] + (['C(lin)'] if 'lineages' in lab else [])
            f = fit_ols(z_, 'y', xx, cov); rows.append(dict(set=lab, outcome=out, exposure=xx, n=f['n'], n_loss=int(z_[xx].sum()), beta=f['beta'], lo=f['lo95'], hi=f['hi95'], p=f['p']))
# R1 vs R2 in cell lines: CDKN2A lost with IFNB1 lost vs CDKN2A lost with IFNB1 intact
cc['cls'] = np.select([(cc.x==1)&(cc.ifn==1), (cc.x==1)&(cc.ifn==0)], ['C1 CDKN2A+IFN lost','C2 CDKN2A lost, IFN intact'], 'C4 neither')
print('CCLE class counts (all):', cc.cls.value_counts().to_dict(), '| breast:', cc[cc.CANCER_TYPE=='Breast Cancer'].cls.value_counts().to_dict())
for out in ['ISG','MHC_I']:
    cc['y'] = (cc[out]-cc[out].mean())/cc[out].std(); f = smf.ols('y ~ C(cls, Treatment("C4 neither")) + FGA + C(CANCER_TYPE)', cc).fit(cov_type='HC3')
    for k in f.params.index:
        if 'cls' in k: rows.append(dict(set='All lineages (lineage FE)', outcome=out, exposure=k.split('[T.')[1][:-1], n=int(f.nobs), n_loss=int((cc.cls==k.split('[T.')[1][:-1]).sum()), beta=f.params[k], lo=f.conf_int().loc[k,0], hi=f.conf_int().loc[k,1], p=f.pvalues[k]))
x11 = pd.DataFrame(rows); x11.to_csv(R+'x11_ccle.csv', index=False); print(x11.round(3).to_string())
# ================= METABRIC (X12)
g_all = pd.read_csv(BIG+'thr_genes.csv'); g_all['arm'] = g_all.Cytoband.str.extract(r'^(\d+|X|Y)([pq])').apply(lambda r: ''.join(r.dropna()) if r.notna().all() else np.nan, axis=1)
samp = g_all.dropna(subset=['arm']).iloc[::40]; ent_b = {r['Gene Symbol']: int(r['Locus ID']) for _, r in samp.iterrows() if str(r['Locus ID']).isdigit()}; arm_of = dict(zip(samp['Gene Symbol'], samp.arm))
for gnm in ['CDKN2A','MTAP','IFNB1','IFNE','B2M','STAT1']: ent_b[gnm] = entrez([gnm])[gnm]
print('METABRIC CNA genes requested', len(ent_b))
Cm = mol('brca_metabric_cna', 'brca_metabric_all', ent_b, chunk=100).apply(pd.to_numeric, errors='coerce')
Em = mol('brca_metabric_mrna', 'brca_metabric_all', ent).apply(pd.to_numeric, errors='coerce'); Sm = scores(Em, log=False)
cm = pd.DataFrame(req('/studies/brca_metabric/clinical-data?clinicalDataType=SAMPLE&projection=SUMMARY&pageSize=100000')); cmw = cm.pivot_table(index='sampleId', columns='clinicalAttributeId', values='value', aggfunc='first')
cp = pd.DataFrame(req('/studies/brca_metabric/clinical-data?clinicalDataType=PATIENT&projection=SUMMARY&pageSize=100000')); cpw = cp.pivot_table(index='patientId', columns='clinicalAttributeId', values='value', aggfunc='first')
for c_ in cpw.columns:
    if c_ not in cmw.columns: cmw[c_] = cmw.index.map(cpw[c_])
mb = cmw.join(Sm, how='inner').join(Cm, how='inner'); mb['sub'] = mb.CLAUDIN_SUBTYPE; mb = mb[mb['sub'].isin(['LumA','LumB','Her2','Basal'])].copy()
cols = [c for c in Cm.columns if c in arm_of]; alt = (Cm[cols].abs() >= 1); armv = pd.Series({c: arm_of[c] for c in cols})
mb['burden_all'] = alt.loc[mb.index].mean(1); mb['burden_ex9p'] = alt.loc[mb.index, armv[armv!='9p'].index].mean(1)
mb['x'] = (mb.CDKN2A <= -1).astype(int); mb['xifn'] = (mb.IFNB1 <= -1).astype(int); mb['CELL'] = mb.CELLULARITY.fillna('NA')
print('METABRIC n =', len(mb), mb['sub'].value_counts().to_dict(), '| CDKN2A loss freq', round(mb.x.mean(), 3), '| cellularity', mb.CELL.value_counts().to_dict())
rows = []
for out in ['CYT','CD8','IFNG_sig','PTPRC','ISG','MHC_I','cis_CDKN2A','cis_MTAP']:
    for xx in ['x','xifn']:
        f1 = fit_ols(mb, out, xx, ['C(sub)']); f3 = fit_ols(mb, out, xx, ['C(sub)','burden_ex9p','C(CELL)'])
        rows.append(dict(outcome=out, exposure=xx, n=f3['n'], n_loss=int(mb[xx].sum()), M1_beta=f1['beta'], M1_lo=f1['lo95'], M1_hi=f1['hi95'], M3_beta=f3['beta'], M3_lo=f3['lo95'], M3_hi=f3['hi95'], M3_p=f3['p']))
x12 = pd.DataFrame(rows); x12.to_csv(R+'x12_metabric.csv', index=False); print(x12.round(3).to_string())
# subtype split in METABRIC
rows = []
for s in ['LumA','LumB','Her2','Basal']:
    zz = mb[mb['sub']==s]; f = fit_ols(zz, 'CYT', 'x', ['burden_ex9p','C(CELL)']); rows.append(dict(subtype=s, n=f['n'], n_loss=int(zz.x.sum()), beta=f['beta'], lo=f['lo95'], hi=f['hi95'], p=f['p']))
pd.DataFrame(rows).to_csv(R+'x12_metabric_subtype.csv', index=False); print(pd.DataFrame(rows).round(3).to_string())
