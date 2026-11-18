"""Write report/numbers.tex: every number quoted in the text, read from results/*.csv (no hand transcription)."""
import sys, warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'scripts')
import numpy as np, pandas as pd, re
R='results/'; out = {}
def put(k, v, d=2):
    out[k] = (f'{v:,}' if isinstance(v,(int,np.integer)) else f'{v:.{d}f}').replace('-','\\ensuremath{-}')
def ci(k, b, l, h, d=2): put(k+'B', b, d); put(k+'L', l, d); put(k+'H', h, d)
m = pd.read_pickle('data/processed/panc.pkl'); put('nPan', len(m)); put('nCancers', m.ct.nunique())
b = m[(m.ct=='BRCA') & m['sub'].isin(['LumA','LumB','Basal','Her2'])].dropna(subset=['purity','AS','ploidy','WGD','logitLF']); put('nBrca', len(b))
for s in ['LumA','LumB','Basal','Her2']: put('n'+s, int((b['sub']==s).sum()))
put('rhoPurity', b[['LF','purity']].corr(method='spearman').iloc[0,1]); put('rhoFGA', b[['LF','FGA']].corr(method='spearman').iloc[0,1]); put('rhoBurden', b[['LF','burden_ex9p']].corr(method='spearman').iloc[0,1])
sh = pd.read_csv(R+'x1a_shapley.csv', index_col=0); put('shapPurity', sh.loc['Purity','shapley_R2']); put('shapSub', sh.loc['Subtype','shapley_R2']); put('shapAS', sh.loc['AS+ploidy+WGD','shapley_R2']); put('shapBurden', sh.loc['Burden (FGA, 9p-out)','shapley_R2'], 3); put('totR2', sh.total_R2.iloc[0])
a = pd.read_csv(R+'x1b_purity_audit.csv')
for e,t in [('zFGA','FGA'),('zAS','AS')]:
    s = a[a.exposure==e]; neg = s[(s.beta_nopurity<0)&(s.q_nopurity<0.05)]; put('audit'+t+'N', len(s)); put('audit'+t+'Neg', len(neg)); put('audit'+t+'Stay', int(((neg.beta_purity<0)&(neg.q_purity<0.05)).sum())); put('audit'+t+'Rho', s.rho_LF_purity.median())
x2 = pd.read_csv(R+'x2_ladder.csv').set_index(['unit','model'])
for u,k in [('9p21.3 (CDKN2A)','P'),('STAT1','Stat'),('RB1 (control)','Rb'),('B2M','Bt'),('Antigen presentation (B2M/HLA-A/B/C/TAP1/2)','Ag')]:
    for mod,t in [('M0 marginal','0'),('M1 +subtype','1'),('M3 +purity','3')]:
        r = x2.loc[(u,mod)]; ci('lad'+k+t, r.beta, r.lo95, r.hi95, 2)
    put('lad'+k+'N', int(x2.loc[(u,'M3 +purity')].n_loss))
rv = pd.read_csv(R+'x2b_rivals.csv'); base = rv.iloc[0]
for i,k in [(1,'Cdk'),(2,'Arm'),(3,'Both'),(4,'Ifn')]: put('dBIC'+k, rv.BIC.iloc[i]-base.BIC, 1); put('cv'+k, rv.cvR2.iloc[i], 3)
put('cvBase', base.cvR2, 3)
dr = pd.read_csv(R+'x3_dose_cdkn2a.csv', index_col=0)
for key,k in [('[T.-2]','Deep'),('[T.-1]','Shal'),('[T.1]','Gain')]:
    r = dr[dr.index.str.contains(key, regex=False)].iloc[0]; ci('dose'+k, r.coef, r.lo, r.hi, 2); put('dose'+k+'N', int(r.n))
gl = pd.read_csv(R+'x3_gain_loss.csv')
for u,k in [('9p21.3 (CDKN2A)','P'),('PTEN (control)','Pten')]:
    r = gl[(gl.unit==u)&(gl.kind=='gain')].iloc[0]; ci('gain'+k, r.beta, r.lo95, r.hi95, 2)
fa = pd.read_csv(R+'x6_focal_arm.csv'); r = fa[(fa.pair=='9p21.3 / 9p')&(fa.model=='gene loss | arm loss')].iloc[0]; ci('brcaGeneGivenArm', r.beta, r.lo95, r.hi95); put('brcaGeneGivenArmP', r.p, 3); put('nFocal', int(r.n_focal)); put('nArm', int(r.n_arm)); put('nGene', int(r.n_gene))
a9 = pd.read_csv(R+'x9_aipw_brca.csv', index_col=0).iloc[:,0]
put('aipw', float(a9['aipw']), 3); put('aipwL', float(a9['lo95']), 3); put('aipwH', float(a9['hi95']), 3); put('aipwLn', float(a9['lo90']), 3); put('aipwHn', float(a9['hi90']), 3); put('att', float(a9['matched_att']), 3); put('attSE', float(a9['matched_se']), 3); put('attN', int(float(a9['matched_n']))); put('evalue', float(a9['evalue_est']), 2); put('residSD', float(a9['resid_sd']), 2)
x7 = pd.read_csv(R+'x7_subtype.csv')
for s in ['LumA','LumB','Basal','Her2']:
    r = x7[(x7.unit=='9p21.3 (CDKN2A)')&(x7.subtype==s)].iloc[0]; ci('sub'+s, r.beta, r.lo95, r.hi95); put('sub'+s+'P', r.p, 3)
put('intP4', x7[(x7.unit=='9p21.3 (CDKN2A)')&(x7.subtype.str.startswith('interaction x'))].p.iloc[0], 3); put('intPlb', x7[(x7.unit=='9p21.3 (CDKN2A)')&(x7.subtype.str.startswith('interaction luminal'))].p.iloc[0], 3)
x8 = pd.read_csv(R+'x8_percancer_9p.csv').set_index('ct'); put('nModelledCancers', int(x8.b3.notna().sum()))
for c in ['HNSC','SKCM','LUAD','LUSC','BRCA','KIRC','LGG']:
    r = x8.loc[c]; ci('pc'+c, r.b3, r.lo3, r.hi3); put('pc'+c+'P', r.p3, 3); put('pc'+c+'Q', r.q3, 3); put('pc'+c+'N', int(r.n)); put('pc'+c+'Loss', int(r.n_loss))
me = pd.read_csv(R+'x8_meta.csv', index_col=0)
for key,k in [('all cancers','All'),('positive controls','Pos'),('all except positive controls & BRCA','Oth')]:
    r = me.loc[key]; ci('meta'+k, r.beta, r.lo, r.hi, 3); put('meta'+k+'K', int(r.k)); put('meta'+k+'I', r.I2*100, 0)
c1 = pd.read_csv(R+'x4_class_counts_thr1.csv', index_col=0).sum(); c2 = pd.read_csv(R+'x4_class_counts_thr2.csv', index_col=0).sum()
put('clsOneShal', int(c1['C1 IFN+MTAP+CDKN2A lost'])); put('clsTwoShal', int(c1['C2 CDKN2A lost, IFN intact'])); put('clsThreeShal', int(c1['C3 IFN lost, CDKN2A intact'])); put('clsNoneShal', int(c1['C4 neither']))
put('clsOneDeep', int(c2['C1 IFN+MTAP+CDKN2A lost'])); put('clsTwoDeep', int(c2['C2 CDKN2A lost, IFN intact'])); put('clsThreeDeep', int(c2['C3 IFN lost, CDKN2A intact']))
cp = pd.read_csv(R+'x4b_classes_pos.csv'); cp = cp[cp.cls.notna()]
for thr,t in [(-1,'S'),(-2,'D')]:
    for c in ['C1','C2']:
        for st,k in [('pooled positive controls','Pos'),('pooled all other cancers','Oth')]:
            r = cp[(cp.thr==thr)&(cp.set==st)&(cp.cls==c)].iloc[0]; ci(f'cls{c}{k}{t}', r.beta, r.lo, r.hi, 3); put(f'cls{c}{k}{t}N', int(r.n_cls))
bp = pd.read_csv(R+'x6b_arm_pan.csv')
for st,k in [('pooled positive controls','Pos'),('pooled other cancers (no BRCA)','Oth')]:
    for mod,t in [('gene | arm','GA'),('arm | gene','AG')]:
        r = bp[(bp.set==st)&(bp.model==mod)].iloc[0]; ci('arm'+k+t, r.beta, r.lo, r.hi, 3)
sc = pd.read_csv(R+'x5_scan.csv'); L = sc[sc.kind=='loss']; G = sc[sc.kind=='gain']
put('scanLossN', len(L)); put('scanLossP05', int((L.perm_p<0.05).sum())); put('scanLossFDR', int((L.q_bh<0.1).sum())); put('scanLossFWER', int((L.fwer_p<0.05).sum()))
put('scanGainN', len(G)); put('scanGainFDR', int((G.q_bh<0.1).sum())); put('scanGainFWER', int((G.fwer_p<0.05).sum()))
r = L[L.band=='9p21.3'].iloc[0]; put('scan9pP', r.perm_p, 3); put('scan9pT', r.t_M3, 2); put('scan9pRank', int((L.perm_p <= r.perm_p).sum())); put('scan9pBeta', r.beta_M3, 3)
g1 = G[G.q_bh<0.1]; out['gainTopArm'] = str(g1.arm.value_counts().index[0]) if len(g1) else '--'; put('gainTopArmN', int((g1.arm==g1.arm.value_counts().index[0]).sum()) if len(g1) else 0)
jb = pd.read_csv(R+'x5b_joint_loss_gain.csv'); from scipy.stats import spearmanr
put('jointBands', len(jb)); put('jointRho', spearmanr(jb.bL, jb.bG)[0], 3)
x10 = pd.read_csv(R+'x10_rna_outcomes.csv').set_index('outcome')
for o,k in [('cis_CDKN2A','CisC'),('cis_MTAP','CisM'),('CYT','Cyt'),('CD8','Cd8'),('IFNG_sig','Ifng'),('PTPRC','Ptprc'),('ISG','Isg'),('MHC_I','Mhc'),('CXCL9_10','Cxcl'),('lf_share','Share')]:
    r = x10.loc[o]; ci('rna'+k, r.M3_beta, r.M3_lo, r.M3_hi); put('rna'+k+'P', r.M3_p, 3)
x11 = pd.read_csv(R+'x11_ccle.csv')
for st,k in [('Breast lines','Br'),('All lineages (lineage FE)','All')]:
    for o,t in [('ISG','Isg'),('MHC_I','Mhc')]:
        r = x11[(x11.set==st)&(x11.outcome==o)&(x11.exposure=='x')].iloc[0]; ci('ccle'+k+t, r.beta, r.lo, r.hi); put('ccle'+k+t+'P', r.p, 3); put('ccle'+k+'N', int(r.n)); put('ccle'+k+'L', int(r.n_loss))
for cl,k in [('C1 CDKN2A+IFN lost','One'),('C2 CDKN2A lost, IFN intact','Two')]:
    r = x11[(x11.exposure==cl)&(x11.outcome=='ISG')].iloc[0]; ci('ccleCls'+k, r.beta, r.lo, r.hi); put('ccleCls'+k+'P', r.p, 3); put('ccleCls'+k+'N', int(r.n_loss))
x12 = pd.read_csv(R+'x12_metabric.csv'); put('mbN', int(x12.n.iloc[0])); put('mbLoss', int(x12.n_loss.iloc[0]))
for o,k in [('CYT','Cyt'),('CD8','Cd8'),('IFNG_sig','Ifng'),('PTPRC','Ptprc'),('ISG','Isg'),('MHC_I','Mhc')]:
    r = x12[(x12.outcome==o)&(x12.exposure=='x')].iloc[0]; ci('mb'+k, r.M3_beta, r.M3_lo, r.M3_hi); put('mb'+k+'P', r.M3_p, 3)
mbs = pd.read_csv(R+'x12_metabric_subtype.csv').set_index('subtype'); ci('mbBasal', mbs.loc['Basal','beta'], mbs.loc['Basal','lo'], mbs.loc['Basal','hi'])
rb = pd.read_csv(R+'x14_robust.csv').set_index('spec'); put('nRobust', len(rb)); put('robEquiv', int(rb['equiv_0.10'].sum()))
lo_, hi_ = rb.lo95.min(), rb.hi95.max(); put('robLoMin', rb.beta.min(), 3); put('robBetaMax', rb.beta.max(), 3)
r = rb.loc['Deep deletion only (<= -2)']; ci('robDeep', r.beta, r.lo95, r.hi95, 3)
sy = pd.read_csv(R+'x5_symmetry.csv'); put('sepRho', spearmanr(sy.beta_M3_loss, sy.beta_M3_gain)[0], 2)
for mod,t in [('M2 +burden(arm-out)','2')]:
    r = x2.loc[('9p21.3 (CDKN2A)',mod)]; ci('ladP'+t, r.beta, r.lo95, r.hi95, 2)
put('ladPThreeAbsL', abs(x2.loc[('9p21.3 (CDKN2A)','M3 +purity')].lo95), 2)
W = ['Zero','One','Two','Three','Four','Five','Six','Seven','Eight','Nine']
name = lambda k: ''.join(W[int(c)] if c.isdigit() else c for c in k)
open('report/numbers.tex','w').write('% auto-generated by scripts/make_numbers.py\n' + ''.join('\\newcommand{\\n%s}{%s}\n' % (name(k), v) for k,v in out.items()))
print(len(out), 'macros')
