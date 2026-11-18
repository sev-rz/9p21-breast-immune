import sys, warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'scripts')
import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import *
R = 'results/'; F = 'report/figures/'
plt.rcParams.update({'font.size': 8, 'axes.titlesize': 8.5, 'axes.labelsize': 8, 'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5, 'legend.fontsize': 7.5, 'font.family': 'DejaVu Sans', 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42, 'axes.linewidth': 0.7})
OI = {'blue':'#0072B2','orange':'#E69F00','green':'#009E73','verm':'#D55E00','sky':'#56B4E9','pink':'#CC79A7','yellow':'#F0E442','grey':'#7f7f7f'}
SUBC = {'LumA':OI['blue'],'LumB':OI['sky'],'Basal':OI['verm'],'Her2':OI['green']}
def panel(ax, s): ax.text(-0.16, 1.14, s, transform=ax.transAxes, fontsize=10, fontweight='bold')
def save(fig, name): fig.savefig(F+name+'.pdf', bbox_inches='tight'); fig.savefig(F+name+'.png', dpi=160, bbox_inches='tight'); plt.close(fig)
def forest(ax, df, lab, b, lo, hi, color=None, sesoi=True, xlab='Effect on logit(LF)', ypos=None, markers=None):
    y = np.arange(len(df))[::-1] if ypos is None else ypos
    if sesoi: ax.axvspan(-0.10, 0.10, color='#e5e5e5', zorder=0, lw=0)
    ax.axvline(0, color='k', lw=0.6, zorder=1)
    cols = color if color is not None else [OI['blue']]*len(df)
    for i,(yy,(_,r)) in enumerate(zip(y, df.iterrows())):
        ax.plot([r[lo], r[hi]], [yy,yy], color=cols[i], lw=1.4, zorder=2); ax.plot(r[b], yy, 'o', color=cols[i], ms=4, zorder=3)
    ax.set_yticks(y); ax.set_yticklabels(df[lab]); ax.set_xlabel(xlab); ax.set_ylim(min(y)-0.7, max(y)+0.7)
m = pd.read_pickle(P+'panc.pkl'); brca = m[(m.ct=='BRCA') & m['sub'].isin(['LumA','LumB','Basal','Her2'])].dropna(subset=['purity','AS','ploidy','WGD','logitLF'])
# ---------------- Fig 1 confounding structure
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.5), gridspec_kw={'width_ratios':[1.15,0.9,1.0]})
for s,c in SUBC.items():
    z = brca[brca['sub']==s]; ax[0].scatter(z.purity, z.LF, s=4, color=c, alpha=0.55, lw=0, label=s)
rho = brca[['purity','LF']].corr(method='spearman').iloc[0,1]; rho2 = brca[['burden_ex9p','LF']].corr(method='spearman').iloc[0,1]
ax[0].set_xlabel('ABSOLUTE tumour purity'); ax[0].set_ylabel('Leukocyte fraction (LF)'); ax[0].legend(frameon=False, markerscale=2, loc='upper right', handlelength=0.5); ax[0].text(0.03,0.04,f'Spearman rho = {rho:.2f}', transform=ax[0].transAxes); panel(ax[0],'a')
sh = pd.read_csv(R+'x1a_shapley.csv', index_col=0); sh = sh.sort_values('shapley_R2'); ax[1].barh(sh.index, sh.shapley_R2, color=[OI['verm'] if 'Purity' in i else OI['grey'] for i in sh.index]); ax[1].set_xlabel('Shapley share of $R^2$ (BRCA)'); ax[1].text(0.97,0.06,f"total $R^2$ = {sh.total_R2.iloc[0]:.2f}", transform=ax[1].transAxes, ha='right'); panel(ax[1],'b')
a = pd.read_csv(R+'x1b_purity_audit.csv'); a = a[a.exposure=='zFGA']
ax[2].plot([-0.6,0.2],[-0.6,0.2], color='k', lw=0.6, ls='--'); sig = a.q_purity<0.05
ax[2].scatter(a.beta_nopurity[sig], a.beta_purity[sig], s=14, color=OI['blue'], label='FDR<0.05 after purity'); ax[2].scatter(a.beta_nopurity[~sig], a.beta_purity[~sig], s=14, color='white', edgecolor=OI['blue'], label='not significant')
for _,r in a[a.ct.isin(['BRCA','HNSC'])].iterrows(): ax[2].annotate(r.ct, (r.beta_nopurity, r.beta_purity), fontsize=6.5, xytext=(3,2), textcoords='offset points')
ax[2].set_xlabel('FGA effect on logit(LF), no purity'); ax[2].set_ylabel('same, + purity'); ax[2].legend(frameon=False, loc='upper left', fontsize=6.5); panel(ax[2],'c')
fig.tight_layout(w_pad=1.2); save(fig, 'fig1_confounding')
# ---------------- Fig 2 attenuation forest
x2 = pd.read_csv(R+'x2_ladder.csv'); units = ['9p21.3 (CDKN2A)','MTAP','IFN cluster (any)','Antigen presentation (B2M/HLA-A/B/C/TAP1/2)','IFN signalling (JAK1/2,STAT1,IRF1,IFNGR1/2)','B2M','STAT1','JAK1','HLA-A','TAP1','JAK2','PTEN (control)','RB1 (control)']
short = {'Antigen presentation (B2M/HLA-A/B/C/TAP1/2)':'Antigen presentation (any of 6)','IFN signalling (JAK1/2,STAT1,IRF1,IFNGR1/2)':'IFN signalling (any of 6)','9p21.3 (CDKN2A)':'9p21.3 (CDKN2A)'}
fig, axs = plt.subplots(1, 3, figsize=(7.2, 3.6), sharey=True)
for ax, (mod, t, c) in zip(axs, [('M0 marginal','Marginal',OI['grey']),('M1 +subtype','+ subtype',OI['sky']),('M3 +purity','+ subtype, burden, purity',OI['blue'])]):
    d = x2[(x2.model==mod)].set_index('unit').loc[units].reset_index(); d['lab'] = d.unit.map(lambda u: short.get(u,u))
    forest(ax, d, 'lab', 'beta', 'lo95', 'hi95', color=[c]*len(d), xlab='logit(LF) effect of loss'); ax.set_title(t); ax.set_xlim(-0.5,0.4)
panel(axs[0],'a'); panel(axs[1],'b'); panel(axs[2],'c'); fig.tight_layout(w_pad=0.8); save(fig, 'fig2_attenuation')
# ---------------- Fig 3 rivals, dose, gain/loss
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.6), gridspec_kw={'width_ratios':[1.0,0.9,1.25]})
rv = pd.read_csv(R+'x2b_rivals.csv'); lab = ['R4 base','+ CDKN2A loss\n(R1/R2)','+ 9p arm loss\n(R3)','+ both','+ IFN-cluster\nloss (R1)']
dB = (rv.BIC - rv.BIC.iloc[0]).values; ax[0].barh(lab[::-1], dB[::-1], color=[OI['blue']]*4+[OI['grey']]); ax[0].set_xlabel('$\\Delta$BIC vs. base (higher = worse)'); panel(ax[0],'a')
dr = pd.read_csv(R+'x3_dose_cdkn2a.csv', index_col=0); order = ['[T.-2]','[T.-1]','[T.1]']; names = ['Deep loss (-2)','Shallow loss (-1)','Gain (+1/+2)']
rows = [dr[dr.index.str.contains(o, regex=False)].iloc[0] for o in order]; d = pd.DataFrame(rows); d['lab'] = [f'{n}\n(n={int(k)})' for n,k in zip(names, d.n)]
forest(ax[1], d, 'lab', 'coef', 'lo', 'hi', color=[OI['verm'],OI['orange'],OI['green']], xlab='vs. neutral (call 0)'); ax[1].set_xlim(-0.35,0.35); panel(ax[1],'b'); ax[1].set_title('CDKN2A copy state', fontsize=8)
gl = pd.read_csv(R+'x3_gain_loss.csv'); names = ['9p21.3 (CDKN2A)','B2M','STAT1','JAK1','HLA-A','PTEN (control)','RB1 (control)']; yy = np.arange(len(names))[::-1]
ax[2].axvspan(-0.10,0.10,color='#e5e5e5',lw=0); ax[2].axvline(0,color='k',lw=0.6)
for k,(kind,c,off) in enumerate([('loss',OI['blue'],0.17),('gain',OI['orange'],-0.17)]):
    for y_,n_ in zip(yy, names):
        r = gl[(gl.unit==n_)&(gl.kind==kind)]
        if len(r): r = r.iloc[0]; ax[2].plot([r.lo95,r.hi95],[y_+off]*2,color=c,lw=1.3); ax[2].plot(r.beta,y_+off,'o',color=c,ms=3.5, label=kind if y_==yy[0] else None)
ax[2].set_yticks(yy); ax[2].set_yticklabels([n.replace(' (control)','*') for n in names]); ax[2].set_xlabel('logit(LF) effect (M3)'); ax[2].legend(frameon=False, loc='upper center', bbox_to_anchor=(0.5,-0.28), ncol=2); panel(ax[2],'c'); ax[2].set_xlim(-0.3,0.3)
fig.tight_layout(w_pad=1.0); save(fig, 'fig3_rivals_dose_gain')
# ---------------- Fig 4 pan-cancer positive control
x8 = pd.read_csv(R+'x8_percancer_9p.csv'); x8 = x8[x8.b3.notna()].sort_values('b3'); meta = pd.read_csv(R+'x8_meta.csv', index_col=0)
fig, ax = plt.subplots(1, 2, figsize=(7.2, 4.4), gridspec_kw={'width_ratios':[1.5,1]})
cols = [OI['verm'] if r.positive_control else (OI['green'] if r.ct=='BRCA' else OI['blue']) for _,r in x8.iterrows()]
d = x8.copy(); d['lab'] = [f'{c} ({int(l)}/{int(n)})' for c,l,n in zip(d.ct, d.n_loss, d.n)]; forest(ax[0], d, 'lab', 'b3', 'lo3', 'hi3', color=cols, xlab='logit(LF) effect of 9p21.3 loss (M3)'); ax[0].set_xlim(-0.8,0.8); panel(ax[0],'a')
ax[0].text(0.02,0.02,'orange: a-priori positive controls   green: BRCA', transform=ax[0].transAxes, fontsize=6.5, color=OI['grey'])
mm = meta.loc[['all cancers','positive controls','all except positive controls & BRCA']]; mm['lab'] = ['All cancers (k=26)','Positive controls (k=4)','Other cancers (k=21)']
forest(ax[1], mm.reset_index(drop=True), 'lab', 'beta', 'lo', 'hi', color=[OI['grey'],OI['verm'],OI['blue']], xlab='random-effects pooled effect'); ax[1].set_xlim(-0.25,0.25); ax[1].set_title('DerSimonian-Laird meta-analysis', fontsize=8); panel(ax[1],'b')
fig.tight_layout(w_pad=1.2); save(fig, 'fig4_pancancer')
# ---------------- Fig 5 breakpoint classes
c1 = pd.read_csv(R+'x4_class_counts_thr1.csv', index_col=0) if False else pd.read_csv(R+'x4_class_counts_thr1.csv', index_col=0)
c2 = pd.read_csv(R+'x4_class_counts_thr2.csv', index_col=0)
fig, ax = plt.subplots(1, 3, figsize=(7.2, 3.0), gridspec_kw={'width_ratios':[1.5,1.0,1.1]})
tab = c2[['C1 IFN+MTAP+CDKN2A lost','C2 CDKN2A lost, IFN intact']].copy(); tab['tot'] = tab.sum(axis=1); tab = tab.sort_values('tot', ascending=True).tail(20)
ax[0].barh(tab.index, tab.iloc[:,0], color=OI['blue'], label='C1: IFN + MTAP + CDKN2A lost'); ax[0].barh(tab.index, tab.iloc[:,1], left=tab.iloc[:,0], color=OI['orange'], label='C2: CDKN2A lost, IFN intact'); ax[0].set_xlabel('Tumours with deep deletion (GISTIC -2)'); ax[0].legend(frameon=False, loc='center right', fontsize=6.5); panel(ax[0],'a')
ce = pd.read_csv(R+'x4b_classes_pos.csv'); ce = ce[ce.cls.notna()]
for k,(thr,tt) in enumerate([(-1,'shallow or deeper'),(-2,'deep only')]):
    a_ = ax[1+k]; d = ce[(ce.thr==thr)&(ce.set.isin(['pooled positive controls','pooled all other cancers']))].copy(); d['lab'] = d.set.map({'pooled positive controls':'Positive\ncontrols','pooled all other cancers':'Other\ncancers'}) + ' ' + d.cls.map({'C1':'C1','C2':'C2'})
    forest(a_, d, 'lab', 'beta', 'lo', 'hi', color=[OI['blue'] if c=='C1' else OI['orange'] for c in d.cls], xlab='logit(LF) effect vs. neither'); a_.set_title(f'Calls: {tt}', fontsize=8); a_.set_xlim(-0.3,0.3)
    if k==1: a_.set_yticklabels([]) if False else None
panel(ax[1],'b'); panel(ax[2],'c'); fig.tight_layout(w_pad=1.0); save(fig, 'fig5_breakpoint_classes')
# ---------------- Fig 6 genome scan
sc = pd.read_csv(R+'x5_scan.csv'); fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.7), gridspec_kw={'width_ratios':[2.1,0.9,0.9]})
s = sc[sc.kind=='loss'].sort_values('pos'); ar = s.arm.str.extract(r'^(\d+|X)')[0]; ar = ar.replace('X','23').astype(int)
col = np.where(ar%2==0, OI['grey'], '#bdbdbd'); ax[0].scatter(np.arange(len(s)), s.t_M3, s=6, c=col, lw=0); h = s.band=='9p21.3'; ax[0].scatter(np.arange(len(s))[h.values], s.t_M3[h], s=22, color=OI['verm'], zorder=3); ax[0].annotate('9p21.3', (np.arange(len(s))[h.values][0], s.t_M3[h].iloc[0]), xytext=(4,8), textcoords='offset points', color=OI['verm'], fontsize=7.5)
ax[0].axhline(0,color='k',lw=0.5); ax[0].set_xlabel('Cytoband (genomic order, chr1 to X)'); ax[0].set_ylabel('t, loss effect (M3)'); panel(ax[0],'a')
obs = np.sort(s.perm_p.values); exp = (np.arange(1,len(obs)+1))/(len(obs)+1); ax[1].plot([0,4],[0,4],'k',lw=0.6); ax[1].scatter(-np.log10(exp), -np.log10(obs), s=5, color=OI['blue'], lw=0); ax[1].scatter(-np.log10(exp[obs==s.loc[h,'perm_p'].iloc[0]][:1]), -np.log10([s.loc[h,'perm_p'].iloc[0]]), s=18, color=OI['verm'], zorder=3); ax[1].set_xlabel('expected $-\\log_{10}$ p'); ax[1].set_ylabel('observed (perm.)'); panel(ax[1],'b')
sy = pd.read_csv(R+'x5b_joint_loss_gain.csv'); ax[2].scatter(sy.bL, sy.bG, s=4, color=OI['grey'], alpha=0.6, lw=0); r9 = sy[sy.band=='9p21.3']; ax[2].scatter(r9.bL, r9.bG, s=22, color=OI['verm'], zorder=3); ax[2].axhline(0,color='k',lw=0.5); ax[2].axvline(0,color='k',lw=0.5); ax[2].set_xlabel('loss $\\beta$ (joint)'); ax[2].set_ylabel('gain $\\beta$ (joint)'); panel(ax[2],'c')
fig.tight_layout(w_pad=1.0); save(fig, 'fig6_genome_scan')
# ---------------- Fig 7 focal vs arm
fa = pd.read_csv(R+'x6_focal_arm.csv'); fp = pd.read_csv(R+'x6b_arm_pan.csv')
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9), sharey=False)
d = fa[fa.pair=='9p21.3 / 9p'].copy(); d['lab'] = d.model; forest(ax[0], d, 'lab', 'beta', 'lo95', 'hi95', color=[OI['blue'],OI['orange'],OI['blue'],OI['orange']], xlab='logit(LF) effect'); ax[0].set_title('TCGA-BRCA (n=915)', fontsize=8); ax[0].set_xlim(-0.3,0.35); panel(ax[0],'a')
d = fp[fp.set.isin(['pooled positive controls','pooled other cancers (no BRCA)'])].copy(); d['lab'] = d.set.map({'pooled positive controls':'Pos. ctrl','pooled other cancers (no BRCA)':'Other'}) + ': ' + d.model
forest(ax[1], d, 'lab', 'beta', 'lo', 'hi', color=[OI['blue'] if 'gene' in m_ else OI['orange'] for m_ in d.model], xlab='logit(LF) effect'); ax[1].set_title('Pooled across cancers (cancer fixed effects)', fontsize=8); ax[1].set_xlim(-0.2,0.2); panel(ax[1],'b')
fig.tight_layout(w_pad=1.5); save(fig, 'fig7_focal_vs_arm')
# ---------------- Fig 8 subtype heterogeneity
x7 = pd.read_csv(R+'x7_subtype.csv'); mb = pd.read_csv(R+'x12_metabric_subtype.csv')
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.6), gridspec_kw={'width_ratios':[1.15,1.15,1.0]})
for a_, un, t in [(ax[0],'9p21.3 (CDKN2A)','9p21.3 loss'), (ax[1],'Antigen presentation (B2M/HLA-A/B/C/TAP1/2)','Antigen-presentation loss')]:
    d = x7[(x7.unit==un)&(x7.subtype.isin(['LumA','LumB','Basal','Her2']))].drop_duplicates('subtype').copy(); d['lab'] = [f'{s} (n={int(n)})' for s,n in zip(d.subtype, d.n_loss)]
    forest(a_, d, 'lab', 'beta', 'lo95', 'hi95', color=[SUBC[s] for s in d.subtype], xlab='logit(LF) effect (M3)'); a_.set_title(t+' (TCGA)', fontsize=8); a_.set_xlim(-0.55,0.55)
mb['lab'] = [f'{s} (n={int(n)})' for s,n in zip(mb.subtype, mb.n_loss)]; forest(ax[2], mb, 'lab', 'beta', 'lo', 'hi', color=[SUBC[s] for s in mb.subtype], sesoi=False, xlab='CYT effect, SD units'); ax[2].set_title('9p21.3 loss (METABRIC)', fontsize=8); ax[2].set_xlim(-0.6,0.6)
panel(ax[0],'a'); panel(ax[1],'b'); panel(ax[2],'c'); fig.tight_layout(w_pad=1.0); save(fig, 'fig8_subtype')
# ---------------- Fig 9 mechanism / external
x10 = pd.read_csv(R+'x10_rna_outcomes.csv'); x11 = pd.read_csv(R+'x11_ccle.csv'); x12 = pd.read_csv(R+'x12_metabric.csv')
fig, ax = plt.subplots(1, 3, figsize=(7.2, 3.0), gridspec_kw={'width_ratios':[1.2,1.0,1.0]})
names = {'cis_CDKN2A':'CDKN2A expr. (cis)','cis_MTAP':'MTAP expr. (cis)','CYT':'Cytolytic score','CD8':'CD8A/B','IFNG_sig':'IFN-g signature','PTPRC':'PTPRC (CD45)','ISG':'ISG score','MHC_I':'MHC-I score','CXCL9_10':'CXCL9/10'}
d = x10[x10.outcome.isin(names)].copy(); d['lab'] = d.outcome.map(names); d = d.set_index('outcome').loc[list(names)].reset_index(); d['lab'] = d.outcome.map(names)
forest(ax[0], d, 'lab', 'M3_beta', 'M3_lo', 'M3_hi', color=[OI['verm']]*2+[OI['blue']]*7, sesoi=False, xlab='SD units (TCGA-BRCA, M3)'); ax[0].set_title('9p21.3 loss, TCGA-BRCA RNA', fontsize=8); panel(ax[0],'a')
d = x11[(x11.outcome.isin(['ISG','MHC_I']))&(x11.exposure=='x')].copy(); d['lab'] = d.set.map({'Breast lines':'Breast lines','All lineages (lineage FE)':'All lineages'}) + ': ' + d.outcome.map({'ISG':'ISG','MHC_I':'MHC-I'}) + f'\n(n={0})'; d['lab'] = [f"{s}: {o}\n(n lost {int(nl)}/{int(n)})" for s,o,nl,n in zip(d.set.map({'Breast lines':'Breast','All lineages (lineage FE)':'All'}), d.outcome.map({'ISG':'ISG','MHC_I':'MHC-I'}), d.n_loss, d.n)]
forest(ax[1], d, 'lab', 'beta', 'lo', 'hi', color=[OI['blue'],OI['orange'],OI['blue'],OI['orange']], sesoi=False, xlab='SD units (CCLE, adj. FGA)'); ax[1].set_title('Cell lines: CDKN2A deep loss', fontsize=8); panel(ax[1],'b')
d = x12[(x12.exposure=='x')&x12.outcome.isin(['CYT','CD8','IFNG_sig','PTPRC','ISG','MHC_I'])].copy(); d['lab'] = d.outcome.map(names)
forest(ax[2], d, 'lab', 'M3_beta', 'M3_lo', 'M3_hi', color=[OI['green']]*len(d), sesoi=False, xlab='SD units (METABRIC, M3)'); ax[2].set_title('9p21.3 loss, METABRIC (n=1,608)', fontsize=8); panel(ax[2],'c')
fig.tight_layout(w_pad=1.0); save(fig, 'fig9_mechanism_external')
# ---------------- Fig 10 robustness
rb = pd.read_csv(R+'x14_robust.csv'); fig, ax = plt.subplots(figsize=(5.2, 3.6)); d = rb.copy(); d['lab'] = d.spec; forest(ax, d, 'lab', 'beta', 'lo95', 'hi95', color=[OI['verm']]+[OI['blue']]*(len(d)-1)); ax.set_xlim(-0.2,0.35); save(fig, 'fig10_robustness')
print('figures done')
