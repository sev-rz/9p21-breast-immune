import sys, warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'scripts')
import numpy as np, pandas as pd
R='results/'; T='report/tables/'
def f(x, d=3):
    return '--' if pd.isna(x) else f'{x:.{d}f}'.replace('-','$-$')
def ci(b,l,h,d=3): return f'{f(b,d)} [{f(l,d)}, {f(h,d)}]'
def pf(p): return '--' if pd.isna(p) else ('<0.001' if p<0.001 else f'{p:.3f}')
def tex(rows, header, align, name, caption, label, size='\\small', note=None):
    s = '\\begin{table}[t]\n\\centering\n' + size + '\n\\begin{tabular}{' + align + '}\n\\toprule\n' + ' & '.join(header) + ' \\\\\n\\midrule\n' + '\n'.join(' & '.join(r) + ' \\\\' for r in rows) + '\n\\bottomrule\n\\end{tabular}\n'
    s += '\\caption{' + caption + '}\n\\label{' + label + '}\n' + ('\\vspace{2pt}\n{\\footnotesize ' + note + '}\n' if note else '') + '\\end{table}\n'
    open(T+name+'.tex','w').write(s.replace('_','\\_').replace('&amp;','\\&') if False else s)
esc = lambda t: str(t).replace('_','\\_').replace('%','\\%').replace('&','\\&').replace('#','\\#')
# T1: ladder primary rows
x2 = pd.read_csv(R+'x2_ladder.csv'); rows=[]
for u in ['9p21.3 (CDKN2A)','MTAP','IFN cluster (any)','Antigen presentation (B2M/HLA-A/B/C/TAP1/2)','IFN signalling (JAK1/2,STAT1,IRF1,IFNGR1/2)','B2M','STAT1','JAK1','HLA-A','TAP1','JAK2','PTEN (control)','RB1 (control)']:
    d = x2[x2.unit==u].set_index('model'); r3 = d.loc['M3 +purity']
    rows.append([esc(u), str(int(r3.n_loss)), ci(d.loc['M0 marginal'].beta, d.loc['M0 marginal'].lo95, d.loc['M0 marginal'].hi95,2), ci(d.loc['M1 +subtype'].beta, d.loc['M1 +subtype'].lo95, d.loc['M1 +subtype'].hi95,2), ci(r3.beta, r3.lo95, r3.hi95,2), esc(r3.cls)])
tex(rows, ['Locus unit','$n$ lost','Marginal','+ subtype','+ subtype, burden, purity (M3)','Pre-registered class'], 'lrllll', 'tab_ladder', 'Attenuation ladder in TCGA-BRCA ($n=915$). Coefficients are effects of loss (GISTIC $\\le -1$) on logit(LF), with 95\\% confidence intervals (HC3). Burden is the leave-arm-out gene-count fraction altered.', 'tab:ladder', '\\footnotesize', 'Classes follow the pre-registration (SESOI $=0.10$). Asterisked controls are loci with no expected immune function.')
# T2: join counts
m = pd.read_pickle('data/processed/panc.pkl'); b = m[m.ct=='BRCA']
rows = [['PanImmune leukocyte fraction estimates for TCGA-BRCA (all samples)','1,101'],['\\quad primary tumours (sample type 01), one per sample','1,084'],['\\quad with GISTIC calls','%d' % len(b)],['\\quad with ABSOLUTE purity','%d' % b.purity.notna().sum()],['\\quad with PAM50 subtype in \\{LumA, LumB, Basal, Her2\\}','%d' % b['sub'].isin(['LumA','LumB','Basal','Her2']).sum()],['\\quad modelled (subtype, purity, AS, ploidy, WGD all present)','915'],['Pan-cancer analysis set (29 cancer types with LF and GISTIC)','%s' % f'{len(m):,}']]
tex(rows, ['Step','$n$'], 'lr', 'tab_joins', 'Sample counts at each join.', 'tab:joins')
# T3: per-cancer 9p
x8 = pd.read_csv(R+'x8_percancer_9p.csv').sort_values('ct'); rows=[]
for _,r in x8.iterrows():
    if pd.isna(r.b3): rows.append([r.ct, str(int(r.n)), str(int(r.n_loss)), '--','--','--','too few loss/neutral']); continue
    rows.append([r.ct + ('$^\\dagger$' if r.positive_control else ''), str(int(r.n)), str(int(r.n_loss)), f(r.b0,2), ci(r.b3,r.lo3,r.hi3,2), pf(r.q3), esc(r.cls)])
tex(rows, ['Cancer','$n$','$n$ lost','Marginal','Adjusted (M3)','FDR $q$','Class'], 'lrrrlll', 'tab_percancer', 'Effect of 9p21.3 (CDKN2A) shallow-or-deeper loss on logit(LF) per cancer type, adjusted for purity and leave-9p-out burden. $\\dagger$: a-priori positive-control cancers. SKCM uses metastatic samples.', 'tab:percancer', '\\footnotesize')
# T4: rivals
rv = pd.read_csv(R+'x2b_rivals.csv'); rows=[[esc(r.model), str(int(r.k)), f(r.AIC,1), f(r.BIC,1), f(r.R2,4), f(r.cvR2,4)] for _,r in rv.iterrows()]
tex(rows, ['Model','$k$','AIC','BIC','$R^2$','CV $R^2$'], 'lrrrrr', 'tab_rivals', 'Comparison of rival explanations in TCGA-BRCA. The base model contains subtype, purity, leave-9p-out burden, aneuploidy score, ploidy and genome doubling.', 'tab:rivals')
# T5: robustness
rb = pd.read_csv(R+'x14_robust.csv'); rows=[[esc(r.spec), str(int(r.n)), ci(r.beta,r.lo95,r.hi95,3), pf(r.p), 'yes' if r['equiv_0.10'] else 'no'] for _,r in rb.iterrows()]
tex(rows, ['Specification','$n$','Effect [95\\% CI]','$p$','Equivalent to null?'], 'lrllc', 'tab_robust', 'Robustness of the primary BRCA contrast (9p21.3 loss on logit(LF)). Equivalent means the 90\\% CI lies within $\\pm0.10$.', 'tab:robust', '\\footnotesize')
# T6: class counts
c1 = pd.read_csv(R+'x4_class_counts_thr1.csv', index_col=0).sum(); c2 = pd.read_csv(R+'x4_class_counts_thr2.csv', index_col=0).sum()
rows=[[esc(k), f'{int(c1[k]):,}', f'{int(c2[k]):,}'] for k in c1.index if not k.startswith('C5')]
tex(rows, ['Deletion class','Shallow or deeper ($\\le-1$)','Deep only ($=-2$)'], 'lrr', 'tab_classes', 'Pan-cancer counts of 9p21.3 deletion classes (10,128 tumours).', 'tab:classes')
# T7: subtype
x7 = pd.read_csv(R+'x7_subtype.csv'); rows=[]
for u in ['9p21.3 (CDKN2A)','Antigen presentation (B2M/HLA-A/B/C/TAP1/2)','IFN signalling (JAK1/2,STAT1,IRF1,IFNGR1/2)','STAT1','B2M']:
    d = x7[x7.unit==u]
    for s in ['LumA','LumB','Basal','Her2']:
        r = d[d.subtype==s].iloc[0]; rows.append([esc(u) if s=='LumA' else '', s, str(int(r.n_loss)), ci(r.beta,r.lo95,r.hi95,2), pf(r.p)])
    rr = d[d.subtype.str.startswith('interaction luminal')].iloc[0]; r4 = d[d.subtype.str.startswith('interaction x')].iloc[0]
    rows.append(['', 'interaction (4-way F / lum.\\ vs basal)', '', 'p = ' + pf(r4.p) + ' / ' + pf(rr.p), ''])
tex(rows, ['Locus unit','Subtype','$n$ lost','Effect [95\\% CI]','$p$'], 'llrll', 'tab_subtype', 'Subtype-specific adjusted effects in TCGA-BRCA (M3 within subtype). No correction for multiple testing; the Basal 9p21.3 estimate is exploratory.', 'tab:subtype', '\\footnotesize')
print('tables ok')
