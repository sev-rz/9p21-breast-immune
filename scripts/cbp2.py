import json,urllib.request,pandas as pd
B='https://www.cbioportal.org/api'
S='brca_tcga_pan_can_atlas_2018'
pt=pd.DataFrame(json.load(urllib.request.urlopen(urllib.request.Request(f'{B}/studies/{S}/clinical-data?clinicalDataType=PATIENT&projection=SUMMARY&pageSize=100000',headers={'Accept':'application/json'}),timeout=120)))
pw=pt.pivot_table(index='patientId',columns='clinicalAttributeId',values='value',aggfunc='first'); pw.to_csv('data/raw/cbp_patient.csv')
print(pw.columns.tolist()); print(pw['SUBTYPE'].value_counts(dropna=False).to_dict())
m=pd.read_csv('data/raw/cbp_gistic.csv',index_col=0)
sub=pd.Series(m.index.str[:12],index=m.index).map(pw['SUBTYPE'])
out=pd.DataFrame({g:(m[g]<=-1).groupby(sub).mean() for g in m.columns}).round(3); out['n']=sub.value_counts()
print(out.T.to_string())
print('ALL <=-1:',(m<=-1).mean().round(3).to_dict()); print('ALL deep -2:',(m<=-2).mean().round(3).to_dict())
# co-deletion: 9p21 any IFN-cluster loss
ifn=(m[['IFNA1','IFNB1','IFNE','IFNW1']]<=-1).any(axis=1); cd=(m[['CDKN2A','CDKN2B','MTAP']]<=-1).any(axis=1)
print('IFN-cluster loss n',ifn.sum(),'CDKN2A/B/MTAP loss n',cd.sum(),'both',(ifn&cd).sum(),'IFN only',(ifn&~cd).sum())
print(pd.crosstab(ifn,sub))
