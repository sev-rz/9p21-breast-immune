import json,urllib.request,pandas as pd
B='https://www.cbioportal.org/api'
def post(u,b):
    r=urllib.request.Request(B+u,json.dumps(b).encode(),{'Content-Type':'application/json','Accept':'application/json'});return json.load(urllib.request.urlopen(r,timeout=120))
def get(u):
    return json.load(urllib.request.urlopen(urllib.request.Request(B+u,headers={'Accept':'application/json'}),timeout=120))
S='brca_tcga_pan_can_atlas_2018'
cl=pd.DataFrame(get(f'/studies/{S}/clinical-data?clinicalDataType=SAMPLE&projection=SUMMARY&pageSize=100000'))
w=cl.pivot_table(index='sampleId',columns='clinicalAttributeId',values='value',aggfunc='first')
print(w.columns.tolist()); w.to_csv('data/raw/cbp_clinical.csv')
print(w.get('SUBTYPE').value_counts(dropna=False).to_dict() if 'SUBTYPE' in w else '')
genes={'CDKN2A':1029,'CDKN2B':1030,'MTAP':4507,'IFNA1':3439,'IFNB1':3456,'IFNE':338376,'IFNW1':3467,'B2M':567,'HLA-A':3105,'HLA-B':3106,'TAP1':6890,'JAK1':3716,'JAK2':3717,'STAT1':6772,'IFNGR1':3459,'PTEN':5728,'RB1':5925,'MYC':4609}
d=post(f'/molecular-profiles/{S}_gistic/molecular-data/fetch',{'sampleListId':f'{S}_all','entrezGeneIds':list(genes.values())})
df=pd.DataFrame(d); df['g']=df.entrezGeneId.map({v:k for k,v in genes.items()})
m=df.pivot_table(index='sampleId',columns='g',values='value',aggfunc='first'); m.to_csv('data/raw/cbp_gistic.csv')
print(len(m)); 
sub=w['SUBTYPE'].reindex(m.index)
out=pd.DataFrame({g:(m[g]<=-1).groupby(sub).mean() for g in m.columns}).round(3); out['n']=sub.value_counts(); print(out.T.to_string())
print('all deep(-2):',(m<=-2).mean().round(3).to_dict()); print('all <=-1:',(m<=-1).mean().round(3).to_dict())
