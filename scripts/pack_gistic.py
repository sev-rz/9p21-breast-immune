import pandas as pd, numpy as np
df=pd.read_csv('data/raw/gistic_thr.tsv',sep='\t',dtype={**{'Gene Symbol':str,'Locus ID':str,'Cytoband':str}})
meta=df.iloc[:,:3].copy(); X=df.iloc[:,3:].to_numpy(dtype=np.int8)
samples=np.array(df.columns[3:])
np.save('data/processed/thr.npy',X); meta.to_csv('data/processed/thr_genes.csv',index=False); pd.Series(samples).to_csv('data/processed/thr_samples.csv',index=False,header=['barcode'])
print(X.shape, np.unique(X,return_counts=True))
print(meta[meta['Gene Symbol'].isin(['CDKN2A','MTAP','IFNA1','IFNB1','IFNE','B2M','JAK2','STAT1'])])
