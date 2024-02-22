import scanpy as sc

adata = sc.read_h5ad('covid_GSE149689_modified.h5ad')

print(adata.obs['sample'].unique())
print(adata.obs.type.unique())