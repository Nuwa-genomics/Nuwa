import scvi
import scanpy as sc
import pandas as pd


class solo_model():

    def __init__(self, adata):
        self.adata = adata
        self.adata.var_names_make_unique()

    def train_vae(self):
        scvi.model.SCVI.setup_anndata(self.adata)
        self.vae = scvi.model.SCVI(self.adata)
        self.vae.train(max_epochs=400, use_gpu=True)

    def train_solo(self):
        self.solo = scvi.external.SOLO.from_scvi_model(self.vae)
        self.solo.train(max_epochs=400, use_gpu=True)

    def predict_solo(self):
        df = self.solo.predict()
        df['prediction'] = self.solo.predict(soft=False)
        df.index = df.index.map(lambda x: x[:-2])
        self.adata.obs['prediction'] = df.prediction.values

    def cluster_analysis(self):
        sc.pp.normalize_total(self.adata, target_sum=1e4)
        sc.pp.log1p(self.adata)
        sc.tl.pca(self.adata)
        sc.pp.neighbors(self.adata)
        sc.tl.umap(self.adata)
        sc.tl.leiden(self.adata, resolution=0.5)

    def train(self):
        self.train_vae()
        self.train_solo()
        self.predict_solo()
        self.cluster_analysis()

    def get_umap_plt(self):
        self.umap_plt = sc.pl.umap(self.adata, color = ['leiden', 'prediction'])
        return self.umap_plt
    
    def get_adata(self):
        return self.adata