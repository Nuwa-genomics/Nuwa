import scvi
import scanpy as sc
import streamlit as st
import pandas as pd

class solo_model():
    
    def __init__(self, adata, epochs=400, lr=1e-3, train_size=0.9, use_gpu=True):
        self.adata = adata
        self.epochs = epochs
        self.lr = lr
        self.train_size = train_size
        self.use_gpu=use_gpu
        self.adata.var_names_make_unique()

    def train_vae(self, callbacks=None):
        scvi.model.SCVI.setup_anndata(self.adata)
        self.vae = scvi.model.SCVI(self.adata)
        self.vae.train(max_epochs=self.epochs, use_gpu=self.use_gpu, train_size=self.train_size, callbacks=callbacks)

    def train_solo(self, callbacks=None):
        self.solo = scvi.external.SOLO.from_scvi_model(self.vae)
        self.solo.train(max_epochs=self.epochs, lr=self.lr, use_gpu=self.use_gpu, train_size=self.train_size, callbacks=callbacks)


    def predict_solo(self):
        df = self.solo.predict()
        df['solo_prediction'] = self.solo.predict(soft=False)
        df.index = df.index.map(lambda x: x[:-2])
        self.adata.obs['solo_prediction'] = df.solo_prediction.values

    def train(self, callbacks):
        self.train_vae(callbacks=callbacks)
        self.train_solo(callbacks=callbacks)
        self.predict_solo()
        return self

    def get_umap_plt(self):
        self.umap_plt = sc.pl.umap(self.adata, color = ['leiden', 'solo_prediction'])
        return self.umap_plt
    
    def get_adata(self):
        return self.adata