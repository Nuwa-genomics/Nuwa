import os
import tempfile

import matplotlib.pyplot as plt
import scanpy as sc
import scvi
import seaborn as sns
import torch
from lightning.pytorch.callbacks import Callback
import streamlit as st
import pandas as pd
from anndata import AnnData

class LDVAE:

    def __init__(self, adata: AnnData, max_epochs=250, lr=5e-3, check_val_every_n_epoch=10, use_gpu=True):

        self.adata = adata
        self.max_epochs = max_epochs
        self.lr = lr
        self.check_val_every_n_epoch = check_val_every_n_epoch
        self.use_gpu = use_gpu
        scvi.settings.seed = 42

        self.adata.layers["counts"] = adata.X.copy()  # preserve counts

        #sc.pp.normalize_total(adata, target_sum=10e4)
        #sc.pp.log1p(adata)
        #adata.raw = adata  # freeze the state in `.raw`

    def train(self, callbacks=None):
        scvi.model.LinearSCVI.setup_anndata(self.adata, layer="counts")
        self.model = scvi.model.LinearSCVI(self.adata, n_latent=10)

        self.model.train(
            max_epochs=self.max_epochs, 
            plan_kwargs={"lr": self.lr}, 
            check_val_every_n_epoch=self.check_val_every_n_epoch, 
            use_gpu=self.use_gpu, 
            callbacks=callbacks
        )

        return self
    
    def get_embedding(self):
        self.SCVI_LATENT_KEY = "X_scVI"
        self.SCVI_CLUSTERS_KEY = "leiden_scVI"

        Z_hat = self.model.get_latent_representation()

        self.adata.obsm[self.SCVI_LATENT_KEY] = Z_hat
        sc.pp.neighbors(self.adata, use_rep=self.SCVI_LATENT_KEY, n_neighbors=20)
        sc.tl.umap(self.adata, min_dist=0.3)
        sc.tl.leiden(self.adata, key_added=self.SCVI_CLUSTERS_KEY, resolution=0.8)

        df = pd.DataFrame({'UMAP1': self.adata.obsm['X_umap'][:,0], 'UMAP2': self.adata.obsm['X_umap'][:,1]})

        return df
    






# extract latent features

# Z_hat = model.get_latent_representation()
# for i, z in enumerate(Z_hat.T):
#     adata.obs[f"Z_{i}"] = z


# fig = plt.figure(figsize=(12, 8))

# for f in range(0, 9, 2):
#     plt.subplot(2, 3, int(f / 2) + 1)

#     plt.scatter(
#         adata.obs[f"Z_{f}"], adata.obs[f"Z_{f + 1}"], marker=".", s=4, label="Cells"
#     )

#     plt.xlabel(f"Z_{f}")
#     plt.ylabel(f"Z_{f + 1}")

# plt.subplot(2, 3, 6)
# plt.scatter(
#     adata.obs[f"Z_{f}"], adata.obs[f"Z_{f + 1}"], marker=".", label="Cells", s=4
# )
# plt.scatter(adata.obs[f"Z_{f}"], adata.obs[f"Z_{f + 1}"], c="w", label=None)
# plt.gca().set_frame_on(False)
# plt.gca().axis("off")

# lgd = plt.legend(scatterpoints=3, loc="upper left")
# for handle in lgd.legendHandles:
#     handle.set_sizes([200])


# plt.tight_layout()


# SCVI_LATENT_KEY = "X_scVI"
# SCVI_CLUSTERS_KEY = "leiden_scVI"

# Z_hat = model.get_latent_representation()

# adata.obsm[SCVI_LATENT_KEY] = Z_hat
# sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY, n_neighbors=20)
# sc.tl.umap(adata, min_dist=0.3)
# sc.tl.leiden(adata, key_added=SCVI_CLUSTERS_KEY, resolution=0.8)

# sc.pl.umap(adata, color=[SCVI_CLUSTERS_KEY])