import numpy as np
import matplotlib.pyplot as plt
from torch.utils.data import DataLoader
import time

from ml.citeseq.dataset import TabularDataset
from ml.citeseq.train import get_encodings
from ml.citeseq.model import CiteAutoencoder
from ml.solo_scvi.solo_model import solo_model
from ml.DeepST.deepst.main import *

import pickle
import pandas as pd

import scanpy as sc
import streamlit as st
import umap.umap_ as umap

from components.sidebar import *
from models.AdataModel import AdataModel
from database.schemas import schemas
from utils.AdataState import AdataState
from database.database import SessionLocal

st.set_page_config(layout="wide", page_title='Nuwa', page_icon='🧬')

with open('css/cluster.css') as f:
    st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)

with open('css/common.css') as f:
    st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)

common_style = """
    <style>
    footer {visibility: hidden;}
    .st-emotion-cache-1cypcdb {background: linear-gradient(180deg, rgb(5, 39, 103) 0%, #3a0647 70%); box-shadow: 1px 0 10px -2px #000;}
    </style>
"""
st.markdown(common_style, unsafe_allow_html=True)


class Analysis:
    def __init__(self, adata):

        st.title("Analysis")

        self.col1, self.col2 = st.columns(2, gap="large")

        self.adata = adata
        
        self.columns = adata.to_df().columns

        self.conn: SessionLocal = SessionLocal()


    def save_adata(self, name):
        self.save_adata_to_session(name)
        self.save_adata_to_db(name)

    def save_adata_to_session(self, name):
        for i, adata in enumerate(st.session_state.adata):
            if adata.adata_name == name:
                st.session_state.adata[i] = AdataModel(work_id=st.session_state.current_workspace.id, id=i, adata_name=name, filename=f"{name}.h5ad", adata=self.adata)
                return
        st.session_state.adata.append(AdataModel(work_id=st.session_state.current_workspace.id, id=len(st.session_state.adata), adata_name=name, filename=f"{name}.h5ad", adata=self.adata))
        

    def save_adata_to_db(self, name):
        try:
            new_adata = schemas.Adata(
                work_id=int(st.session_state.current_workspace.id),
                adata_name=f"{name}",
                filename=f"/streamlit-volume/{st.session_state.current_workspace.id}/{name}.h5ad",
                notes="noteeesss"
            )
            self.conn.add(new_adata)
            self.conn.commit()
            self.conn.refresh(new_adata)
        except Exception as e:
            print(e)


    def autoencoder_cluster_plot(self):
        with self.col1:
            try:
                with st.spinner(text='Preparing data from model'):

                    rna_df_original = self.adata.to_df()
                    rna_df = rna_df_original.reset_index(drop=True)

                    test_ds = TabularDataset(rna_df.to_numpy(dtype=np.float32))
                    test_dl = DataLoader(test_ds, batch_size=64, shuffle=False)

                    if "trained_model" not in st.session_state:
                        st.error("Model hasn't been trained to provide data")
                        st.link_button(label="Train", url='/Create_model')
                    else:
                        
                        trained_model = st.session_state["trained_model"]
                        device = st.session_state["device"]

                        if isinstance(trained_model, CiteAutoencoder):
                            st.subheader("Autoencoder Clusters")
                            sc.pp.neighbors(self.adata, n_neighbors=10, n_pcs=40)
                            sc.tl.leiden(self.adata)

                            encodings = get_encodings(trained_model, test_dl, device)
                            encodings = encodings.cpu().numpy()
                            metadata_df = self.adata.obs
                            cell_ids = rna_df_original.index.values
                            plot_df = metadata_df.loc[metadata_df.index.isin(cell_ids)]
                            embedding = umap.UMAP(random_state=0).fit_transform(encodings)
                            
                            plot_df["UMAP1"] = embedding[:, 0]
                            plot_df["UMAP2"] = embedding[:, 1]

                            
                            colors_var = reversed(self.adata.obs.columns)

                            def update_autoencoder_colors():
                                st.session_state.update()
                            
                            st.selectbox(label="Colour", options=(colors_var), key="sb_auto_colors", on_change=update_autoencoder_colors)
                            st.scatter_chart(plot_df, x="UMAP1", y="UMAP2", color=st.session_state['sb_auto_colors'], size=18)

                        elif(isinstance(trained_model, solo_model)):
                            st.subheader("Doublet prediction")
                            def filter_out_doublets():
                                self.adata = self.adata[self.adata.obs.prediction == 'singlet']
                                self.save_adata(name="adata_solo")
                                st.toast("Removed doublet predictions", icon='✅')
                            
                            df_solo = pd.DataFrame({'umap1': self.adata.obsm['X_umap'][:,0], 'umap2': self.adata.obsm['X_umap'][:,1], 'color': self.adata.obs['prediction']})
                          
                            st.scatter_chart(data=df_solo, x='umap1', y='umap2', color='color', size=18)
                            st.button(label="Filter out doublets", key="btn_filter_doublets", on_click=filter_out_doublets)
                            st.divider()

                        elif(isinstance(trained_model, DeepSTModel)):
                            st.subheader("DeepST model")
                            ax_df = trained_model.get_adata_df()
                            st.scatter_chart(ax_df, x='fr1', y='fr2', color='color', height=600)


                        else:
                            st.error("Unknown model")

            except Exception as e:
                st.error(e)


    def pca_graph(self):
        with self.col2:
            try:
                st.subheader("PCA")
                plt.style.use('dark_background')
                st.multiselect(label="Gene", options=self.columns, default=(self.columns[0], self.columns[1]), key="ms_pca_gene", max_selections=24)
                sc.tl.pca(self.adata, svd_solver='arpack')
                ax_pca = sc.pl.pca(self.adata, color=st.session_state.ms_pca_gene or self.columns[0])
                with st.expander(label="Show Figure"):
                    st.pyplot(ax_pca)
                st.divider()
            except Exception as e:
                st.error(e)


    def tsne_graph(self):
        with self.col2:
            try:
                st.subheader("tSNE")
                with st.spinner(text="Running tSNE"):
                    perplexity = st.slider(label="Perplexity", min_value=1, max_value=100, value=30)
                    sc.tl.tsne(self.adata, perplexity=perplexity)
                    colors_var = reversed(self.adata.obs.columns)
                    tsne_color = st.selectbox(label="Colour", options=(colors_var), key="sb_tsne_colors")
                    df_tsne = pd.DataFrame({'tsne1': self.adata.obsm['X_tsne'][:,0], 'tsne2': self.adata.obsm['X_tsne'][:,1], 'color': self.adata.obs[f'{tsne_color}']})  
                    st.scatter_chart(data=df_tsne, x='tsne1', y='tsne2', color='color', size=18)
                st.divider()
            except Exception as e:
                st.error(e)


    def variance_ratio_graph(self):
        with self.col2:
            try:
                st.markdown("### Variance ratio")
                plt.style.use('default')
                ax_variance_ratio = sc.pl.pca_variance_ratio(self.adata, log=True)
                subcol1, _ = st.columns(2)
                with subcol1:
                    with st.expander(label="Show Figure"):
                        st.pyplot(ax_variance_ratio)
                st.divider()
            except Exception as e:
                st.error(e)


    def neighbourhood_graph(self):
        plt.style.use('dark_background')
        with self.col2:
            try:
                st.subheader("Neighbourhood graph")
                with st.spinner(text="Computing neighbourhood graph"):
                    sc.pp.neighbors(self.adata, n_neighbors=10, n_pcs=40)
                    sc.tl.leiden(self.adata) 
                    sc.tl.paga(self.adata)
                    sc.pl.paga(self.adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
                    sc.tl.umap(self.adata, init_pos='paga')

                    umap_options = np.append(self.columns, 'leiden')
                    selected_genes_umap = []
                    st.multiselect(label='Gene', options=(umap_options), key="ms_umap_select_gene", default=["leiden", self.columns[0]], max_selections=24)  
                    
                    ax_umap = sc.pl.umap(self.adata, color=st.session_state.ms_umap_select_gene or 'leiden')

                    with st.expander(label="Show Figure"):
                        st.pyplot(ax_umap)

            except Exception as e:
                st.error(e)


    def find_marker_genes(self):
        with self.col1:
            try:
                st.subheader("Find marker genes")
                st.radio(label="Algorithm", options=(['logreg', 't-test', 'wilcoxon']), key='rb_algo')
                with st.spinner(text="Grouping marker genes"):
                    with st.expander(label="Show Figure"): 
                        sc.tl.rank_genes_groups(self.adata, groupby='leiden', method = str.lower(st.session_state.rb_algo))
                        ax = sc.pl.rank_genes_groups(self.adata, n_genes=25, sharey=False)
                        st.pyplot(ax)

            except Exception as e:
                st.error(e)


try:
    adata_state = AdataState(workspace_id=st.session_state.current_workspace.id)

    sidebar = Sidebar()
    sidebar.show()

    analysis = Analysis(adata_state.current.adata)

    analysis.autoencoder_cluster_plot()
    analysis.neighbourhood_graph()
    analysis.tsne_graph()
    analysis.pca_graph()
    analysis.variance_ratio_graph()
    analysis.find_marker_genes()

    sidebar.show_preview()
    sidebar.delete_experiment_btn()

except KeyError as ke:
    print("KeyError: ", ke)
    st.error("Couldn't find adata object in session, have you uploaded one?")



