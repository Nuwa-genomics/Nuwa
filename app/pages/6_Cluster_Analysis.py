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
import os

st.set_page_config(layout="wide", page_title='Nuwa', page_icon='🧬')

os.chdir('/app')

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

        self.col1, self.col2 = st.columns(2, gap="medium")

        self.adata = adata
        
        self.columns = adata.to_df().columns

        self.conn: SessionLocal = SessionLocal()


    def save_adata(self, name):
        try:
            new_adata = schemas.Adata(
                work_id=int(st.session_state.current_workspace.id),
                adata_name=f"{name}",
                filename=os.path.join(os.getenv('WORKDIR'), "adata", f"{name}.h5ad"),
                notes=st.session_state.adata_state.current.notes
            )
            self.conn.add(new_adata)
            self.conn.commit()
            self.conn.refresh(new_adata)
        except Exception as e:
            print(e)

    def autoencoder_cluster_plot(self):
        with self.col1:
            try:
                if "trained_model" not in st.session_state:
                    st.error("Model hasn't been trained to provide data")
                else:
                    trained_model = st.session_state["trained_model"]
                    device = st.session_state["device"]

                    if isinstance(trained_model, CiteAutoencoder):
                        with st.form(key="citeseq_form"):
                            st.subheader("Autoencoder Clusters")
                            
                            with st.spinner(text="Preparing embeddings"): 
                            
                                rna_df_original = self.adata.to_df()
                                rna_df = rna_df_original.reset_index(drop=True)

                                test_ds = TabularDataset(rna_df.to_numpy(dtype=np.float32))
                                test_dl = DataLoader(test_ds, batch_size=64, shuffle=False)
                                
                                
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
                                
                                citeseq_container = st.empty()

                                color = st.selectbox(label="Colour", options=(colors_var), key="sb_auto_colors")
                                citeseq_container.scatter_chart(plot_df, x="UMAP1", y="UMAP2", color='leiden', size=18)
                                
                                subcol1, _, _, _ = st.columns(4)
                                submit_btn = subcol1.form_submit_button(label="Update colours")
                                
                                if submit_btn:
                                    with st.spinner(text="Replotting chart"):
                                        citeseq_container.scatter_chart(plot_df, x="UMAP1", y="UMAP2", color=color, size=18)

                    elif(isinstance(trained_model, solo_model)):
                        with st.form(key="solo_doublet_form"):
                            st.subheader("Doublet prediction")
                            
                            with st.spinner("Loading Solo doublet predictions"):
                                
                                df_solo = pd.DataFrame({'umap1': self.adata.obsm['X_umap'][:,0], 'umap2': self.adata.obsm['X_umap'][:,1], 'color': self.adata.obs['prediction']})
                            
                                st.scatter_chart(data=df_solo, x='umap1', y='umap2', color='color', size=18)
                                
                                subcol1, _, _, _ = st.columns(4)
                                submit_btn = subcol1.form_submit_button(label="Filter out doublets", use_container_width=True)
                                
                                if submit_btn:
                                    self.adata = self.adata[self.adata.obs.prediction == 'singlet']
                                    self.save_adata(name="adata_solo")
                                    st.toast("Removed doublet predictions", icon='✅')

                    elif(isinstance(trained_model, DeepSTModel)):
                        with st.spinner(text="Preparing embeddings"):
                            st.subheader("DeepST model")
                            ax_df = trained_model.get_adata_df()
                            st.scatter_chart(ax_df, x='fr1', y='fr2', color='color', height=600)


                    else:
                        st.error("Unknown model")

            except Exception as e:
                st.error(e)


    def pca_graph(self):
        with self.col1:
            try:
                with st.form(key="pca_cluster_form"):
                    st.subheader("PCA")
                    plt.style.use('dark_background')
                    genes = st.multiselect(label="Gene", options=self.columns, default=(self.columns[0], self.columns[1]), key="ms_pca_gene", max_selections=24)
                    pca_container = st.empty()
                    
                    with st.spinner(text="Computing PCA"):
                        sc.tl.pca(self.adata, svd_solver='arpack')
                        ax_pca = sc.pl.pca(self.adata, color=st.session_state.ms_pca_gene or self.columns[0])
                        pca_container.pyplot(ax_pca)
                        #TODO: add variance ratio
                        # with st.expander(label="Show PCA variance ratio"):
                        #     ax_variance_ratio = sc.pl.pca_variance_ratio(self.adata, log=True)
                        #     st.pyplot(ax_variance_ratio)
                        
                    subcol1, _, _, _ = st.columns(4)
                    submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
                    if submit_btn:
                        with st.spinner(text="Computing PCA"):
                            sc.tl.pca(self.adata, svd_solver='arpack')
                            ax_pca = sc.pl.pca(self.adata, color=genes or self.columns[0])
                            pca_container.pyplot(ax_pca)
                            # with st.expander(label="Show PCA variance ratio"):
                            #     ax_variance_ratio = sc.pl.pca_variance_ratio(self.adata, log=True)
                            #     pca_container.pyplot(ax_variance_ratio)
            except Exception as e:
                st.error(e)


    def tsne_graph(self):
        with self.col2:
            try:
                        
                with st.form(key="tsne_form"):
                            
                             
                    st.subheader("tSNE")
                    perplexity = st.slider(label="Perplexity", min_value=1, max_value=100, value=30)
                    colors_var = reversed(self.adata.obs.columns)
                    
                    tsne_color = st.selectbox(label="Colour", options=(colors_var), key="sb_tsne_colors")
                    tsne_container = st.empty()
                    
                    with st.spinner(text="Computing tSNE"):
                        sc.pp.neighbors(self.adata, n_neighbors=10, n_pcs=40)
                        sc.tl.leiden(self.adata) 
                        sc.tl.tsne(self.adata, perplexity=30)   
                        df_tsne = pd.DataFrame({'tsne1': self.adata.obsm['X_tsne'][:,0], 'tsne2': self.adata.obsm['X_tsne'][:,1], 'color': self.adata.obs[f'leiden']})  
                        tsne_container.scatter_chart(data=df_tsne, x='tsne1', y='tsne2', color='color', size=18)  
                    
                    subcol1, _, _, _ = st.columns(4)
                    submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)

                    if submit_btn:
                        sc.tl.tsne(self.adata, perplexity=perplexity)   
                        df_tsne = pd.DataFrame({'tsne1': self.adata.obsm['X_tsne'][:,0], 'tsne2': self.adata.obsm['X_tsne'][:,1], 'color': self.adata.obs[f'{tsne_color}']})  
                        tsne_container.scatter_chart(data=df_tsne, x='tsne1', y='tsne2', color='color', size=18)  
                        
            except Exception as e:
                st.error(e)



    def neighbourhood_graph(self):
        plt.style.use('dark_background')
        with self.col2:
            try:
                with st.form(key="nhood_graph_form"):
                    st.subheader("Neighbourhood graph")
                    
                    umap_options = np.append(self.columns, 'leiden')
                    genes = st.multiselect(label='Gene', options=(umap_options), key="ms_umap_select_gene", default=["leiden", self.columns[0]], max_selections=24) 
                        
                    nhood_container = st.empty() 
                    
                    with st.spinner(text="Computing neighbourhood graph"):
                                
                        sc.pp.neighbors(self.adata, n_neighbors=10, n_pcs=40)
                        sc.tl.leiden(self.adata) 
                        sc.tl.paga(self.adata)
                        sc.pl.paga(self.adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
                        sc.tl.umap(self.adata, init_pos='paga')
                            
                        ax_umap = sc.pl.umap(self.adata, color=genes or 'leiden')

                        nhood_container.pyplot(ax_umap)
                        
                    subcol1, _, _, _ = st.columns(4)
                        
                    submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
                        
                    if submit_btn:
                        with st.spinner(text="Computing neighbourhood graph"):
                                
                            sc.pp.neighbors(self.adata, n_neighbors=10, n_pcs=40)
                            sc.tl.leiden(self.adata) 
                            sc.tl.paga(self.adata)
                            sc.pl.paga(self.adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
                            sc.tl.umap(self.adata, init_pos='paga')
                            
                            ax_umap = sc.pl.umap(self.adata, color=genes or 'leiden')

                            nhood_container.pyplot(ax_umap)

            except Exception as e:
                st.error(e)



try:
    adata = st.session_state.adata_state.current.adata

    sidebar = Sidebar()
    sidebar.show()

    analysis = Analysis(adata)

    analysis.autoencoder_cluster_plot()
    analysis.tsne_graph()
    analysis.neighbourhood_graph()
    analysis.pca_graph()

    sidebar.show_preview()
    sidebar.export_script()
    sidebar.delete_experiment_btn()
    sidebar.show_version()


except KeyError as ke:
    print("KeyError: ", ke)
    st.error("Couldn't find adata object in session, have you uploaded one?")



