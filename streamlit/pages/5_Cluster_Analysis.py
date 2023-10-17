import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from torch.utils.data import DataLoader
import seaborn as sns

from ml.citeseq.dataset import TabularDataset
from ml.citeseq.train import get_encodings
import pickle

import scanpy as sc
import streamlit as st
import umap.umap_ as umap

st.set_page_config(layout="wide", page_title='Nuwa', page_icon='🧬')

common_style = """
            <style>
            footer {visibility: hidden;}
            .st-emotion-cache-1cypcdb {background: linear-gradient(180deg, rgb(5, 39, 103) 0%, #3a0647 70%); box-shadow: 1px 0 10px -2px #000;}
            .st-emotion-cache-86cver {rgba(250, 250, 250, 0.6)}
            </style>
            """
st.markdown(common_style, unsafe_allow_html=True)

if 'adata' not in st.session_state:
    tmp_file = open("./tmp/adata.pkl",'rb')
    cached_adata = pickle.load(tmp_file)
    st.session_state["adata"] = cached_adata


class Analysis:
    def __init__(self, adata):

        st.title("Analysis")

        self.col1, self.col2 = st.columns(2, gap="large")

        self.adata = adata

        self.columns = adata.to_df().columns
        

    def autoencoder_cluster_plot(self):
        with self.col1:
            st.subheader("Autoencoder cluster plot")
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

                    encodings = get_encodings(trained_model, test_dl, device)
                    encodings = encodings.cpu().numpy()

                    metadata_df = self.adata.obs

                    cell_ids = rna_df_original.index.values

                    plot_df = metadata_df.loc[metadata_df.index.isin(cell_ids)]

                    embedding = umap.UMAP(random_state=0).fit_transform(encodings)
                    
                    plot_df["UMAP1"] = embedding[:, 0]
                    plot_df["UMAP2"] = embedding[:, 1]

                    colors_var = self.adata.obs.columns

                    def update_autoencoder_colors():
                        st.session_state.update()
                    

                    st.selectbox(label="Colour", options=(colors_var), key="sb_auto_colors", on_change=update_autoencoder_colors)

                    
                    

                    st.scatter_chart(plot_df, x="UMAP1", y="UMAP2", color=st.session_state['sb_auto_colors'])

    def pca_graph(self):
        with self.col2:
            st.subheader("PCA")
            plt.style.use('dark_background')
            st.multiselect(label="Gene", options=self.columns, default=(self.columns[0], self.columns[1]), key="ms_pca_gene", max_selections=24)
            sc.tl.pca(self.adata, svd_solver='arpack')
            ax_pca = sc.pl.pca(self.adata, color=st.session_state.ms_pca_gene or self.columns[0])
            sc.pl.pca(self.adata, color=st.session_state.ms_pca_gene or self.columns[0], save=True)
            with st.expander(label="Show Figure"):
                st.pyplot(ax_pca)
            st.button(label="Save figure", key="btn_pca")


    def variance_ratio_graph(self):
        with self.col2:
            st.markdown("### Variance ratio")
            plt.style.use('default')
            ax_variance_ratio = sc.pl.pca_variance_ratio(self.adata, log=True)
            self.save_figure_variance_ratio = lambda: sc.pl.pca_variance_ratio(self.adata, log=True, save=True)
            subcol1, subcol2 = st.columns(2)
            with subcol1:
                with st.expander(label="Show Figure"):
                    st.pyplot(ax_variance_ratio)
            with subcol2:
                st.button(label="Save figure", key="btn_variance_ratio", on_click=self.save_figure_variance_ratio)


    def neighbourhood_graph(self):
        plt.style.use('dark_background')
        with self.col2:
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

            self.save_neighbourhood_graph = lambda: sc.pl.umap(self.adata, color=st.session_state.ms_umap_select_gene or 'leiden', save=True)
            st.button(label="Save figure", key="btn_neighbours", on_click=self.save_neighbourhood_graph)


    def find_marker_genes(self):
        with self.col1:
            st.subheader("Find marker genes")
            st.radio(label="Algorithm", options=(['logreg', 't-test', 'wilcoxon']), key='rb_algo')
            with st.spinner(text="Grouping marker genes"):
                with st.expander(label="Show Figure"): 
                    sc.tl.rank_genes_groups(self.adata, groupby='leiden', method = str.lower(st.session_state.rb_algo))
                    ax = sc.pl.rank_genes_groups(self.adata, n_genes=25, sharey=False)
                    st.pyplot(ax)
            st.button(label="Save figure", key="btn_marker_gene")
            
            
adata = st.session_state.adata

analysis = Analysis(adata)

analysis.autoencoder_cluster_plot()
analysis.pca_graph()
analysis.variance_ratio_graph()
analysis.neighbourhood_graph()
analysis.find_marker_genes()
