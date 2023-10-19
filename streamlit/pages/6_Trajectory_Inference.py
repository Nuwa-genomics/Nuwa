import scanpy as sc
import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time

st.set_page_config(layout="wide")

with open('css/common.css') as f:
    common_style = f"""
                <style>
                {f.read()}
                </style>
                """
    st.markdown(common_style, unsafe_allow_html=True)

class Trajectory_Inference:
    def __init__(self, adata):
        self.adata = adata
        self.adata.raw = adata
        #convert to float64
        self.adata.X = self.adata.X.astype('float64')
        #compute pca
        sc.tl.pca(self.adata, svd_solver='arpack')
        #compute neighbours
        sc.pp.neighbors(self.adata, n_neighbors=4, n_pcs=20)
        sc.tl.draw_graph(self.adata)
        

    def draw_graph(self):
        with self.col1:
            st.subheader("Cluster chart")
            with st.spinner(text="Plotting clusters"):
                subcol1, _, _, _ = st.columns(4, gap='large')
                subcol1.selectbox(label='color', options=self.adata.obs.columns, key='sb_plt_colors')
                ax = sc.pl.draw_graph(self.adata, color=st.session_state.sb_plt_colors)
                #denoise
                # sc.tl.diffmap(self.adata)
                # sc.pp.neighbors(self.adata, n_neighbors=10, use_rep='X_diffmap')
                # sc.tl.draw_graph(self.adata)
                ax_df = pd.DataFrame({'fr1': self.adata.obsm['X_draw_graph_fr'][:,0], 'fr2': self.adata.obsm['X_draw_graph_fr'][:,1], 'color': self.adata.obs[st.session_state.sb_plt_colors]})
                st.scatter_chart(ax_df, x='fr1', y='fr2', color='color', height=600)

    def draw_graph_paga(self):
        plt.style.use('dark_background')
        with self.col2:
            st.subheader("Louvain with PAGA embedding")
            with st.spinner(text="Recomputing with PAGA initialisation"):
                sc.tl.draw_graph(self.adata, init_pos='paga')
                st.multiselect(label="Gene", options=np.append(self.adata.to_df().columns, 'louvain'), default=['louvain', self.adata.to_df().columns[0]], key="ms_louvain_colors_paga")
                with st.expander(label="Show figure", expanded=True):
                    ax = sc.pl.draw_graph(self.adata, color=st.session_state.ms_louvain_colors, legend_loc='on data')
                    st.pyplot(ax)



    def louvain_cluster(self):
        #plt.style.use('dark_background')
        with self.col1:
            st.subheader("PAGA")
            plt.style.use('classic')
            st.multiselect(label="Gene", options=np.append(self.adata.to_df().columns, 'louvain'), default=['louvain', self.adata.to_df().columns[0]], key="ms_louvain_colors")
            #louvain cluster
            sc.tl.louvain(self.adata, resolution=1.0)

            #paga
            sc.tl.paga(self.adata, groups='louvain')
            ax = sc.pl.paga(self.adata, colors=st.session_state.ms_louvain_colors, cmap='viridis', node_size_scale=3)
            with st.expander(label="Show figure", expanded=True):
                st.pyplot(ax)

    def draw_page(self):
        st.title("Trajectory Inference")
        #columns
        self.col1, self.col2 = st.columns(2, gap="large")


    def dpt(self):
        with self.col2:
            st.subheader("Diffusion Pseudotime")
            st.selectbox(label='Root cell', options=(self.adata.obs['louvain']), key='sb_root_cell')
            with st.expander(label="Show figure", expanded=True):
                self.adata.uns['iroot'] = np.flatnonzero(self.adata.obs['louvain']  == st.session_state.sb_root_cell)[0]
                

                gene_names = ['Gata2', 'Gata1', 'Klf1', 'Epor', 'Hba-a2', 'Elane', 'Cebpe', 'Gfi1', 'Irf8', 'Csf1r', 'Ctsg']    

                adata_raw = self.adata
                sc.tl.dpt(adata_raw)
                sc.pp.log1p(adata_raw)
                sc.pp.scale(adata_raw)

                ax = sc.pl.draw_graph(adata_raw, color=['louvain', 'dpt_pseudotime'], legend_loc='on data', cmap='viridis')
                st.pyplot(ax)





adata = st.session_state.adata

tji = Trajectory_Inference(adata)

tji.draw_page()

tji.draw_graph()

tji.louvain_cluster()

tji.draw_graph_paga()

tji.dpt()