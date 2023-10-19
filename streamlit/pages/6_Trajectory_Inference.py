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
        #convert to float64
        self.adata.X = self.adata.X.astype('float64')
        #compute pca
        sc.tl.pca(self.adata, svd_solver='arpack')
        #compute neighbours
        sc.pp.neighbors(self.adata, n_neighbors=4, n_pcs=20)
        sc.tl.draw_graph(self.adata)
        

    def draw_graph(self):
        
        with self.col1:
            with st.spinner(text="Plotting clusters"):
                st.selectbox(label='color', options=self.adata.obs.columns, key='sb_plt_colors')
                ax = sc.pl.draw_graph(self.adata, color=st.session_state.sb_plt_colors)
                ax_df = pd.DataFrame({'fr1': self.adata.obsm['X_draw_graph_fr'][:,0], 'fr2': self.adata.obsm['X_draw_graph_fr'][:,1], 'color': self.adata.obs[st.session_state.sb_plt_colors]})
   
                st.scatter_chart(ax_df, x='fr1', y='fr2', color='color', height=600)

    def draw_graph_paga(self):
        plt.style.use('dark_background')
        with self.col2:
            with st.spinner(text="Recomputing with PAGA initialisation"):
                with st.expander(label="Show figure"):
                    sc.tl.draw_graph(self.adata, init_pos='paga')
                    st.multiselect(label="Gene", options=np.append(self.adata.to_df().columns, 'louvain'), default=['louvain', self.adata.to_df().columns[0]], key="ms_louvain_colors_paga")
                    ax = sc.pl.draw_graph(self.adata, color=st.session_state.ms_louvain_colors, legend_loc='on data')
                    st.pyplot(ax)



    def louvain_cluster(self):
        with self.col2:
            st.subheader("PAGA")
            plt.style.use('classic')
            st.multiselect(label="Gene", options=np.append(self.adata.to_df().columns, 'louvain'), default=['louvain', self.adata.to_df().columns[0]], key="ms_louvain_colors")
            #louvain cluster
            sc.tl.louvain(self.adata, resolution=1.0)

            #paga
            sc.tl.paga(self.adata, groups='louvain')
            ax = sc.pl.paga(self.adata, colors=st.session_state.ms_louvain_colors, cmap='viridis')
            with st.expander(label="Show figure"):
                st.pyplot(ax)

    def draw_page(self):
        st.title("Trajectory Inference")
        #columns
        self.col1, self.col2 = st.columns(2, gap="large")


adata = st.session_state.adata

tji = Trajectory_Inference(adata)

tji.draw_page()

tji.draw_graph()

tji.louvain_cluster()

tji.draw_graph_paga()