import pickle
import scanpy as sc
import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time

from models.AdataModel import AdataModel
from components.sidebar import *

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
                st.slider(label="Point size", min_value=1, max_value=100, key="sl_traj_graph_size", value=50)
                #denoise
                # sc.tl.diffmap(self.adata)
                # sc.pp.neighbors(self.adata, n_neighbors=10, use_rep='X_diffmap')
                # sc.tl.draw_graph(self.adata)
                self.ax_df = pd.DataFrame({'fr1': self.adata.obsm['X_draw_graph_fr'][:,0], 'fr2': self.adata.obsm['X_draw_graph_fr'][:,1], 'color': self.adata.obs[st.session_state.sb_plt_colors]})
                st.scatter_chart(self.ax_df, x='fr1', y='fr2', color='color', height=600, size=st.session_state.sl_traj_graph_size)

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
        self.draw_graph()
        self.louvain_cluster()
        self.draw_graph_paga()
        self.dpt()
        #self.show_path()


    def dpt(self):
        with self.col2:
            st.subheader("Diffusion Pseudotime")
            st.selectbox(label='Root cell', options=(self.adata.obs['louvain'].unique()), key='sb_root_cell')
            with st.expander(label="Show figure", expanded=True):
                self.adata.uns['iroot'] = np.flatnonzero(self.adata.obs['louvain']  == st.session_state.sb_root_cell)[0]
                
                adata_raw = self.adata.copy()
                sc.tl.dpt(adata_raw)
                sc.pp.log1p(adata_raw)
                sc.pp.scale(adata_raw)

                ax = sc.pl.draw_graph(adata_raw, color=['louvain', 'dpt_pseudotime'], legend_loc='on data', cmap='viridis')
                st.pyplot(ax)

    def show_path(self):
        with self.col2:
            st.subheader("Paths")

            paths = [('erythrocytes', [16, 12, 7, 13, 18, 6, 5, 10]),
                    ('neutrophils', [16, 0, 4, 2, 14, 19]),
                    ('monocytes', [16, 0, 4, 11, 1, 9, 24])]
            
            self.adata.obs['distance'] = self.adata.obs['dpt_pseudotime']
            self.adata.obs['clusters'] = self.adata.obs['louvain']  # just a cosmetic change
            self.adata.uns['clusters_colors'] = self.adata.uns['louvain_colors']

            _, axs = plt.subplots(ncols=3, figsize=(6, 2.5), gridspec_kw={'wspace': 0.05, 'left': 0.12})
            plt.subplots_adjust(left=0.05, right=0.98, top=0.82, bottom=0.2)

            gene_names = ['Gata2', 'Gata1', 'Klf1', 'Epor', 'Hba-a2',  # erythroid
              'Elane', 'Cebpe', 'Gfi1',                    # neutrophil
              'Irf8', 'Csf1r', 'Ctsg'] 

      

            for ipath, (descr, path) in enumerate(paths):
                
                _, data = sc.pl.paga_path(
                    self.adata, path, gene_names,                         
                    show_node_names=False,
                    ax=axs[ipath],
                    ytick_fontsize=12,
                    left_margin=0.15,
                    n_avg=50,
                    annotations=['distance'],
                    show_yticks=True if ipath==0 else False,
                    show_colorbar=False,
                    color_map='Greys',
                    groups_key='clusters',
                    color_maps_annotations={'distance': 'viridis'},
                    title='{} path'.format(descr),
                    return_data=True,
                    show=False)
                
            st.pyplot(axs)
                


try:
    adata_model: AdataModel = st.session_state["adata"]
    show_sidebar(adata_model)
except KeyError as ke:
    print('Key Not Found in Employee Dictionary:', ke)

adata_bytes = get_adata(adataList=adata_model, name=st.session_state.sb_adata_selection).adata
st.session_state["current_adata"] = pickle.loads(adata_bytes)

tji = Trajectory_Inference(adata=st.session_state.current_adata)

tji.draw_page()

show_preview()

