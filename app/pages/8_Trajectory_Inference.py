import pickle
import scanpy as sc
import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time

from models.AdataModel import AdataModel
from components.sidebar import *
import os

st.set_page_config(layout="wide", page_title='Nuwa', page_icon='🧬')

os.chdir('/app')

with open('css/common.css') as f:
    common_style = f"""
                <style>
                {f.read()}
                </style>
                """
    st.markdown(common_style, unsafe_allow_html=True)

class Trajectory_Inference:
    def __init__(self, adata):
        try:
            self.adata = adata
            self.adata.raw = adata
            #convert to float64
            self.adata.X = self.adata.X.astype('float64')
            #compute pca
            sc.tl.pca(self.adata, svd_solver='arpack')
            #compute neighbours
            sc.pp.neighbors(self.adata, n_neighbors=4, n_pcs=20)
            sc.tl.draw_graph(self.adata)
            #compute louvain clusters
            sc.tl.louvain(self.adata, resolution=1.0)

        except Exception as e:
            st.error(e)
        

    def draw_graph(self):
        with self.col1:
            try:
                st.subheader("Cluster chart")
                with st.spinner(text="Plotting clusters"):
                    subcol1, _, _, _ = st.columns(4, gap='large')
                    key_index = self.adata.obs_keys().index('louvain') #created in constructor
                    if key_index == -1: #in case louvain was not created
                        key_index = 0
                    subcol1.selectbox(label='color', options=self.adata.obs.columns, key='sb_plt_colors', index=key_index)
                    st.slider(label="Point size", min_value=1, max_value=100, key="sl_traj_graph_size", value=50)
                    #denoise
                    #sc.tl.diffmap(self.adata)
                    #sc.pp.neighbors(self.adata, n_neighbors=10, use_rep='X_diffmap')
                    #sc.tl.draw_graph(self.adata)
                    self.ax_df = pd.DataFrame({'fr1': self.adata.obsm['X_draw_graph_fr'][:,0], 'fr2': self.adata.obsm['X_draw_graph_fr'][:,1], 'color': self.adata.obs[st.session_state.sb_plt_colors]})
                    st.scatter_chart(self.ax_df, x='fr1', y='fr2', color='color', height=600, size=st.session_state.sl_traj_graph_size)

            except Exception as e:
                st.error(e)

    def draw_graph_paga(self):
        with self.col2:
            try:
                with st.form(key="form_draw_graph_paga"):
                    st.subheader("Louvain with PAGA embedding")
                    plt.style.use('dark_background')
                    graph_paga_options = st.multiselect(label="Gene", options=np.append(self.adata.to_df().columns, 'louvain'), default=['louvain', self.adata.to_df().columns[0]], key="ms_louvain_colors_paga")
                    subcol1, _, _, _ = st.columns(4)
                    empty = st.container()
                    draw_paga_graph_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
                    if draw_paga_graph_btn:
                        with st.spinner(text="Recomputing with PAGA initialisation"):
                            sc.tl.draw_graph(self.adata, init_pos='paga')
                            ax = sc.pl.draw_graph(self.adata, color=graph_paga_options, legend_loc='on data')
                            empty.pyplot(ax)
                            #write to script
                            st.session_state["script_state"].add_script("sc.tl.draw_graph(adata, init_pos='paga')")
                            st.session_state["script_state"].add_script(f"sc.pl.draw_graph(adata, color={graph_paga_options}, legend_loc='on data')")
                            st.session_state["script_state"].add_script("plt.show()")
            
            except Exception as e:
                st.error(e)



    def louvain_cluster(self):
        with self.col2:
            try:
                with st.form(key="louvain_cluster_form"):
                    st.subheader("PAGA")
                    plt.style.use('classic')
                    options = st.multiselect(label="Gene", options=np.append(self.adata.to_df().columns, 'louvain'), default=['louvain', self.adata.to_df().columns[0]], key="ms_louvain_colors")
                    empty = st.container()
                    subcol1, _, _, _ = st.columns(4)
                    submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
                    
                    if submit_btn:
                        with st.spinner(text="Computing clusters"):
                            #louvain cluster
                            sc.tl.louvain(self.adata, resolution=1.0)

                            #paga
                            sc.tl.paga(self.adata, groups='louvain')
                            ax = sc.pl.paga(self.adata, colors=st.session_state.ms_louvain_colors, cmap='viridis', node_size_scale=3)
                            empty.pyplot(ax)
                            
                            #write to script state
                            st.session_state["script_state"].add_script("sc.tl.louvain(adata, resolution=1.0)")
                            st.session_state["script_state"].add_script("sc.tl.paga(adata, groups='louvain')")
                            st.session_state["script_state"].add_script(f"sc.pl.paga(adata, colors={options}, cmap='viridis', node_size_scale=3)")
                            st.session_state["script_state"].add_script("plt.show()")
                        

            except Exception as e:
                st.error(e)

    def draw_page(self):
        st.title("Trajectory Inference")
        #columns
        self.col1, self.col2 = st.columns(2, gap="medium")
        self.draw_graph()
        self.louvain_cluster()
        self.dpt()
        self.draw_graph_paga()
        #self.show_path()


    def dpt(self):
        with self.col2:
            try:
                with st.form(key="dpt_form"):
                    st.subheader("Diffusion Pseudotime")
                    plt.style.use('dark_background')
                    root_cell = st.selectbox(label='Root cell', options=(self.adata.obs['louvain'].unique()), key='sb_root_cell')
                    empty = st.container()
                    subcol1, _, _, _ = st.columns(4)
                    dpt_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
                    if dpt_btn:
                        with st.spinner(text="Computing dpt"):
                            self.adata.uns['iroot'] = np.flatnonzero(self.adata.obs['louvain']  == st.session_state.sb_root_cell)[0]
                            
                            sc.tl.dpt(self.adata)
                            sc.pp.log1p(self.adata)
                            sc.pp.scale(self.adata)

                            ax1 = sc.pl.draw_graph(self.adata, color=['louvain', 'dpt_pseudotime'], legend_loc='on data', cmap='viridis')
                            empty.pyplot(ax1)
                            
                            #add to script state
                            
                            st.session_state["script_state"].add_script("adata.uns['iroot'] = np.flatnonzero(adata.obs['louvain']  == st.session_state.sb_root_cell)[0]")
                            st.session_state["script_state"].add_script("sc.tl.dpt(adata)")
                            st.session_state["script_state"].add_script("sc.pp.log1p(adata)")
                            st.session_state["script_state"].add_script("sc.pp.scale(adata)")
                            st.session_state["script_state"].add_script("sc.pl.draw_graph(adata, color=['louvain', 'dpt_pseudotime'], legend_loc='on data', cmap='viridis')")
                            st.session_state["script_state"].add_script("plt.show()")

                            #     sc.tl.dpt(adata_raw, n_branchings=1, n_dcs=10)
                            #     adata_raw.obs.dpt_groups = adata_raw.obs.dpt_groups[:60]
                            #     ax2 = sc.pl.dpt_timeseries(adata_raw)
                            #     st.pyplot(ax2)

            except Exception as e:
                st.error(e)

    def show_path(self):
        with self.col2:
            try:
                st.subheader("Paths")


                paths = [('erythrocytes', [16, 12, 7, 13, 18, 6, 5, 10]),
                        ('neutrophils', [16, 0, 4, 2, 14, 19]),
                        ('monocytes', [16, 0, 4, 11, 1, 9, 24])]


                self.adata.obs['distance'] = self.adata.obs['dpt_pseudotime']
                self.adata.obs['clusters'] = self.adata.obs['louvain']  # just a cosmetic change
                self.adata.uns['clusters_colors'] = self.adata.uns['louvain_colors']

                gene_names = ['Gata2', 'Gata1', 'Klf1', 'Epor', 'Hba-a2',  # erythroid
                'Elane', 'Cebpe', 'Gfi1',                    # neutrophil
                'Irf8', 'Csf1r', 'Ctsg'] 

        
                _, axs = plt.subplots(ncols=3, figsize=(6, 2.5), gridspec_kw={'wspace': 0.05, 'left': 0.12})
                plt.subplots_adjust(left=0.05, right=0.98, top=0.82, bottom=0.2)


                for ipath, (descr, path) in enumerate(paths):
                    ax = sc.pl.paga_path(
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
                    #data.to_csv('./write/paga_path_{}.csv'.format(descr))
                
                st.pyplot(axs)
                plt.savefig('./figures/paga_path_paul15.pdf')

            except Exception as e:
                st.error(e)
                

try:
    adata = st.session_state.adata_state.current.adata
    sidebar = Sidebar()

    sidebar.show()

    tji = Trajectory_Inference(adata)

    tji.draw_page()

    sidebar.show_preview()
    sidebar.export_script()
    sidebar.delete_experiment_btn()
    sidebar.show_version()

except KeyError as ke:
    print("KeyError: ", ke)
    st.error("Couldn't find adata object in session, have you uploaded one?")





