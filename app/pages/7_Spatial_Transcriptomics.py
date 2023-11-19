import pickle
import streamlit as st
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import pandas as pd
import threading
import numpy as np

from models.AdataModel import AdataModel
from components.sidebar import *
from utils.AdataState import AdataState

st.set_page_config(layout="wide")

with open('css/common.css') as f:
    common_style = f"""
                <style>
                {f.read()}
                </style>
                """
    st.markdown(common_style, unsafe_allow_html=True)

class Spatial_Transcriptomics:
    def __init__(self, adata):
        self.adata = adata

    def draw_page(self):
        st.title("Spatial Transcriptomics")
        plt.style.use('dark_background')
        self.col1, self.col2, self.col3 = st.columns(3, gap="large")
        self.spatial_scatter()
        self.neighbourhood_enrichment()
        self.interaction_matrix()
        self.ripley_score()
        self.co_occurance_score()
        self.centrality_score()
        self.ligand_receptor_interaction()
        

    def spatial_scatter(self):
        with st.spinner(text="Loading graph"):
            with self.col1:
                st.subheader("Spatial plot")
                try:
                    plt.style.use('seaborn-v0_8-pastel')
                    st.selectbox(label="Colour", options=reversed(self.adata.obs.columns), key="sb_spatial_colours")
                    ax_df = pd.DataFrame({'spatial1': self.adata.obsm['spatial'][:,0], 'spatial2': self.adata.obsm['spatial'][:,1], 'color': self.adata.obs[st.session_state.sb_spatial_colours]})
                    st.slider(label="Point size", min_value=1, max_value=100, value=10, key="sl_point_size")
                    tab_plt, tab_img, tab_both = st.tabs(['plt', 'Image', 'both'])
                    with tab_plt:
                        st.scatter_chart(ax_df, x='spatial1', y='spatial2', color='color', height=480, size=st.session_state.sl_point_size)
                    with tab_img:
                        fig, ax = plt.subplots(nrows=1, ncols=1)
                        ax_sp_img: plt.Axes = sq.pl.spatial_scatter(self.adata, color="cluster", size=0, legend_fontsize=0.0, legend_fontweight="regular", colorbar=False, ax=ax)
                        st.pyplot(ax_sp_img)
                    with tab_both:
                        ax_sp = sq.pl.spatial_scatter(self.adata, color="cluster")
                        st.pyplot(ax_sp)
                    
                except Exception as e:
                    st.error(e)

    def neighbourhood_enrichment(self):
        with self.col2:
            st.subheader("Neighbourhood Enrichment")
            try:
                st.selectbox(label="Cluster Key", options=reversed(self.adata.obs.columns), key="sb_cluster_key_n_enrich")
                with st.spinner(text="Running neighbourhood enrichment"):
                    sq.gr.spatial_neighbors(self.adata)
                    sq.gr.nhood_enrichment(self.adata, cluster_key=st.session_state.sb_cluster_key_n_enrich)
                    ax_n_enrich = sq.pl.nhood_enrichment(self.adata, cluster_key=st.session_state.sb_cluster_key_n_enrich)
                    with st.expander(label="Show figure"):
                        st.pyplot(ax_n_enrich)
            except Exception as e:
                st.error(e)
            st.divider()

    def ripley_score(self):
        with self.col2:
            st.subheader("Ripley score")
            try:
                st.selectbox(label="Cluster Key", options=reversed(self.adata.obs.columns), key="sb_cluster_key_ripley")
                with st.spinner(text="Calculating Ripley score"):
                    sq.gr.ripley(self.adata, cluster_key=st.session_state.sb_cluster_key_ripley, mode="L", max_dist=500)
                    ax_ripley = sq.pl.ripley(self.adata, cluster_key=st.session_state.sb_cluster_key_ripley, mode="L")
                    with st.expander(label="Show Figure"):
                        st.pyplot(ax_ripley)
            except Exception as e:
                st.error(e)
            st.divider()

    def co_occurance_score(self):
        with self.col3:
            st.subheader("Co-occurance score")
            try:
                st.selectbox(label="Cluster Key", options=reversed(self.adata.obs.columns.unique()), key="sb_cluster_key_cooc")
                options = self.adata.obs[st.session_state.sb_cluster_key_cooc].unique()
                st.multiselect(label="Clusters", options=options, key="ms_clusters_cooc", default=options[0])
                with st.spinner(text="Calculating co-occurance score"):
                    sq.gr.co_occurrence(self.adata, cluster_key=st.session_state.sb_cluster_key_cooc)
                    ax_cooc = sq.pl.co_occurrence(self.adata, cluster_key=st.session_state.sb_cluster_key_cooc, clusters=st.session_state.ms_clusters_cooc)
                    with st.expander(label="Show Figure"):
                        st.pyplot(ax_cooc)
            except Exception as e:
                st.error(e)


    def interaction_matrix(self):
        with self.col3:
            st.subheader("Interaction matrix")
            try:
                st.selectbox(label="Cluster Key", options=reversed(self.adata.obs.columns), key="sb_cluster_key_inter_matrix")
                with st.spinner(text="Computing interaction matrix"):
                    sq.gr.interaction_matrix(self.adata, cluster_key=st.session_state.sb_cluster_key_inter_matrix)
                    ax_inter_matrix = sq.pl.interaction_matrix(self.adata, cluster_key=st.session_state.sb_cluster_key_inter_matrix)
                    with st.expander(label="Show figure"):
                        st.pyplot(ax_inter_matrix)
            except Exception as e:
                st.error(e)
            st.divider()

    def centrality_score(self):
        with self.col2:
            st.subheader("Centrality score")
            try:
                st.selectbox(label="Cluster Key", options=reversed(self.adata.obs.columns), key="sb_cluster_key_centrality_score")
                with st.spinner(text="Calculating centrality score"):
                    sq.gr.centrality_scores(self.adata, cluster_key=st.session_state.sb_cluster_key_centrality_score)
                    ax_cent = sq.pl.centrality_scores(self.adata, cluster_key=st.session_state.sb_cluster_key_centrality_score, s=500)
                    with st.expander(label="Show figure"):
                        st.pyplot(ax_cent)
            except Exception as e:
                st.error(e)

    def ligand_receptor_interaction(self):
        with self.col3:
            st.subheader("Ligand receptor interaction")
            try:
                st.selectbox(label="Cluster Key", options=reversed(self.adata.obs.columns), key="sb_cluster_key_lri")
                options = self.adata.obs[st.session_state.sb_cluster_key_lri].unique()
                st.multiselect(label="Source groups", options=options, default=options[0], key="ms_lri_source_groups")
                st.multiselect(label="Target groups", options=options, default=options[0], key="ms_lri_target_groups")
                with st.spinner(text="Computing ligand-receptor interaction matrix"):
                    sq.gr.ligrec(self.adata, n_perms=100, cluster_key=st.session_state.sb_cluster_key_lri)
                    ax_lri = sq.pl.ligrec(self.adata, cluster_key=st.session_state.sb_cluster_key_lri, 
                    source_groups=st.session_state.ms_lri_source_groups, target_groups=st.session_state.ms_lri_target_groups,
                    means_range=(0.3, np.inf), alpha=1e-4, swap_axes=True)
                    with st.expander(label="Show figure"):
                        st.pyplot(ax_lri)
            except Exception as e:
                st.error(e)
                

try:
    adata_state = AdataState(workspace_id=st.session_state.current_workspace.id)
    sidebar = Sidebar()
    sidebar.show()

    spatial_t = Spatial_Transcriptomics(adata_state.current.adata)
    spatial_t.draw_page()
    sidebar.show_preview()
    sidebar.delete_experiment_btn()

except KeyError as ke:
    print("KeyError: ", ke)
    st.error("Couldn't find adata object in session, have you uploaded one?")

