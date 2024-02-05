import streamlit as st
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import pandas as pd
import threading
import numpy as np
import os
import plotly.express as px

from models.AdataModel import AdataModel
from components.sidebar import *
from utils.AdataState import AdataState

st.set_page_config(layout="wide", page_title='Nuwa', page_icon='🧬')

os.chdir('/app')

with open('css/common.css') as f:
    common_style = f"""
                <style>
                {f.read()}
                </style>
                """
    st.markdown(common_style, unsafe_allow_html=True)

class Spatial_transcriptomics:
    def __init__(self, adata):
        self.adata = adata

    def draw_page(self):
        st.title("Spatial Transcriptomics")
        plt.style.use('dark_background')
        self.col1, self.col2, self.col3 = st.columns(3, gap="medium")
        if 'spatial_plots' not in st.session_state:
            st.session_state["spatial_plots"] = dict(nhood_enrichment=None, interaction_matrix=None)
        self.spatial_scatter()
        self.neighbourhood_enrichment()
        self.interaction_matrix()
        self.ripley_score()
        self.co_occurance_score()
        self.centrality_score()
        self.render_plots()
        #self.ligand_receptor_interaction()

    def render_plots(self):
        st.subheader("Plots")
        nhood_enrichment, interaction_matrix = st.tabs(['Neighbourhood enrichment', 'Interaction matrix'])
        with nhood_enrichment:
            if st.session_state['spatial_plots']["nhood_enrichment"] == None:
                st.info("You must run neighbourhood enrichment first.")
            else:
                fig = st.session_state.spatial_plots["nhood_enrichment"]
                fig_plt = px.imshow(
                    self.adata.uns[f"{fig['cluster_key']}_nhood_enrichment"][fig['mode']], 
                    labels=dict(x=fig['cluster_key'], y=fig['cluster_key'], color=fig['mode']), 
                    title="Neighbourhood enrichment", aspect="auto", height=800,
                    x=adata.obs[fig['cluster_key']].unique(), y=adata.obs[fig['cluster_key']].unique()
                )
                st.plotly_chart(fig_plt, use_container_width=True)

        with interaction_matrix:
            if st.session_state['spatial_plots']['interaction_matrix'] == None:
                st.info("You must run interaction matrix first.")
            else:
                fig = st.session_state.spatial_plots["interaction_matrix"]
                fig_plt = px.imshow(
                    self.adata.uns[f"{fig['cluster_key']}_interactions"], 
                    labels=dict(x=fig['cluster_key'], y=fig['cluster_key'], color='Score'), 
                    title="Neighbourhood enrichment", aspect="auto", height=800,
                    x=adata.obs[fig['cluster_key']].unique(), y=adata.obs[fig['cluster_key']].unique()
                )
                st.plotly_chart(fig_plt, use_container_width=True)
        

    def spatial_scatter(self):
        """
        Plot a spatial graph. Optionally overlay on histology image if available.

        Parameters
        ----------
        color: str
            Obs value to color samples.

        point_size: int
            Size of markers in scatter graph.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/spatial_plot.png
        
        Example
        -------
        import scanpy as sc

        sc.pl.spatial(adata, color="cell_type")
        """
        with st.spinner(text="Loading graph"):
            with self.col1:
                st.subheader("Spatial plot")
                try:
                    plt.style.use('seaborn-v0_8-pastel')
                    st.selectbox(label="Colour", options=reversed(self.adata.obs.columns), key="sb_spatial_colours")
                    ax_df = pd.DataFrame({'spatial1': self.adata.obsm['spatial'][:,0], 'spatial2': self.adata.obsm['spatial'][:,1], 'color': self.adata.obs[st.session_state.sb_spatial_colours]})
                    st.slider(label="Point size", min_value=1, max_value=100, value=10, key="sl_point_size")
                    st.scatter_chart(ax_df, x='spatial1', y='spatial2', color='color', height=480, size=st.session_state.sl_point_size)
                    
                except Exception as e:
                    st.error(e)

    def neighbourhood_enrichment(self):
        """
        Plots a matrix plot using an enrichment score on spatial proximity of clusters. If spots belonging to two different clusters are often close to each other, then they will have a high score and can be defined as being enriched. On the other hand, if they are far apart, the score will be low and they can be defined as depleted.
        
        Parameters
        ----------
        cluster_key: str
            Key by which to compute neighbouhood enrichment scores.

        n_perms: int
            Number of permutations when computing the score.

        Exmaple
        -------
        import squidpy as sq

        sq.gr.spatial_neighbors(adata)
        sq.gr.nhood_enrichment(adata, n_perms=1000, cluster_key="cell_type")
        sq.pl.nhood_enrichment(adata, cluster_key="cell_type")
        """
        with self.col2:
            with st.form(key="n_enrich_form"):
                st.subheader("Neighbourhood Enrichment")
                try:
                    col1, col2 = st.columns(2, gap="large")
                    cluster_key = col1.selectbox(label="Cluster Key", options=self.adata.obs_keys())
                    n_perms = col2.number_input(label="n_perms", min_value=1, value=1000)
                    mode = col1.selectbox(label="mode", options=['zscore', 'count'])
                    subcol1, _, _, _ = st.columns(4)
                    submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
                    if submit_btn:
                        with st.spinner(text="Running neighbourhood enrichment"):
                            sq.gr.spatial_neighbors(self.adata)
                            sq.gr.nhood_enrichment(self.adata, n_perms=n_perms, cluster_key=cluster_key)
                            
                            st.session_state.spatial_plots['nhood_enrichment'] = dict(cluster_key=cluster_key, mode=mode)
                            #write to script state
                            # st.session_state["script_state"].add_script("sq.gr.spatial_neighbors(adata)")
                            # st.session_state["script_state"].add_script(f"sq.gr.nhood_enrichment(adata, cluster_key={st.session_state.sb_cluster_key_n_enrich})")
                            # st.session_state["script_state"].add_script("")
                            # st.session_state["script_state"].add_script("")
                except Exception as e:
                    st.error(e)

    def ripley_score(self):
        with self.col2:
            with st.form(key="ripley_score_form"):
                st.subheader("Ripley score")
                try:
                    st.selectbox(label="Cluster Key", options=reversed(self.adata.obs.columns), key="sb_cluster_key_ripley")
                    empty = st.empty()
                    subcol1, _, _, _ = st.columns(4)
                    submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
                    if submit_btn:
                        with st.spinner(text="Calculating Ripley score"):
                            sq.gr.ripley(self.adata, cluster_key=st.session_state.sb_cluster_key_ripley, mode="L", max_dist=500)
                            ax_ripley = sq.pl.ripley(self.adata, cluster_key=st.session_state.sb_cluster_key_ripley, mode="L")
                            empty.pyplot(ax_ripley)
                        #write to script state
                        st.session_state["script_state"].add_script(f"sq.gr.ripley(adata, cluster_key={st.session_state.sb_cluster_key_ripley}, mode='L', max_dist=500)")
                        st.session_state["script_state"].add_script(f"sq.pl.ripley(adata, cluster_key={st.session_state.sb_cluster_key_ripley}, mode='L')")
                        st.session_state["script_state"].add_script("plt.show()")
                except Exception as e:
                    st.error(e)

    def co_occurance_score(self):
        with self.col3:
            with st.form(key="cooccurance_form"):
                st.subheader("Co-occurance score")
                try:
                    st.selectbox(label="Cluster Key", options=reversed(self.adata.obs.columns.unique()), key="sb_cluster_key_cooc")
                    options = self.adata.obs[st.session_state.sb_cluster_key_cooc].unique()
                    st.multiselect(label="Clusters", options=options, key="ms_clusters_cooc", default=options[0])
                    empty = st.empty()
                    subcol1, _, _, _ = st.columns(4)
                    submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
                    if submit_btn:
                        with st.spinner(text="Calculating co-occurance score"):
                            sq.gr.co_occurrence(self.adata, cluster_key=st.session_state.sb_cluster_key_cooc)
                            ax_cooc = sq.pl.co_occurrence(self.adata, cluster_key=st.session_state.sb_cluster_key_cooc, clusters=st.session_state.ms_clusters_cooc)
                            empty.empty()
                            empty.pyplot(ax_cooc)
                            #write to script state
                            st.session_state["script_state"].add_script(f"sq.gr.co_occurrence(adata, cluster_key={st.session_state.sb_cluster_key_cooc})")
                            st.session_state["script_state"].add_script(f"sq.pl.co_occurrence(adata, cluster_key={st.session_state.sb_cluster_key_cooc}, clusters={st.session_state.ms_clusters_cooc}")
                            st.session_state["script_state"].add_script("plt.show()")
                except Exception as e:
                    st.error(e)


    def interaction_matrix(self):
        with self.col3:
            with st.form(key="interaction_matrix_form"):
                st.subheader("Interaction matrix")
                try:
                    col1, col2 = st.columns(2, gap="medium")
                    cluster_key = col1.selectbox(label="Cluster Key", options=self.adata.obs_keys())
                    subcol1, _, _, _ = st.columns(4)
                    submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
                    if submit_btn:
                        with st.spinner(text="Computing interaction matrix"):
                            sq.gr.interaction_matrix(self.adata, cluster_key=cluster_key)
                            sq.pl.interaction_matrix(self.adata, cluster_key=cluster_key)
                            st.session_state.spatial_plots['interaction_matrix'] = dict(cluster_key=cluster_key)

                            #write to script state
                            # st.session_state["script_state"].add_script(f"sq.gr.interaction_matrix(adata, cluster_key={st.session_state.sb_cluster_key_inter_matrix})")
                            # st.session_state["script_state"].add_script(f"sq.pl.interaction_matrix(adata, cluster_key={st.session_state.sb_cluster_key_inter_matrix})")
                            # st.session_state["script_state"].add_script("plt.show()")
                except Exception as e:
                    st.error(e)

    def centrality_score(self):
        with self.col3:
            with st.form(key="centrality_score_form"):
                st.subheader("Centrality score")
                try:
                    st.selectbox(label="Cluster Key", options=reversed(self.adata.obs.columns), key="sb_cluster_key_centrality_score")
                    empty = st.empty()
                    subcol1, _, _, _ = st.columns(4)
                    submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
                    if submit_btn:
                        with st.spinner(text="Calculating centrality score"):
                            sq.gr.centrality_scores(self.adata, cluster_key=st.session_state.sb_cluster_key_centrality_score)
                            ax_cent = sq.pl.centrality_scores(self.adata, cluster_key=st.session_state.sb_cluster_key_centrality_score, s=500)
                            empty.empty()
                            empty.pyplot(ax_cent)
                            #write to script state
                            st.session_state["script_state"].add_script(f"sq.gr.centrality_scores(adata, cluster_key={st.session_state.sb_cluster_key_centrality_score})")
                            st.session_state["script_state"].add_script(f"sq.pl.centrality_scores(adata, cluster_key={st.session_state.sb_cluster_key_centrality_score}, s=500)")
                            st.session_state["script_state"].add_script("plt.show()")
                except Exception as e:
                    st.error(e)

    def ligand_receptor_interaction(self):
        with self.col2:
            with st.form(key="ligand_receptor_form"):
                st.subheader("Ligand receptor interaction")
                try:
                    st.selectbox(label="Cluster Key", options=reversed(self.adata.obs.columns), key="sb_cluster_key_lri")
                    options = self.adata.obs[st.session_state.sb_cluster_key_lri].unique()
                    st.multiselect(label="Source groups", options=options, default=options[0], key="ms_lri_source_groups")
                    st.multiselect(label="Target groups", options=options, default=options[0], key="ms_lri_target_groups")
                    empty = st.empty()
                    subcol1, _, _, _ = st.columns(4)
                    submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
                    if submit_btn:
                        with st.spinner(text="Computing ligand-receptor interaction matrix"):
                            sq.gr.ligrec(self.adata, n_perms=100, cluster_key=st.session_state.sb_cluster_key_lri)
                            ax_lri = sq.pl.ligrec(self.adata, cluster_key=st.session_state.sb_cluster_key_lri, 
                            source_groups=st.session_state.ms_lri_source_groups, target_groups=st.session_state.ms_lri_target_groups,
                            means_range=(0.3, np.inf), alpha=1e-4, swap_axes=True)
                            empty.empty()
                            empty.pyplot(ax_lri)
                            #write to script state
                            st.session_state["script_state"].add_script(f"sq.gr.ligrec(adata, n_perms=100, cluster_key={st.session_state.sb_cluster_key_lri})")
                            st.session_state["script_state"].add_adata(f"sq.pl.ligrec(adata, cluster_key={st.session_state.sb_cluster_key_lri},  \
                                source_groups={st.session_state.ms_lri_source_groups}, target_groups={st.session_state.ms_lri_target_groups}, \
                                means_range=(0.3, np.inf), alpha=1e-4, swap_axes=True)")
                except Exception as e:
                    st.error(e)
                

try:
    adata = st.session_state.adata_state.current.adata
    sidebar = Sidebar()
    sidebar.show()

    spatial_t = Spatial_transcriptomics(adata)
    spatial_t.draw_page()
    sidebar.show_preview()
    sidebar.export_script()
    sidebar.delete_experiment_btn()
    sidebar.show_version()

except KeyError as ke:
    st.error("KeyError: ", ke)

