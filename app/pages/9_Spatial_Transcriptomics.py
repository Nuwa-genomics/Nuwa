import streamlit as st
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import pandas as pd
import threading
import numpy as np
import os
import plotly.express as px
from utils.plotting import plot_ripley, plot_co_occurrence, plot_centrality_scores

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
    """
    View expression profiles while retaining spatial information.

    Notes
    -----
    .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/spatial_page.png
    """

    def __init__(self, adata):
        self.adata = adata

    def draw_page(self):
        st.title("Spatial Transcriptomics")
        plt.style.use('dark_background')
        self.col1, self.col2, self.col3 = st.columns(3, gap="medium")
        if 'spatial_plots' not in st.session_state:
            st.session_state["spatial_plots"] = dict(nhood_enrichment=None, interaction_matrix=None, ripley_score=None, co_occurance_score=None, centrality_score=None)
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
        nhood_enrichment, interaction_matrix, ripley_score, co_occurance_score, centrality_score = st.tabs(
            ['Neighbourhood enrichment', 'Interaction matrix', 'Ripley score', 'Co-occurance score', 'Centrality score']
        )

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

        with ripley_score:
            if st.session_state['spatial_plots']['ripley_score'] == None:
                st.info("You must run Ripley scoring first.")
            else:
                fig = st.session_state.spatial_plots["ripley_score"]
                fig_plt = plot_ripley(adata=self.adata, cluster_key=fig['cluster_key'], plot_sims=fig['plot_sims'], mode=fig['mode'], height=800)
                st.plotly_chart(fig_plt, use_container_width=True)

        with co_occurance_score:
            if st.session_state['spatial_plots']['co_occurance_score'] == None:
                st.info("You must run co-occurance scoring first.")
            else:
                fig = st.session_state.spatial_plots["co_occurance_score"]
                for fig_plt in plot_co_occurrence(self.adata, cluster_key=fig['cluster_key'], clusters=fig['clusters'], height=800):
                    st.plotly_chart(fig_plt, use_container_width=True)

        with centrality_score:
            if st.session_state.spatial_plots["centrality_score"] == None:
                st.info("You must run centrality scoring first.")
            else:
                fig = st.session_state.spatial_plots["centrality_score"]
                col1, col2, col3 = st.columns(3, gap="large")
                with col1:
                    fig_plt = plot_centrality_scores(self.adata, cluster_key=fig['cluster_key'], height=800)[0]
                    st.plotly_chart(fig_plt, use_container_width=True)
                with col2:
                    fig_plt = plot_centrality_scores(self.adata, cluster_key=fig['cluster_key'], height=800)[1]
                    st.plotly_chart(fig_plt, use_container_width=True)
                with col3:
                    fig_plt = plot_centrality_scores(self.adata, cluster_key=fig['cluster_key'], height=800)[2]
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
        import squidpy as sq

        sq.pl.spatial_scatter(adata, color="celltype_mapped_refined", shape=None, figsize=(10, 10))
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

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/nhood_enrichment1.png
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/nhood_enrichment2.png

        Example
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
                    cluster_key = st.selectbox(label="Cluster Key", options=self.adata.obs_keys())
                    col1, col2 = st.columns(2)
                    n_perms = col1.number_input(label="n_perms", min_value=1, value=1000)
                    mode = col2.selectbox(label="mode", options=['zscore', 'count'])
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

                        st.toast("Ran neighbourhood enrichment", icon="✅")
                except Exception as e:
                    st.error(e)

    def ripley_score(self):
        """
        In addition to the neighbor enrichment score, we can further investigate spatial organization of cell types in tissue by means of the Ripley’s statistics. Ripley’s statistics allow analyst to evaluate whether a discrete annotation (e.g. cell-type) appears to be clustered, dispersed or randomly distributed on the area of interest.
        
        Parameters
        ----------
        cluster_key: str
            Key by which to compute neighbouhood enrichment scores.

        max_dist: int
            Maximum distance between spatial clusters.

        n_neighbours: int
            Number of neighbours.

        n_simulations: int
            Number of simulations.

        mode: char
            Ripley mode (F, G or L).

        plot_sims: bool
            Plot simulations.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/ripley_score1.png
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/ripley_score2.png

        Example
        -------
        import squidpy as sq

        mode = "L"
        sq.gr.ripley(adata, cluster_key="cluster", mode=mode, max_dist=500)
        sq.pl.ripley(adata, cluster_key="cluster", mode=mode)
        """
        with self.col2:
            with st.form(key="ripley_score_form"):
                st.subheader("Ripley score")
                try:
                    col1, col2 = st.columns(2, gap="medium")
                    cluster_key = col1.selectbox(label="Cluster Key", options=self.adata.obs_keys(), key="sb_cluster_key_ripley")
                    max_dist = col2.number_input(label="max_dist", value=500, step=1)
                    n_neighbours = col1.number_input(label="n_neighbours", value=2, min_value=2, step=1)
                    n_simulations = col2.number_input(label="n_simulations", value=100)
                    mode = col1.radio(label="mode", options=['F', 'G', 'L'])
                    plot_sims = st.toggle(label="plot_sims", value=False)
                    subcol1, _, _, _ = st.columns(4)
                    submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
                    if submit_btn:
                        with st.spinner(text="Calculating Ripley score"):
                            sq.gr.ripley(self.adata, cluster_key=cluster_key, mode=mode, max_dist=max_dist, n_neigh=n_neighbours, n_simulations=n_simulations)
                            st.session_state.spatial_plots['ripley_score'] = dict(cluster_key=cluster_key, mode=mode, plot_sims=plot_sims)
                        #write to script state
                        st.session_state["script_state"].add_script(f"sq.gr.ripley(adata, cluster_key={st.session_state.sb_cluster_key_ripley}, mode='L', max_dist=500)")
                        st.session_state["script_state"].add_script(f"sq.pl.ripley(adata, cluster_key={st.session_state.sb_cluster_key_ripley}, mode='L')")
                        st.session_state["script_state"].add_script("plt.show()")

                        st.toast(f"Ran Ripley's {mode} scoring", icon="✅")
                except Exception as e:
                    st.error(e)

    def co_occurance_score(self):
        """
        Calculate co-occurance scoring between original spatial coordinates. It can be defined as $$ \\frac{p(exp|cond)}{p(exp)} $$.

        Parameters
        ----------
        cluster_key: str
            Key by which to compute neighbouhood enrichment scores.

        clusters: str
            Cluster annotation to use in co-occurance scoring.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/cooccurance1.png
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/cooccurance2.png

        Example
        -------
        import squidpy as sq

        sq.gr.co_occurrence(adata, cluster_key="celltype_mapped_refined")
        sq.pl.co_occurrence(adata, cluster_key="celltype_mapped_refined", clusters="Lateral plate mesoderm", figsize=(10, 5))
        """
        with self.col3:
            st.subheader("Co-occurance score")
            st.selectbox(label="Cluster Key", options=self.adata.obs_keys(), key="sb:spatial:co_occurance:cluster_key")
            with st.form(key="cooccurance_form"):
                
                try:
                    col1, col2 = st.columns(2)
                    clusters = st.multiselect(label="Clusters", options=self.adata.obs[f"{st.session_state['sb:spatial:co_occurance:cluster_key']}"].unique())
                   
                    subcol1, _, _, _ = st.columns(4)
                    submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
                    if submit_btn:
                        with st.spinner(text="Calculating co-occurance score"):
                            cluster_key = st.session_state['sb:spatial:co_occurance:cluster_key']
                            sq.gr.co_occurrence(self.adata, cluster_key=cluster_key)  
                            st.session_state.spatial_plots['co_occurance_score'] = dict(cluster_key=cluster_key, clusters=clusters)                
                            
                            #write to script state
                            # st.session_state["script_state"].add_script(f"sq.gr.co_occurrence(adata, cluster_key={st.session_state.sb_cluster_key_cooc})")
                            # st.session_state["script_state"].add_script(f"sq.pl.co_occurrence(adata, cluster_key={st.session_state.sb_cluster_key_cooc}, clusters={st.session_state.ms_clusters_cooc}")
                            # st.session_state["script_state"].add_script("plt.show()")

                        st.toast("Ran co-occurance scoring", icon="✅")
                except Exception as e:
                    st.error(e)


    def interaction_matrix(self):
        """
        Count the number of edges that each cluster share with all the others and plot in a heatmap.
        
        Parameters
        ----------
        cluster_key: str
            Key by which to compute neighbouhood enrichment scores.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/interaction_matrix1.png
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/interaction_matrix2.png
        
        Example
        -------
        import squidpy as sq

        sq.gr.spatial_neighbors(adata)
        sq.gr.interaction_matrix(adata, cluster_key="cell type")
        sq.pl.interaction_matrix(adata, cluster_key="cell type")
        """
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
                            sq.gr.spatial_neighbors(self.adata)
                            sq.gr.interaction_matrix(self.adata, cluster_key=cluster_key)
                            sq.pl.interaction_matrix(self.adata, cluster_key=cluster_key)
                            st.session_state.spatial_plots['interaction_matrix'] = dict(cluster_key=cluster_key)

                            #write to script state
                            # st.session_state["script_state"].add_script(f"sq.gr.interaction_matrix(adata, cluster_key={st.session_state.sb_cluster_key_inter_matrix})")
                            # st.session_state["script_state"].add_script(f"sq.pl.interaction_matrix(adata, cluster_key={st.session_state.sb_cluster_key_inter_matrix})")
                            # st.session_state["script_state"].add_script("plt.show()")

                        st.toast("Successfully ran interaction matrix", icon="✅")
                except Exception as e:
                    st.error(e)

    def centrality_score(self):
        """
        Compute degree, average and closeness centralities of the spatial graph.

        Parameters
        ----------
        cluster_key: str
            Key by which to compute neighbouhood enrichment scores.
        
        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/centrality_score1.png
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/centrality_score2.png

        Example
        -------
        import squidpy as sq

        sq.gr.spatial_neighbors(adata)
        sq.gr.centrality_scores(adata, cluster_key="cell type")
        sq.pl.centrality_scores(adata, cluster_key="cell type", figsize=(20, 5), s=500)
        """
        with self.col3:
            with st.form(key="centrality_score_form"):
                st.subheader("Centrality score")
                try:
                    cluster_key = st.selectbox(label="Cluster Key", options=self.adata.obs_keys())
                    
                    subcol1, _, _, _ = st.columns(4)
                    submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
                    if submit_btn:
                        with st.spinner(text="Calculating centrality score"):
                            sq.gr.spatial_neighbors(self.adata)
                            sq.gr.centrality_scores(self.adata, cluster_key=cluster_key)
                            
                            ax_cent = sq.pl.centrality_scores(self.adata, cluster_key=cluster_key, s=500)
                            st.session_state.spatial_plots['centrality_score'] = dict(cluster_key=cluster_key)
                            
                            #write to script state
                            # st.session_state["script_state"].add_script(f"sq.gr.centrality_scores(adata, cluster_key={st.session_state.sb_cluster_key_centrality_score})")
                            # st.session_state["script_state"].add_script(f"sq.pl.centrality_scores(adata, cluster_key={st.session_state.sb_cluster_key_centrality_score}, s=500)")
                            # st.session_state["script_state"].add_script("plt.show()")

                        st.toast("Successfully ran centrality scoring", icon="✅")
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

