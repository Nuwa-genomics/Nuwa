import streamlit as st
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import pandas as pd

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
        self.col1, self.col2, self.col3 = st.columns(3, gap="large")
        self.spatial_scatter()
        self.neighbourhood_enrichment()
        self.interaction_matrix()
        self.ripley_score()
        self.co_occurance_score()
        self.centrality_score()

    def spatial_scatter(self):
        with st.spinner(text="Loading graph"):
            with self.col1:
                st.subheader("Spatial plot")
                st.selectbox(label="Colour", options=reversed(self.adata.obs.columns), key="sb_spatial_colours")
                ax_df = pd.DataFrame({'spatial1': self.adata.obsm['spatial'][:,0], 'spatial2': self.adata.obsm['spatial'][:,1], 'color': self.adata.obs[st.session_state.sb_spatial_colours]})
                st.scatter_chart(ax_df, x='spatial1', y='spatial2', color='color', height=480, size=10)

    def neighbourhood_enrichment(self):
        with self.col2:
            st.subheader("Neighbourhood Enrichment")
            st.selectbox(label="Cluster Key", options=reversed(self.adata.obs.columns), key="sb_cluster_key_n_enrich")
            with st.spinner(text="Running neighbourhood enrichment"):
                sq.gr.spatial_neighbors(self.adata)
                sq.gr.nhood_enrichment(self.adata, cluster_key=st.session_state.sb_cluster_key_n_enrich)
                ax_n_enrich = sq.pl.nhood_enrichment(self.adata, cluster_key=st.session_state.sb_cluster_key_n_enrich)
                with st.expander(label="Show figure"):
                    st.pyplot(ax_n_enrich)
            st.divider()

    def ripley_score(self):
        with self.col2:
            st.subheader("Ripley score")
            st.selectbox(label="Cluster Key", options=reversed(self.adata.obs.columns), key="sb_cluster_key_ripley")
            with st.spinner(text="Calculating Ripley score"):
                sq.gr.ripley(self.adata, cluster_key=st.session_state.sb_cluster_key_ripley, mode="L", max_dist=500)
                ax_ripley = sq.pl.ripley(self.adata, cluster_key=st.session_state.sb_cluster_key_ripley, mode="L")
                with st.expander(label="Show Figure"):
                    st.pyplot(ax_ripley)
            st.divider()

    def co_occurance_score(self):
        with self.col3:
            st.subheader("Co-occurance score")
            st.selectbox(label="Cluster Key", options=reversed(self.adata.obs.columns.unique()), key="sb_cluster_key_cooc")
            options = self.adata.obs[st.session_state.sb_cluster_key_cooc].unique()
            st.multiselect(label="Clusters", options=options, key="ms_clusters_cooc", default=options[0])
            with st.spinner(text="Calculating co-occurance score"):
                sq.gr.co_occurrence(self.adata, cluster_key=st.session_state.sb_cluster_key_cooc)
                ax_cooc = sq.pl.co_occurrence(self.adata, cluster_key=st.session_state.sb_cluster_key_cooc, clusters=st.session_state.ms_clusters_cooc)
                with st.expander(label="Show Figure"):
                    st.pyplot(ax_cooc)


    def interaction_matrix(self):
        with self.col3:
            st.subheader("Interaction matrix")
            st.selectbox(label="Cluster Key", options=reversed(self.adata.obs.columns), key="sb_cluster_key_inter_matrix")
            with st.spinner(text="Computing interaction matrix"):
                sq.gr.interaction_matrix(self.adata, cluster_key=st.session_state.sb_cluster_key_inter_matrix)
                ax_inter_matrix = sq.pl.interaction_matrix(self.adata, cluster_key=st.session_state.sb_cluster_key_inter_matrix)
                with st.expander(label="Show figure"):
                    st.pyplot(ax_inter_matrix)
            st.divider()

    def centrality_score(self):
        with self.col2:
            st.subheader("Centrality score")
            st.selectbox(label="Cluster Key", options=reversed(self.adata.obs.columns), key="sb_cluster_key_centrality_score")
            with st.spinner(text="Calculating centrality score"):
                sq.gr.centrality_scores(self.adata, cluster_key=st.session_state.sb_cluster_key_centrality_score)
                ax_cent = sq.pl.centrality_scores(self.adata, cluster_key=st.session_state.sb_cluster_key_centrality_score, s=500)
                with st.expander(label="Show figure"):
                    st.pyplot(ax_cent)




adata = st.session_state.adata

spatial_t = Spatial_Transcriptomics(adata)

spatial_t.draw_page()

