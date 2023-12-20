import streamlit as st
import scanpy as sc
from anndata import AnnData
from components.sidebar import *
import os
import numpy as np
import pandas as pd
import scanpy as sc
import gseapy
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

st.set_page_config(layout="wide", page_title='Nuwa', page_icon='🧬')

with open('css/common.css') as f:
    common_style = f"""
                <style>
                {f.read()}
                </style>
                """
    st.markdown(common_style, unsafe_allow_html=True)


class DGE:
    def __init__(self, adata):
        self.adata = adata
        self.draw_page()

    def save_adata(self):
        sc.write(filename=os.path.join(os.getenv('WORKDIR'), 'adata', st.session_state.adata_state.current.adata_name), adata=self.adata)
        st.session_state.adata_state.current.adata = self.adata

    def draw_page(self):
        col1, col2, col3 = st.columns(3)
        with col1:
            self.add_embeddings()
        with col2:
            self.stat_tests()
            
        self.visualize()

    def stat_tests(self):
        try:
            with st.form(key="stat_tests_form"):
                st.subheader("Run statistical tests")
                method = st.selectbox(label="Method", options=(['logreg', 't-test', 't-test_overestim_var', 'wilcoxon']), key='rb_method')
                group_by = st.selectbox(label="Group by", options=st.session_state.adata_state.current.adata.obs_keys())
                marker_genes_container = st.empty()
                subcol1, _, _, _ = st.columns(4)
                submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
                if submit_btn:
                    marker_genes_container.empty()
                    with st.spinner(text="Computing tests"):
                        sc.tl.rank_genes_groups(self.adata, groupby=group_by, method = str.lower(method), key_added=method)
                        #sc.pl.rank_genes_groups(st.session_state.adata_state.current.adata, n_genes=25, sharey=False, save=True)
                        self.save_adata()

        except Exception as e:
            st.error(e)

    def add_embeddings(self):
        try:
            with st.form(key="add_embeddings_form"):
                st.subheader("Cluster embeddings", help="Embed clusters in dataset for DGE.")
                algorithm = st.selectbox(label="Algorithm", options=['Leiden', 'Louvain'])
                resolution = st.number_input(label="Resolution", min_value=0.1, value=0.6, format="%.1f")
                subcol1, _, _ = st.columns(3)
                submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
                if submit_btn:
                    with st.spinner(text=f"Computing {algorithm} clusters"):
                        sc.pp.neighbors(self.adata, n_neighbors=10, n_pcs=40)
                        if algorithm == "Leiden":
                            sc.tl.leiden(self.adata, resolution=resolution, key_added=f"leiden_{resolution}") 
                        if algorithm == "Louvain":
                            sc.tl.louvain(self.adata, resolution=resolution, key_added=f"louvain_{resolution}")
                    
                    self.save_adata()

        except Exception as e:
            st.error(e)




    def compare_tests(self):
        try:
            with st.form(key="Compare tests"):
                st.subheader("Compare tests")
                compare_tests_container = st.empty()
                subcol1, _, _ = st.columns(3)
                compare_tests_submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
                if compare_tests_submit_btn:
                    wc = sc.get.rank_genes_groups_df(self.adata, group='0', key='wilcoxon', pval_cutoff=0.01, log2fc_min=0)['names']
                    tt = sc.get.rank_genes_groups_df(self.adata, group='0', key='t-test', pval_cutoff=0.01, log2fc_min=0)['names']
                    tt_ov = sc.get.rank_genes_groups_df(self.adata, group='0', key='t-test_overestim_var', pval_cutoff=0.01, log2fc_min=0)['names']
                    fig, venn_ax = plt.subplots()
                    venn3([set(wc),set(tt),set(tt_ov)], ('Wilcox','T-test','T-test_ov'), ax=venn_ax)
                    st.pyplot(venn_ax)
           


        except Exception as e:
            st.error(e)

    def visualize(self):
        try:
            with st.form(key="visualize_form"):

                st.subheader("Visualize plots")

                input_col1, input_col2, input_col3, _, _ = st.columns(5, gap="large")
                method = input_col1.radio(label="Method", options=['t-test', 't-test_overestim_var', 'wilcoxon', 'logreg'])
                n_genes = input_col2.text_input(label="n_genes", value=5, help="Number of genes to display in each cluster.")
                group_by = input_col2.selectbox(label="Group by", options=st.session_state.adata_state.current.adata.obs_keys())
                subcol1, _, _, _, _, _ = st.columns(6)
                viz_submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
                
                if viz_submit_btn:
                    plt.style.use('dark_background')
                    n_genes = int(n_genes)
                    heatmap, dotplot, stacked_violins, matrix_plot = st.tabs(['Heatmap', 'Dotplot', 'Stacked violins', 'Matrix plot'])
                    with st.spinner(text="Generating plots"):
                        with heatmap:
                            heatmap_ax = sc.pl.rank_genes_groups_heatmap(st.session_state.adata_state.current.adata, n_genes=n_genes, key=method, groupby=group_by, show_gene_labels=True)
                            st.pyplot(heatmap_ax)
                        with dotplot:
                            dotplot_ax = sc.pl.rank_genes_groups_dotplot(st.session_state.adata_state.current.adata, n_genes=n_genes, key=method, groupby=group_by)
                            st.pyplot(dotplot_ax)
                        with stacked_violins:
                            stacked_violins_ax = sc.pl.rank_genes_groups_stacked_violin(st.session_state.adata_state.current.adata, n_genes=n_genes, key=method, groupby=group_by)
                            st.pyplot(stacked_violins_ax)
                        with matrix_plot:
                            matrix_plot_ax = sc.pl.rank_genes_groups_matrixplot(st.session_state.adata_state.current.adata, n_genes=n_genes, key=method, groupby=group_by)
                            st.pyplot(matrix_plot_ax)
   

        except Exception as e:
            st.error(e)


try:
    sidebar = Sidebar()
    sidebar.show()

    st.title("Differential Gene Expression")

    sidebar.show_preview()

    dge = DGE(st.session_state.adata_state.current.adata.copy())

    sidebar.export_script()
    sidebar.delete_experiment_btn()
    sidebar.show_version()
    
except KeyError as ke:
    print("KeyError: ", ke)
    st.error("Couldn't find adata object in session, have you uploaded one?")
except Exception as e:
    print("Error: ", e)
    st.error(e)