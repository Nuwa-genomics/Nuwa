import streamlit as st
import scanpy as sc
import pickle
import pandas as pd
import warnings
import numpy as np

from models.AdataModel import AdataModel
from utils.plotting import highest_expr_genes_box_plot
from components.sidebar import *
from datetime import datetime

from database.database import SessionLocal
from sqlalchemy.orm import Session

from database.schemas import schemas
from utils.AdataState import AdataState
from time import sleep
import os
import re
import plotly.graph_objects as go


st.set_page_config(layout="wide", page_title='Nuwa', page_icon='🧬')

os.chdir('/app')

with open('css/common.css') as f:
    common_style = f"""
                <style>
                {f.read()}
                </style>
                """
    st.markdown(common_style, unsafe_allow_html=True)



st.set_option('deprecation.showPyplotGlobalUse', False)


class Preprocess:
    def __init__(self, adata):
        self.adata = adata
        self.conn: Session = SessionLocal()
        st.title("Preprocess")
    
        
    def save_adata(self):
        sc.write(filename=os.path.join(os.getenv('WORKDIR'), 'adata', st.session_state.adata_state.current.adata_name), adata=self.adata)
        st.session_state.adata_state.current.adata = self.adata
         

    def filter_highest_expr_genes(self):
        """
        Fraction of counts assigned to each gene over all cells. Computes, for each gene, the fraction of counts assigned to that gene within a cell. The n_top genes with the highest mean fraction over all cells are plotted as boxplots.
        
        Parameters
        ----------
        n_top_genes : int
            Number of top gene symbols to show.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/highest_expr_genes.png

        Example
        -------
        import os
        import scanpy as sc
        """
        with st.form(key="form_highest_expr"):
            st.subheader("Show highest expressed genes")
            n_top_genes = st.number_input(label="Number of genes", min_value=1, max_value=100, value=20, key="ni:pp:highly_variable:n_top_genes")
            subcol1, _, _ = st.columns(3)
            submit_btn = subcol1.form_submit_button(label="Filter", use_container_width=True)

            if submit_btn:
                with st.spinner(text="Calculating highest expressed genes"):
                    fig = highest_expr_genes_box_plot(self.adata, n_top=n_top_genes)
                    st.plotly_chart(fig)
                        
                        
    def remove_genes(self):
        with st.form(key="remove_genes_form"):
            st.subheader("Remove genes")
            remove_genes = st.multiselect(label="Genes", options=self.adata.var_names)
            subcol1, _, _ = st.columns(3)
            submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
            if submit_btn:
                with st.spinner(text="Removing genes"):
                    for gene in remove_genes:
                        remove_genes = self.adata.var_names.str.startswith(gene)
                        remove = np.array(remove_genes)
                        keep = np.invert(remove)
                        self.adata = self.adata[:,keep]
                        self.save_adata()
                    #write to script state
                    st.session_state["script_state"].add_script(f"#remove selected genes\nfor gene in {remove_genes}:\n\tremove_genes = adata.var_names.str.startswith(gene)\nremove = np.array(remove_genes)\nkeep = np.invert(remove)\nadata = self.adata[:,keep]")


    def filter_highly_variable_genes(self):
        try:
            with st.form(key="form_highly_variable"):
                st.subheader("Highly variable genes")
                min_mean = st.number_input(label="min mean", value=0.0125, key="input_highly_variable_min_mean")
                max_mean = st.number_input(label="max mean", value=3.0, key="input_highly_variable_max_mean")
                remove = st.toggle(label="Remove non-variable genes", value=False, help="By default, highly variable genes are only annoted. This option will remove genes without highly variable expression.")
                fn = 'figures/filter_genes_dispersion.pdf'
                subcol1, _, _ = st.columns(3)
                submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)

                if submit_btn:
                    with st.spinner(text="Calculating highly variable genes"):
                        sc.pp.normalize_total(self.adata, target_sum=1e4); 
                        sc.pp.log1p(self.adata)
                        sc.pp.highly_variable_genes(self.adata, min_mean=min_mean, max_mean=max_mean, min_disp=0.5)
                        #write to script state
                        st.session_state["script_state"].add_script("#Filter highly variable genes")
                        st.session_state["script_state"].add_script("sc.pp.normalize_total(adata, target_sum=1e4)")
                        st.session_state["script_state"].add_script("sc.pp.log1p(adata)")
                        st.session_state["script_state"].add_script(f"sc.pp.highly_variable_genes(adata, min_mean={min_mean}, max_mean={max_mean}, min_disp=0.5)")
                        
                        if remove:
                            self.adata = self.adata[:, self.adata.var.highly_variable]
                        
                        #make adata
                        self.save_adata()
                        ax = sc.pl.highly_variable_genes(self.adata, save=True)
                        st.pyplot(ax)
                        #add to script state
                        st.session_state["script_state"].add_script("sc.pl.highly_variable_genes(adata)")
        except Exception as e:
            st.toast(f"Failed to normalize data: {e}", icon="❌")

    def normalize_counts(self):
        st.subheader("Normalization")

        with st.form(key="form_normalize_total"):
            target_sum = st.number_input(label="Target sum", value=1, key="ni_target_sum")
            subcol_input1, subcol_input2 = st.columns(2, gap="medium")
            exclude_high_expr = subcol_input1.checkbox(label="Exclude highly_expr", value=False)
            log_transform_total = subcol_input2.checkbox(label="Log transform", value=False)

            subcol1, _, _ = st.columns(3)
            submit_btn = subcol1.form_submit_button(label="Apply", use_container_width=True)

            if submit_btn:
                sc.pp.normalize_total(self.adata, target_sum=target_sum, exclude_highly_expressed=exclude_high_expr)
                if log_transform_total:
                    sc.pp.log1p(self.adata)
                #write to script state
                st.session_state["script_state"].add_script(f"#Normalize counts (total)\nsc.pp.normalize_total(adata, target_sum={target_sum}, exclude_highly_expressed={exclude_high_expr})")
                #make adata
                self.save_adata()
                st.toast("Normalized data", icon='✅')

    def filter_cells(self):
        with st.form(key="form_filter_cells"):
            st.subheader("Filter Cells")
            subcol1, subcol2 = st.columns(2)
            with subcol1:
                min_count = st.number_input(label="Min count", min_value=1, value=None, key="filter_cell_min_count")
                min_genes = st.number_input(label="min genes for cell", min_value=1, value=None, key="filter_cell_min_genes")

            with subcol2:
                max_count = st.number_input(label="Max count", min_value=1, value=None, key="filter_cell_max_count")
                max_genes = st.number_input(label="max genes for cell", min_value=1, value=None, key="filter_cell_max_genes")
      
            subcol1, _, _ = st.columns(3)
            submit_btn = subcol1.form_submit_button(label="Apply", use_container_width=True)

            if submit_btn:
                sc.pp.filter_cells(self.adata, max_genes=max_genes, min_genes=min_genes, max_counts=max_count, min_counts=min_count)
                #write to script state
                st.session_state["script_state"].add_script(f"#Filter cells\nsc.pp.filter_cells(adata, max_genes={max_genes}, min_genes={min_genes}, max_counts={max_count}, min_counts={min_count})")
                #make adata
                self.save_adata()
                st.toast("Filtered cells", icon='✅')


    def filter_genes(self):
        with st.form(key="form_filter_genes"):
            st.subheader("Filter Genes")
            subcol1, subcol2 = st.columns(2)
            with subcol1:
                min_count = st.number_input(label="Min count", min_value=1, value=None, key="filter_gene_min_count")
                min_cells = st.number_input(label="min cells for gene", min_value=1, value=None, key="filter_gene_min_cells")

            with subcol2:
                max_count = st.number_input(label="Max count", min_value=1, value=None, key="filter_gene_max_count")
                max_cells = st.number_input(label="max cells for gene", min_value=1, value=None, key="filter_gene_max_cells")
            subcol1, _, _ = st.columns(3)
            submit_btn = subcol1.form_submit_button(label="Apply", use_container_width=True)
            if submit_btn:
                sc.pp.filter_genes(self.adata, max_cells=max_cells, min_cells=min_cells, max_counts=max_count, min_counts=min_count)
                #write to script state
                st.session_state["script_state"].add_script(f"#Filter genes\nsc.pp.filter_genes(adata, max_cells={max_cells}, min_cells={min_cells}, max_counts={max_count}, min_counts={min_count})")
                #make adata
                self.save_adata()
                st.toast("Filtered genes", icon='✅')


    def recipes(self):
        st.subheader("Preprocess Recipes")
        seurat_tab, weinreb17_tab, zheng17_tab = st.tabs(['Seurat', 'Weinreb17', 'Zheng17'])
        with seurat_tab:
            with st.form(key="form_seurat"):
                st.write("Parameters")
                log = st.checkbox(label="Log", value=False)
                subcol1, _, _ = st.columns(3)
                submit_btn = subcol1.form_submit_button(label='Apply', use_container_width=True)
                
                if submit_btn:
                    sc.pp.recipe_seurat(self.adata, log=log) 
                    st.session_state["script_state"].add_script("#Apply preprocess recipe\nsc.pp.recipe_seurat(adata)")
                    self.save_adata()
                    st.toast(f"Applied recipe: Seurat", icon='✅')
        
        with weinreb17_tab:
            with st.form(key="form_weinreb17"):
                st.write("Parameters")
                col1, col2, col3 = st.columns(3)
                mean_threshold = col1.number_input(label="Mean threshold", value=0.01, step=0.01)
                cv_threshold = col2.number_input(label="CV threshold", value=2.0, step=1.0)
                n_pcs = col3.number_input(label="n_pcs", min_value=1, value=50, step=1, format="%i")
                log = st.checkbox(label="Log", value=False)
                subcol1, _, _ = st.columns(3)
                submit_btn = subcol1.form_submit_button(label='Apply', use_container_width=True)
                if submit_btn:
                    sc.pp.recipe_weinreb17(self.adata, log=log, mean_threshold=mean_threshold, cv_threshold=cv_threshold, n_pcs=n_pcs)
                    st.session_state["script_state"].add_script("#Apply preprocess recipe\nsc.pp.recipe_weinreb17(adata)")
                    self.save_adata()
                    st.toast(f"Applied recipe: Weinreb17", icon='✅')

        with zheng17_tab:
            with st.form(key="form_zheng17"):
                st.write("Parameters")
                n_top_genes = st.number_input(label="n_top_genes", key="ni_zheng17_n_genes", min_value=1, max_value=self.adata.n_vars, value=1000 if self.adata.n_vars >= 1000 else self.adata.n_vars)
                log = st.checkbox(label="Log", value=False)
                subcol1, _, _ = st.columns(3)
                submit_btn = subcol1.form_submit_button(label='Apply', use_container_width=True)
                if submit_btn:
                    sc.pp.recipe_zheng17(self.adata, log=log, n_top_genes=n_top_genes)
                    st.session_state["script_state"].add_script("#Apply preprocess recipe\nsc.pp.recipe_zheng17(adata)")
                    self.save_adata()
                    st.toast(f"Applied recipe: Zheng17", icon='✅')

    
    def annotate_mito(self):
        with st.form(key="form_annotate_mito"):
            st.subheader("Annotate Mitochondrial Genes", help="Filter mitochrondrial gene counts. All mitochrondrial genes \
                        are by default annotated and placed in the 'mt' variable.")
            
            self.adata.var['mt'] = self.adata.var_names.str.startswith(('MT-', 'mt-'))
            sc.pp.calculate_qc_metrics(self.adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
                
            st.text(f"Found {self.adata.var.mt.sum()} mitochondrial genes")
            
            mito_container = st.empty()

            def plot_charts(color=None):
                if color == 'None':
                    color=None
                mito_container.empty()
                subcol1, subcol2 = mito_container.columns(2, gap="small")
                with subcol1:
                    ax_scatter = sc.pl.scatter(self.adata, x='total_counts', y='pct_counts_mt', color=color)
                    with st.expander(label="Scatter"):
                        st.pyplot(ax_scatter)

                with subcol2:
                    ax_violin = sc.pl.violin(self.adata, 'pct_counts_mt', groupby=color)
                    with st.expander(label="Violin"):
                        st.pyplot(ax_violin)

            
            subcol_input1, subcol_input2 = st.columns(2)
            options = np.append('None', st.session_state.adata_state.current.adata.obs_keys())
            color_key_mito = subcol_input1.selectbox(label="Color key", options=options)
            ni_pct_counts_mt = subcol_input2.number_input(label="max pct_counts_mt", key="ni_pct_counts_mt", min_value=0, value=100)

            plot_charts()

            subcol1, _, _ = st.columns(3)
            mito_annot_submit_btn = subcol1.form_submit_button(label="Apply", use_container_width=True)

            if mito_annot_submit_btn:
                plot_charts(color_key_mito)
                self.adata = self.adata[self.adata.obs.pct_counts_mt < ni_pct_counts_mt, :]
                #write to script state
                st.session_state["script_state"].add_script("#Filter mitochondrial genes")
                st.session_state["script_state"].add_script("sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')")
                st.session_state["script_state"].add_script("sc.pl.violin(adata, 'pct_counts_mt')")
                st.session_state["script_state"].add_script("sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)")
                st.session_state["script_state"].add_script(f"adata = adata[adata.obs.pct_counts_mt < {ni_pct_counts_mt}, :]")
                #make adata
                self.save_adata()
                st.toast("Filtered mitochondrial genes", icon="✅")

            

    def annotate_ribo(self):
        with st.form(key="form_annotate_ribo"):
            st.subheader("Annotate Ribosomal Genes", help="Filter ribosomal gene counts. All ribosomal genes \
                        are by default annotated and placed in the 'ribo' variable.")

            with st.spinner(text="Fetching ribosomal genes"):
                ribo_url = "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"
                ribo_genes = pd.read_table(ribo_url, skiprows=2, header=None)
                self.adata.var['ribo'] = self.adata.var_names.isin(ribo_genes[0].values)
                sc.pp.calculate_qc_metrics(self.adata, qc_vars=['ribo'], percent_top=None, log1p=False, inplace=True)
                
                st.text(f"Found {self.adata.var['ribo'].sum()} ribosomal genes")

                ribo_container = st.empty()

            def plot_charts(color=None):
                if color == 'None':
                    color=None
                ribo_container.empty()
                subcol1, subcol2 = ribo_container.columns(2, gap="small")
                with subcol1:
                    ax_scatter = sc.pl.scatter(self.adata, x='total_counts', y='pct_counts_ribo', color=color)
                    with st.expander(label="Scatter"):
                        st.pyplot(ax_scatter)

                with subcol2:
                    ax_violin = sc.pl.violin(self.adata, 'pct_counts_ribo', groupby=color)
                    with st.expander(label="Violin"):
                        st.pyplot(ax_violin)

            
            subcol_input1, subcol_input2 = st.columns(2)
            options = np.append('None', st.session_state.adata_state.current.adata.obs_keys())
            color_key_mito = subcol_input1.selectbox(label="Color key", options=options)
            ni_pct_counts_ribo = subcol_input2.number_input(label="max pct_counts_ribo", key="ni_pct_counts_ribo", min_value=0, value=100)

            plot_charts()

            subcol1, _, _ = st.columns(3)
            mito_annot_submit_btn = subcol1.form_submit_button(label="Apply", use_container_width=True)

            if mito_annot_submit_btn:
                plot_charts(color_key_mito)
                self.adata = self.adata[self.adata.obs.pct_counts_ribo < ni_pct_counts_ribo, :]
                #add to script adata
                st.session_state["script_state"].add_script("#Filter ribosomal genes")
                st.session_state["script_state"].add_script("sc.pl.scatter(adata, x='total_counts', y='pct_counts_ribo')")
                st.session_state["script_state"].add_script("sc.pl.violin(adata, 'pct_counts_ribo')")
                st.session_state["script_state"].add_script("sc.pp.calculate_qc_metrics(adata, qc_vars=['ribo'], percent_top=None, log1p=False, inplace=True)")
                st.session_state["script_state"].add_script(f"adata = adata[adata.obs.pct_counts_ribo < {ni_pct_counts_ribo}, :]")
                #make adata
                self.save_adata()
                st.toast("Filtered ribosomal genes", icon="✅")
                
                
    def annotate_hb(self):
        with st.form(key="form_annotate_hb"):
            st.subheader("Annotate haemoglobin genes", help="Filter haemoglobin gene counts. All haemoglobin genes \
                        are by default annotated and placed in the 'hb' variable.")
            
            # hemoglobin genes.
            self.adata.var['hb'] = self.adata.var_names.str.contains(("^HB[^(P)]"))
            
            sc.pp.calculate_qc_metrics(self.adata, qc_vars=['hb'], percent_top=None, log1p=False, inplace=True)
            
            st.text(f"Found {self.adata.var.hb.sum()} haemoglobin genes")

            hb_container = st.empty()

            def plot_charts(color=None):
                if color == 'None':
                    color=None
                hb_container.empty()
                subcol1, subcol2 = hb_container.columns(2, gap="small")
                with subcol1:
                    ax_scatter = sc.pl.scatter(self.adata, x='total_counts', y='pct_counts_hb', color=color)
                    with st.expander(label="Scatter"):
                        st.pyplot(ax_scatter)

                with subcol2:
                    ax_violin = sc.pl.violin(self.adata, 'pct_counts_hb', groupby=color)
                    with st.expander(label="Violin"):
                        st.pyplot(ax_violin)

            
            subcol_input1, subcol_input2 = st.columns(2)
            options = np.append('None', st.session_state.adata_state.current.adata.obs_keys())
            color_key_mito = subcol_input1.selectbox(label="Color key", options=options)
            ni_pct_counts_hb = subcol_input2.number_input(label="max pct_counts_hb", key="ni_pct_counts_hb", min_value=0, value=100)

            plot_charts()

            subcol1, _, _ = st.columns(3)
            hb_annot_submit_btn = subcol1.form_submit_button(label="Apply", use_container_width=True)

            if hb_annot_submit_btn:
                plot_charts(color_key_mito)
                self.adata = self.adata[self.adata.obs.pct_counts_ribo < ni_pct_counts_hb, :]
                #add to script adata
                st.session_state["script_state"].add_script("#Filter haemoglobin genes")
                st.session_state["script_state"].add_script("sc.pl.scatter(adata, x='total_counts', y='pct_counts_hb')")
                st.session_state["script_state"].add_script("sc.pl.violin(adata, 'pct_counts_hb')")
                st.session_state["script_state"].add_script("sc.pp.calculate_qc_metrics(adata, qc_vars=['hb'], percent_top=None, log1p=False, inplace=True)")
                st.session_state["script_state"].add_script(f"adata = adata[adata.obs.pct_counts_hb < {ni_pct_counts_hb}, :]")
                #make adata
                self.save_adata()
                st.toast("Filtered ribosomal genes", icon="✅")
                
                

    def run_scrublet(self):
        with st.form(key="scrublet_form"):
            st.subheader("Doublet Prediction", help="Use Scrublet to remove cells predicted to be doublets.")
            col1, col2, col3 = st.columns(3)
            sim_doublet_ratio = col1.number_input(label="Sim doublet ratio", value=2.00, key="ni_sim_doublet_ratio")
            expected_doublet_rate = col2.number_input(label="Exp doublet rate", value=0.05, key="ni_expected_doublet_rate")
            stdev_doublet_rate = col3.number_input(label="stdev_doublet_rate", value=0.02, key="ni_stdev_doublet_rate")
            batch_key = st.selectbox(label="Batch key", key="sb_scrublet_batch_key", options=np.append('None', self.adata.obs_keys()))
            subcol1, _, _ = st.columns(3)
            scrublet_submit = subcol1.form_submit_button(label="Filter", use_container_width=True)

            if scrublet_submit:
                with st.spinner("Running scrublet"):
                    if batch_key == 'None':
                        batch_key = None
                    sc.external.pp.scrublet(self.adata, sim_doublet_ratio=sim_doublet_ratio, 
                        expected_doublet_rate=expected_doublet_rate, stdev_doublet_rate=stdev_doublet_rate, batch_key=batch_key, random_state=42)
                    #write to script state
                    st.session_state["script_state"].add_script(f"#Detect doublets with scrublet\nadata_scrublet = sc.external.pp.scrublet(adata, sim_doublet_ratio={sim_doublet_ratio}, expected_doublet_rate={expected_doublet_rate})")
                    #make adata
                    self.save_adata()
                    st.toast("Completed doublet predictions", icon="✅")
                    
                    
    def regress_out(self):
        with st.form(key="regress_out_form"):
            st.subheader("Regress out", help="Uses linear regression to remove unwanted sources of variation.")
            regress_keys = st.multiselect(label="Keys", options=self.adata.obs_keys(), key="ms_regress_out_keys")
            subcol1, _, _ = st.columns(3)
            submit_btn = subcol1.form_submit_button(label="Apply", use_container_width=True)
            if submit_btn:
                if st.session_state.ms_regress_out_keys:
                    sc.pp.regress_out(self.adata, keys=regress_keys)
                    #write to script state
                    st.session_state["script_state"].add_script(f"#Regress out\nsc.pp.regress_out(adata, keys={st.session_state.ms_regress_out_keys})")
                    st.toast("Successfully regressed out data", icon="✅")
                    self.save_adata()
                else:
                    st.toast("No option selected, not regressing data.", icon="ℹ️")
            
            
    def scale_to_unit_variance(self):
        with st.form(key="scale_to_unit_variance_form"):
            st.subheader("Scale to unit variance")
            st.number_input(label="Max value", value=10, key="ni_scale_data_max_value")
            subcol1, _, _ = st.columns(3)
            btn_scale_data_btn = subcol1.form_submit_button(label="Apply", use_container_width=True)
            if btn_scale_data_btn:
                if st.session_state.ni_scale_data_max_value:
                    sc.pp.scale(self.adata, max_value=st.session_state.ni_scale_data_max_value)
                    #write to script state
                    st.session_state["script_state"].add_script(f"#Scale to unit variance\nsc.pp.scale(adata, max_value={st.session_state.ni_scale_data_max_value})")
                    self.save_adata()
                    st.toast("Successfully scaled data", icon="✅")
                else:
                    st.toast("Max value cannot be blank", icon="❌")
            
    

    def downsample_data(self):
        st.subheader("Downsample data")
        counts_per_cell, total_counts = st.tabs(["counts_per_cell", "total_counts"])
        with counts_per_cell:
            with st.form(key="downsample_form_counts_per_cell"):
                counts_per_cell = st.number_input(label="Counts per cell", key="ni_downsample_counts_per_cell", help="Target total counts per cell. If a cell has more than 'counts_per_cell', it will be downsampled to this number. Resulting counts can be specified on a per cell basis by passing an array.Should be an integer or integer ndarray with same length as number of obs.")
                subcol1, _, _ = st.columns(3)
                btn_downsample_counts_per_cell = subcol1.form_submit_button(label="Apply", use_container_width=True)
                if btn_downsample_counts_per_cell:
                    sc.pp.downsample_counts(self.adata, counts_per_cell=counts_per_cell)
                    #write to script state
                    st.session_state["script_state"].add_script(f"#Downsample dataset\nsc.pp.downsample_counts(adata, counts_per_cell={counts_per_cell})")
                    self.save_adata()
                    st.toast("Successfully downsampled data per cell", icon="✅")
        with total_counts:
            with st.form(key="downsample_form_total_counts"):
                total_counts = st.number_input(label="Total counts", key="ni_downsample_total_counts", help="Target total counts. If the count matrix has more than total_counts it will be downsampled to have this number.")
                subcol1, _, _ = st.columns(3)
                btn_downsample_total_counts = subcol1.form_submit_button(label="Apply", use_container_width=True)
                if btn_downsample_total_counts:
                    sc.pp.downsample_counts(self.adata, total_counts=total_counts)
                    #write to script state
                    st.session_state["script_state"].add_script(f"#Downsample dataset\nsc.pp.downsample_counts(adata, total_counts={total_counts})")
                    self.save_adata()
                    st.toast("Successfully downsampled data by total counts", icon="✅")

            
    def subsample_data(self):
        st.subheader("Subsample data")
        n_obs, fraction = st.tabs(["n_obs", "fraction"])
        with n_obs:
            with st.form(key="subsample_form_n_obs"):
                n_obs = st.number_input(label="n obs", key="ni_subsample_n_obs", help="Subsample to this number of observations.")
                subcol1, _, _ = st.columns(3)
                btn_subsample_n_obs = subcol1.form_submit_button(label="Apply", use_container_width=True)
                if btn_subsample_n_obs:
                    sc.pp.subsample(self.adata, n_obs=st.session_state.ni_n_obs)
                    #write to script session
                    st.session_state["script_state"].add_script(f"#Subsample dataset\nsc.pp.subsample(adata, n_obs={st.session_state.ni_n_obs}, fraction={st.session_state.ni_subsample_fraction})")
                    self.save_adata()
                    st.toast("Successfully subsampled data", icon="✅")
        with fraction:
            with st.form(key="subsample_form_fraction"):
                fraction = st.number_input(label="subsample_fraction", key="ni_subsample_fraction", help="Subsample this fraction of the number of observations.")
                subcol1, _, _ = st.columns(3)
                btn_subsample_fraction = subcol1.form_submit_button(label="Apply", use_container_width=True)
                if btn_subsample_fraction:
                    sc.pp.subsample(self.adata, fraction=st.session_state.ni_subsample_fraction)
                    #write to script session
                    st.session_state["script_state"].add_script(f"#Subsample dataset\nsc.pp.subsample(adata, n_obs={st.session_state.ni_n_obs}, fraction={st.session_state.ni_subsample_fraction})")
                    self.save_adata()
                    st.toast("Successfully subsampled data", icon="✅")
                    
    def batch_effect_removal(self):
        with st.form(key="batch_effect_removal_form"):
            st.subheader("Batch effect correction", help="Uses Combat to correct non-biological differences caused by batch effect.")
            index = 0
            for i, obs in enumerate(self.adata.obs_keys()):
                if obs.lower().replace("_", "").__contains__("batch"):
                    index = i
            key = st.selectbox(label="Key", options=self.adata.obs_keys(), key="sb_batch_effect_key", index=index)
            covariates = st.multiselect(placeholder="Optional", label="Covariates", options=self.adata.obs_keys())
            subcol1, _, _ = st.columns(3)
            btn_batch_effect_removal = subcol1.form_submit_button(label="Apply", use_container_width=True)
            if btn_batch_effect_removal:
                with st.spinner(text="Running Combat batch effect correction"):
                    sc.pp.combat(self.adata, key=key, covariates=covariates)
                    #write to script state
                    st.session_state["script_state"].add_script(f"#Correct batch effect\nsc.pp.combat(adata, key={key}, covariates={covariates})")
                self.save_adata()
                st.toast("Batch corrected data", icon='✅')
                
                
    def pca(self):
        with st.form(key="pca_pp_form"):
            st.subheader("PCA")
            
            def run_pca(adata):
                with st.spinner(text="Running PCA"):
                    sc.pp.pca(adata, random_state=42)
                    st.session_state['pp_df_pca'] = pd.DataFrame({'pca1': adata.obsm['X_pca'][:,0], 'pca2': adata.obsm['X_pca'][:,1], 'color': adata.obs[f'{st.session_state.sb_pca_color_pp}']})  
                    pca_empty.empty()
                    pca_empty.scatter_chart(data=st.session_state['pp_df_pca'], x='pca1', y='pca2', color='color', size=18)
                
               
            index = 0      
            for i, item in enumerate(self.adata.obs_keys()):
                  if item.lower().replace("_", "").__contains__("batch"): #give precedence to batch if present since it is relevant to preprocessing
                      index = i           
            pca_color = st.selectbox(label="Color", options=self.adata.obs_keys(), key="sb_pca_color_pp", index=index)
            subcol1, _, _ = st.columns(3)
            pca_pp_btn = subcol1.form_submit_button("Apply", use_container_width=True)
            pca_empty = st.empty()
            
            run_pca(self.adata)
            
            if pca_pp_btn:
                st.session_state["script_state"].add_script("#Run PCA")
                st.session_state["script_state"].add_script("sc.pp.pca(adata, random_state=42)")
                st.session_state["script_state"].add_script("sc.pl.pca(adata, random_state=42)")
                run_pca(self.adata)
                
    
    def predict_sex(self):
        
        st.subheader("Measure gene counts")
        single_dataset, subsample = st.tabs(['Current dataset', 'Subsample'])
        
        with single_dataset:
            with st.form(key="sex_predict_single_dataset"):
                st.subheader("Collective counts across dataset")
                options = np.append(['total_counts', 'n_genes_by_counts'], self.adata.var_names)
                gene = st.selectbox(label="Gene (e.g. XIST for detecting sex)", options=options)
                subcol_btn1, _, _ = st.columns(3)
                submit_btn = subcol_btn1.form_submit_button(label="Run", use_container_width=True)
                if submit_btn:
                    with st.spinner(text="Locating genes"):
                        if gene == "total_counts" or gene == "n_genes_by_counts":
                            arr = np.array([f'{st.session_state.adata_state.current.adata_name}'])
                            df = pd.DataFrame({f"{gene}": self.adata.obs[gene], "Dataset": np.repeat(arr, self.adata.n_obs)})
                            bar_chart, violin_plot = st.tabs(['Bar chart', 'Violin plot'])
                            bar_chart.bar_chart(data=df, x="Dataset", y=f"{gene}")
                            with violin_plot:
                                violin_plot_single_dataset = sc.pl.violin(self.adata, [f'{gene}'], jitter=0.4, rotation = 45)
                                st.pyplot(violin_plot_single_dataset)
                        else:
                            self.adata.obs["gene-counts"] = self.adata.X[:,self.adata.var_names.str.match(f'{gene}')].toarray()
                            arr = np.array([f'{st.session_state.adata_state.current.adata_name}'])
                            df = pd.DataFrame({f'{gene} count': self.adata.obs["gene-counts"], "Dataset": np.repeat(arr, self.adata.n_obs)})
                            st.bar_chart(data=df, x="Dataset", y=f'{gene} count')
                            #write to script state
                            st.session_state["script_state"].add_script("#Measure gene counts in single dataset")
                            st.session_state["script_state"].add_script(f"adata.obs['gene-counts'] = adata.X[:,adata.var_names.str.match('{gene}')].toarray()")
                            st.session_state["script_state"].add_script(f"arr = np.array(['{st.session_state.adata_state.current.adata_name}'])")
                            st.session_state["script_state"].add_script(f"df = pd.DataFrame('{gene} count': adata.obs['gene-counts'], 'Dataset': np.repeat(arr, adata.n_obs))")
        with subsample:
            with st.form(key="sex_predict_multiple_datasets"):
                st.subheader("Subsample counts in dataset")
                gene_options = np.append(['total_counts', 'n_genes_by_counts'], self.adata.var_names)
                batch_key_sex_pred = st.selectbox(label="Obs key", options=self.adata.obs_keys(), key="sb_sex_pred_batch_key")
                gene = st.selectbox(label="Gene (e.g. XIST for detecting sex)", options=gene_options)
                subcol_btn1, _, _ = st.columns(3)
                submit_btn = subcol_btn1.form_submit_button(label="Run", use_container_width=True)
                if submit_btn:
                    with st.spinner(text="Locating genes"):
                        if gene == "total_counts" or gene == "n_genes_by_counts":
                            df = pd.DataFrame({f'{gene}': self.adata.obs[f"{gene}"], f"{batch_key_sex_pred}": self.adata.obs[f"{batch_key_sex_pred}"]})
                            bar_chart, violin_plot = st.tabs(['Bar chart', 'violin plot'])
                            bar_chart.bar_chart(data=df, x=f"{batch_key_sex_pred}", y=f'{gene}')
                            with violin_plot:
                                violin_plot_subsample_dataset = sc.pl.violin(self.adata, [f'{gene}'], jitter=0.4, groupby = f'{batch_key_sex_pred}', rotation = 45)
                                st.pyplot(violin_plot_subsample_dataset)
                        else:
                            self.adata.obs["gene-counts"] = self.adata.X[:,self.adata.var_names.str.match(f'{gene}')].toarray()
                            df = pd.DataFrame({f'{gene} count': self.adata.obs["gene-counts"], f"{batch_key_sex_pred}": self.adata.obs[f"{batch_key_sex_pred}"]})
                            st.bar_chart(data=df, x=f"{batch_key_sex_pred}", y=f'{gene} count')
                            #write to script state
                            st.session_state["script_state"].add_script("#Measure gene counts across datasets")
                            st.session_state["script_state"].add_script(f"adata.obs['gene-counts'] = adata.X[:,adata.var_names.str.match('{gene}')].toarray()")
                            st.session_state["script_state"].add_script(f"df = pd.DataFrame('{gene} count': adata.obs['gene-counts'], '{batch_key_sex_pred}': adata.obs['{batch_key_sex_pred}'])")
                  
    def cell_cycle_scoring(self):
        st.subheader("Cell cycle score")
        col1, col2, col3 = st.columns(3)
        col1.file_uploader(label="Cell cycle genes", type=["csv", "tsv"], accept_multiple_files=False, key="cell_cycle_file_uploader", help="File must be a csv with the gene name and corresponding phase. Fore example: \ngenes,phase\n0,MCM5,s_genes\n1,PCNA,s_genes\n2,TYMS,s_genes\nNote: csv files may also have an index and header.")
     
        with st.form(key="cell_cycle_scoring_form"):

            form_col1, form_col2, form_col3 = st.columns(3, gap="large")
            
            if not st.session_state.cell_cycle_file_uploader:
                info_col1, info_col2, info_col3 = st.columns(3)
                info_col1.info("No file uploaded. CSV must include header values.")
            elif st.session_state.cell_cycle_file_uploader: 
                delim = "\t" if st.session_state.cell_cycle_file_uploader.type == "tsv" else ","
                file = st.session_state.cell_cycle_file_uploader
                st.session_state["pp_cell_cycle_marker_genes_df"] = pd.read_csv(file, delimiter=delim)
                df = st.session_state["pp_cell_cycle_marker_genes_df"]

                form_col1.write("Preview")
                form_col1.dataframe(df.head(4), use_container_width=True)

                form_col2.write("Columns")
                gene_column = form_col2.selectbox(label="Gene column", options=df.columns, help="Specify which column contains the gene names.", key="sb_gene_col_cell_cycle")
                phase_column = form_col2.selectbox(label="Phase column", options=df.columns, help="Specify which column contains the phase (e.g. s phase)", key="sb_phase_col_cell_cycle")
                
                form_col3.write("Plot")
                group_by = form_col3.selectbox(label="Group by", options=np.append('None', self.adata.obs_keys()), key="sb_group_cell_cycle")
                plot_col1, plot_col2 = form_col3.columns(2)
                bandwidth = plot_col1.number_input(label="Bandwidth", min_value=0.1, max_value=1.0, value=0.4, step=0.1)
                jitter = plot_col2.number_input(label="Jitter", min_value=0.1, max_value=1.0, value=0.4, step=0.1)
            
            subcol1, _, _, _, _, _, _, _, _ = st.columns(9)
            submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True, disabled=(not "pp_cell_cycle_marker_genes_df" in st.session_state))
            cell_cycle_container = st.empty()

            if submit_btn:
                with st.spinner(text="Running cell cycle scoring"):
                    delim = "\t" if st.session_state.cell_cycle_file_uploader.type == "tsv" else ","
                    df = st.session_state["pp_cell_cycle_marker_genes_df"]

                    gene_column_index = df.columns.get_loc(gene_column)
                    phase_column_index = df.columns.get_loc(phase_column)

                    s_genes = df.iloc[:, phase_column_index].str.contains("s", flags=re.IGNORECASE, regex=True)
                    g2m_genes = df.iloc[:, phase_column_index].str.contains("g2m", flags=re.IGNORECASE, regex=True)

                    s_genes = df[s_genes].iloc[:, gene_column_index].values
                    g2m_genes = df[g2m_genes].iloc[:, gene_column_index].values

                    #cell_cycle_genes = [x for x in cell_cycle_genes if x in self.adata.var_names]
                    self.adata.raw = self.adata
                    sc.pp.normalize_per_cell(self.adata, counts_per_cell_after=1e4)
                    sc.pp.log1p(self.adata)
                    sc.pp.scale(self.adata)
                    if group_by == 'None':
                        group_by = None
                    sc.tl.score_genes_cell_cycle(self.adata, s_genes=s_genes, g2m_genes=g2m_genes)

                    #using matplotlib
                    #cell_cycle_ax = sc.pl.violin(self.adata, ['S_score', 'G2M_score'], jitter=0.4, groupby = group_by, rotation=45)
                    #cell_cycle_container.pyplot(cell_cycle_ax)

                    fig = go.Figure()

                    s_score_df = pd.DataFrame({'phase': np.repeat('S_score', len(self.adata.obs['S_score'])), 'score': self.adata.obs['S_score']})
                    g2m_score_df = pd.DataFrame({'phase': np.repeat('G2M_score', len(self.adata.obs['G2M_score'])), 'score': self.adata.obs['G2M_score']})

                    violin_df = pd.concat([s_score_df, g2m_score_df])

                    if group_by != None:

                        violin_df["group"] = self.adata.obs[group_by]

                        fig.add_trace(go.Violin(x=violin_df['group'][violin_df['phase'] == 'S_score'], 
                            y=violin_df['score'][violin_df['phase'] == 'S_score'],
                            legendgroup='S', scalegroup='S', name='S',
                            bandwidth=bandwidth, jitter=jitter, line_color='blue')
                        )

                        fig.add_trace(go.Violin(x=violin_df['group'][violin_df['phase'] == 'G2M_score'], 
                            y=violin_df['score'][violin_df['phase'] == 'G2M_score'],
                            legendgroup='G2M', scalegroup='G2M', name='G2M',
                            bandwidth=bandwidth, jitter=jitter, line_color='orange')
                        )

                        fig.update_traces(meanline_visible=True)
                        fig.update_layout(violingap=0, violinmode='group', xaxis_title=group_by, yaxis_title="Score", legend_title="Phase") #add legend title

                    else:

                        fig.add_trace(go.Violin(x=s_score_df['phase'], 
                            y=s_score_df['score'],
                            legendgroup='S', scalegroup='S', name='S',
                            bandwidth=bandwidth, jitter=jitter, line_color='blue')
                        )

                        fig.add_trace(go.Violin(x=g2m_score_df['phase'], 
                            y=g2m_score_df['score'],
                            legendgroup='G2M', scalegroup='G2M', name='G2M',
                            bandwidth=bandwidth, jitter=jitter, line_color='orange')
                        )

                        fig.update_traces(meanline_visible=True)
                        fig.update_layout(violingap=0, violinmode='overlay', xaxis_title="Phase", yaxis_title="Score", legend_title="Phase") #add legend title

                    
                    st.markdown("""<div style='margin-left: 20px; display: flex; align-items: center; justify-content: center;'><h1 style='text-align: center; font-size: 2rem;'>Cell cycle score</h1></div>""", unsafe_allow_html=True)

                    st.plotly_chart(fig, use_container_width=True)


try:
    sidebar = Sidebar()
    
    sidebar.show()
    
    preprocess = Preprocess(st.session_state.adata_state.current.adata.copy())

    col1, col2, col3 = st.columns(3, gap="medium")

    with col1:
        preprocess.filter_highest_expr_genes()
        preprocess.remove_genes()
        preprocess.filter_highly_variable_genes()
        preprocess.normalize_counts()
        preprocess.regress_out()
        preprocess.subsample_data()
        preprocess.downsample_data()
        

    with col2:
        preprocess.filter_cells()
        preprocess.filter_genes()
        preprocess.run_scrublet()
        preprocess.recipes()
        preprocess.batch_effect_removal()
        preprocess.predict_sex()
        
            
    with col3:
        preprocess.annotate_mito()
        preprocess.annotate_ribo()
        preprocess.annotate_hb()
        preprocess.pca()
        preprocess.scale_to_unit_variance()
        
        
    preprocess.cell_cycle_scoring()
        

    sidebar.show_preview()
    sidebar.export_script()
    sidebar.delete_experiment_btn()
    sidebar.show_version()

except KeyError as ke:
    print("KeyError: ", ke)
    st.error("Couldn't find adata object in session, have you uploaded one?")
    st.error(ke)
except Exception as e:
    print('Error: ', e)
    st.error(e)
