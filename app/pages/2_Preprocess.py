import streamlit as st
import scanpy as sc
import pickle
import pandas as pd
import warnings

from models.AdataModel import AdataModel
from components.sidebar import *
from datetime import datetime

from database.database import SessionLocal
from sqlalchemy.orm import Session

from database.schemas import schemas
from utils.AdataState import AdataState
from time import sleep
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



st.set_option('deprecation.showPyplotGlobalUse', False)


class Preprocess:
    def __init__(self, adata):
        self.adata = adata
        self.conn: Session = SessionLocal()
        st.title("Preprocess")
         

    def filter_highest_expr_genes(self):
        with st.form(key="form_highest_expr"):
            st.subheader("Show highest expressed genes")
            num_genes = st.number_input(label="Number of genes", min_value=1, max_value=100, value=20, key="n_top_genes")
            fn = 'figures/highest_expr_genes.pdf'
            subcol1, _, _ = st.columns(3)
            submit_btn = subcol1.form_submit_button(label="Filter", use_container_width=True)

            if submit_btn:
                with st.spinner(text="Calculating highest expressed genes"):
                    with st.expander(label="Show figure"):
                        ax = sc.pl.highest_expr_genes(self.adata, n_top=num_genes)
                        st.pyplot(ax)


    def filter_highly_variable_genes(self):
        try:
            with st.form(key="form_highly_variable"):
                st.subheader("Show highly variable genes")
                min_mean = st.number_input(label="min mean", value=0.0125, key="input_highly_variable_min_mean")
                max_mean = st.number_input(label="max mean", value=3.0, key="input_highly_variable_max_mean")
                fn = 'figures/filter_genes_dispersion.pdf'
                subcol1, _, _ = st.columns(3)
                submit_btn = subcol1.form_submit_button(label="Filter", use_container_width=True)

                if submit_btn:
                    with st.spinner(text="Calculating highly variable genes"):
                        with st.expander(label="Show figure"):
                            sc.pp.normalize_total(self.adata, target_sum=1e4)
                            sc.pp.log1p(self.adata)
                            sc.pp.highly_variable_genes(self.adata, min_mean=min_mean, max_mean=max_mean, min_disp=0.5)
                            #make adata
                            st.session_state.adata_state.add_adata(AdataModel(
                                work_id=st.session_state.current_workspace.id, adata=self.adata, 
                                filename=os.path.join(os.getenv('WORKDIR'), "adata", "adata_pp.h5ad"), adata_name="adata_pp"))
                            st.session_state.adata_state.switch_adata(adata_name="adata_pp")
                            ax = sc.pl.highly_variable_genes(self.adata)
                            st.pyplot(ax)
        except Exception as e:
            st.toast(f"Failed to normalize data: {e}", icon="❌")

    def normalize_counts(self):
        st.subheader("Normalization")
        tab_total, tab_per_cell = st.tabs(['Total', 'Per cell'])

        with tab_total:
            with st.form(key="form_normalize_total"):
                target_sum = st.number_input(label="Target sum", value=1, key="ni_target_sum")
                exclude_high_expr = st.checkbox(label="Exclude highly expressed", value=False)

                subcol1, _, _ = st.columns(3)
                submit_btn = subcol1.form_submit_button(label="Apply", use_container_width=True)

                if submit_btn:
                    sc.pp.normalize_total(self.adata, target_sum=target_sum, exclude_highly_expressed=exclude_high_expr)
                    #make adata
                    st.session_state.adata_state.add_adata(AdataModel(
                        work_id=st.session_state.current_workspace.id, adata=self.adata, 
                        filename=os.path.join(os.getenv('WORKDIR'), "adata", "adata_pp.h5ad"), adata_name="adata_pp"))
                    st.session_state.adata_state.switch_adata(adata_name="adata_pp")
                    st.toast("Normalized data", icon='✅')
        with tab_per_cell:
            with st.form(key="form_normalize_per_cell"):
                counts_per_cell_after = st.number_input(label="Counts per cell after")

                subcol1, _, _ = st.columns(3)
                submit_btn = subcol1.form_submit_button(label="Apply", use_container_width=True)

                if submit_btn:
                    sc.pp.normalize_per_cell(self.adata, counts_per_cell_after=counts_per_cell_after)
                    #make adata
                    st.session_state.adata_state.add_adata(AdataModel(
                        work_id=st.session_state.current_workspace.id, adata=self.adata, 
                        filename=os.path.join(os.getenv('WORKDIR'), "adata", "adata_pp.h5ad"), adata_name="adata_pp"))
                    st.session_state.adata_state.switch_adata(adata_name="adata_pp")
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
                #make adata
                st.session_state.adata_state.add_adata(AdataModel(
                    work_id=st.session_state.current_workspace.id, adata=self.adata, 
                    filename=os.path.join(os.getenv('WORKDIR'), "adata", "adata_pp.h5ad"), adata_name="adata_pp"))
                st.session_state.adata_state.switch_adata(adata_name="adata_pp")
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
                #make adata
                st.session_state.adata_state.add_adata(AdataModel(
                    work_id=st.session_state.current_workspace.id, adata=self.adata, 
                    filename=os.path.join(os.getenv('WORKDIR'), "adata", "adata_pp.h5ad"), adata_name="adata_pp"))
                st.session_state.adata_state.switch_adata(adata_name="adata_pp")
                st.toast("Filtered genes", icon='✅')


    def recipes(self):
        with st.form(key="form_recipes"):
            st.subheader("Preprocess Recipes")
            recipe = st.selectbox(label="Recipe", key="sb_pp_recipe", options=(['Seurat', 'Weinreb17', 'Zheng17']))
            subcol1, _, _ = st.columns(3)
            submit_btn = subcol1.form_submit_button(label='Apply', use_container_width=True)
            
            if submit_btn:
                if recipe == 'Seurat':
                    sc.pp.recipe_seurat(self.adata)
                elif recipe == 'Weinreb17':
                    sc.pp.recipe_weinreb17(self.adata)
                elif recipe == 'Zheng17':
                    sc.pp.recipe_zheng17(self.adata)
                else:
                    st.error("Recipe not found")

                #make adata
                st.session_state.adata_state.add_adata(AdataModel(
                    work_id=st.session_state.current_workspace.id, adata=self.adata, 
                    filename=os.path.join(os.getenv('WORKDIR'), "adata", "adata_pp.h5ad"), adata_name="adata_pp"))
                st.session_state.adata_state.switch_adata(adata_name="adata_pp")
                st.toast(f"Applied recipe: {st.session_state.sb_pp_recipe}", icon='✅')

    
    def annotate_mito(self):
        with st.form(key="form_annotate_mito"):
            st.subheader("Annotate Mitochondrial Genes", help="Filter mitochrondrial gene counts. All mitochrondrial genes \
                        are by default annotated and placed in the 'mt' variable.")
            
            self.adata.var['mt'] = self.adata.var_names.str.startswith(('MT-', 'mt-'))
            sc.pp.calculate_qc_metrics(self.adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
                
            st.text(f"Found {self.adata.var.mt.sum()} mitochondrial genes")
            
            subcol1, subcol2, subcol3 = st.columns(3, gap="small")
            with subcol1:
                ax_scatter = sc.pl.scatter(self.adata, x='total_counts', y='pct_counts_mt')
                with st.expander(label="Scatter"):
                    st.pyplot(ax_scatter)

            with subcol2:
                ax_violin = sc.pl.violin(self.adata, 'pct_counts_mt')
                with st.expander(label="Violin"):
                    st.pyplot(ax_violin)
                    

            max_pct_counts_mt = st.number_input(label="max pct_counts_mt", key="ni_pct_counts_mt", min_value=0)

            subcol1, _, _ = st.columns(3)
            mito_annot_submit_btn = subcol1.form_submit_button(label="Apply", use_container_width=True)

            if mito_annot_submit_btn:
                self.adata = self.adata[self.adata.obs.pct_counts_mt < max_pct_counts_mt, :]
                #make adata
                st.session_state.adata_state.add_adata(AdataModel(
                    work_id=st.session_state.current_workspace.id, adata=self.adata, 
                    filename=os.path.join(os.getenv('WORKDIR'), "adata", "adata_pp.h5ad"), adata_name="adata_pp"))
                st.session_state.adata_state.switch_adata(adata_name="adata_pp")

            

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

                subcol1, subcol2, subcol3 = st.columns(3, gap="small")
                with subcol1:
                    ax_scatter = sc.pl.scatter(self.adata, x='total_counts', y='pct_counts_ribo')
                    with st.expander(label="Scatter"):
                        st.pyplot(ax_scatter)

                with subcol2:
                    ax_violin = sc.pl.violin(self.adata, 'pct_counts_ribo')
                    with st.expander(label="Violin"):
                        st.pyplot(ax_violin)
                        

            max_pct_counts_ribo = st.number_input(label="max pct_counts_ribo", key="ni_pct_counts_ribo", min_value=0)

            subcol1, _, _ = st.columns(3)
            ribo_annot_submit_btn = subcol1.form_submit_button(label="Apply", use_container_width=True)

            if ribo_annot_submit_btn:
                self.adata = self.adata[self.adata.obs.pct_counts_ribo < max_pct_counts_ribo, :]
                #make adata
                st.session_state.adata_state.add_adata(AdataModel(
                    work_id=st.session_state.current_workspace.id, adata=self.adata, 
                    filename=os.path.join(os.getenv('WORKDIR'), "adata", "adata_pp.h5ad"), adata_name="adata_pp"))
                st.session_state.adata_state.switch_adata(adata_name="adata_pp")
                
        

    def run_scrublet(self):
        with st.form(key="scrublet_form"):
            st.subheader("Doublet Prediction")
            st.write("Use Scrublet to remove cells predicted to be doublets.")
            sim_doublet_ratio = st.number_input(label="Sim doublet ratio", value=2, key="ni_sim_doublet_ratio")
            expected_doublet_rate = st.number_input(label="Expected doublet rate", value=0.05, key="ni_expected_doublet_rate")
            scrublet_submit = st.form_submit_button(label="Filter")

            if scrublet_submit:
                with st.spinner("Running scrublet"):
                    adata_scrublet = sc.external.pp.scrublet(self.adata, sim_doublet_ratio=sim_doublet_ratio, expected_doublet_rate=expected_doublet_rate)
                    self.adata = adata_scrublet #TODO: only temporary, change to saving separately
                    #make adata
                    st.session_state.adata_state.add_adata(AdataModel(
                        work_id=st.session_state.current_workspace.id, adata=self.adata, 
                        filename=os.path.join(os.getenv('WORKDIR'), "adata", "adata_pp.h5ad"), adata_name="adata_pp"))
                    st.session_state.adata_state.switch_adata(adata_name="adata_pp")
                    
                    
    def regress_out(self):
        def regress_out_btn():
            if st.session_state.ms_regress_out_keys:
                sc.pp.regress_out(self.adata, keys=st.session_state.ms_regress_out_keys)
                st.toast("Successfully regressed out data", icon="✅")
            else:
                st.toast("No option selected, not regressing data.", icon="ℹ️")
        with st.form(key="regress_out_form"):
            st.subheader("Regress out", help="Regress out (mostly) unwanted sources of variation. Uses simple linear regression. This is inspired by Seurat's regressOut function in R [Satija15]. Note that this function tends to overcorrect in certain circumstances as described in :issue:526.")
            regress_keys = st.multiselect(label="Keys", options=self.adata.obs_keys(), key="ms_regress_out_keys")
            regress_out_btn = st.form_submit_button(label="Regress out", on_click=regress_out_btn)
            
            
    def scale_to_unit_variance(self):
        with st.form(key="scale_to_unit_variance_form"):
            st.subheader("Scale to unit variance")
            st.number_input(label="Max value", value=10, key="ni_scale_data_max_value")
            btn_scale_data_btn = st.form_submit_button(label="Apply")
            if btn_scale_data_btn:
                if st.session_state.ni_scale_data_max_value:
                    sc.pp.scale(self.adata, max_value=st.session_state.ni_scale_data_max_value)
                    st.toast("Successfully scaled data", icon="✅")
                else:
                    st.toast("Max value cannot be blank", icon="❌")
            
    
    def sample_data(self):
        st.subheader("Sample data")
        downsample_tab, subsample_tab = st.tabs(['Downsample', 'Subsample'])
        
        with downsample_tab:
            with st.form(key="downsample_form"):
                st.subheader("Downsample counts")
                counts_per_cell = st.number_input(label="Counts per cell")
                total_counts = st.number_input(label="Total counts")
                btn_downsample = st.form_submit_button(label="Apply")
                if btn_downsample:
                    sc.pp.downsample_counts(self.adata, counts_per_cell=counts_per_cell, total_counts=total_counts)
            
        with subsample_tab:
            with st.form(key="subsample_form"):
                st.subheader("Subsample counts")
                n_obs = st.number_input(label="n obs")
                fraction = st.number_input(label="subsample_fraction")
                btn_subsample = st.form_submit_button(label="Apply")
                if btn_subsample:
                    sc.pp.subsample(self.adata, n_obs=st.session_state.ni_n_obs, fraction=st.session_state.ni_subsample_fraction)
                    
                    
    def batch_effect_removal(self):
        with st.form(key="batch_effect_removal_form"):
            st.subheader("Batch effect removal")
            key = st.text_input(label="Key")
            covariates = st.multiselect(placeholder="Optional", label="Covariates", options=self.adata.obs_keys())
            btn_batch_effect_removal = st.form_submit_button(label="Apply")
            if btn_batch_effect_removal:
                sc.pp.combat(self.adata, key=key, covariates=covariates)
        
            

try:
    #make adata
    st.session_state.adata_state.add_adata(AdataModel(
        work_id=st.session_state.current_workspace.id, adata=st.session_state.adata_state.current.adata.copy(), 
        filename=os.path.join(os.getenv('WORKDIR'), "adata", "adata_pp.h5ad"), adata_name="adata_pp"))
    st.session_state.adata_state.switch_adata(adata_name="adata_pp")
    
    sidebar = Sidebar()

    sidebar.show()

    st.session_state["current_adata"] = st.session_state.adata_state.current.adata
    preprocess = Preprocess(st.session_state.adata_state.current.adata.copy())

    col1, col2, col3 = st.columns(3, gap="medium")

    with col1:
        preprocess.filter_highest_expr_genes()
        preprocess.filter_highly_variable_genes()
        preprocess.normalize_counts()
        preprocess.regress_out()
        preprocess.scale_to_unit_variance()

    with col2:
        preprocess.recipes()
        preprocess.filter_cells()
        preprocess.filter_genes()
        preprocess.sample_data()
        
            
    with col3:
        preprocess.annotate_mito()
        preprocess.annotate_ribo()
        preprocess.run_scrublet()
        preprocess.batch_effect_removal()

    sidebar.show_preview()
    sidebar.delete_experiment_btn()

except KeyError as ke:
    print("KeyError: ", ke)
    st.error("Couldn't find adata object in session, have you uploaded one?")
except Exception as e:
    print('Error: ', e)
    st.error(e)



