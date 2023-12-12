import streamlit as st
import scanpy as sc
import pickle
import pandas as pd
import warnings
import numpy as np

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
        
        
        
    def make_adata(self, adata):
        adata_name = f"{st.session_state.adata_state.current.adata_name}"
        st.session_state.adata_state.insert_record(AdataModel(
            work_id=st.session_state.current_workspace.id, adata=adata, 
            filename=os.path.join(os.getenv('WORKDIR'), "adata", f"{adata_name}.h5ad"), adata_name=adata_name))
        st.session_state.adata_state.switch_adata(adata_name=adata_name)
        
         

    def filter_highest_expr_genes(self):
        with st.form(key="form_highest_expr"):
            st.subheader("Show highest expressed genes")
            num_genes = st.number_input(label="Number of genes", min_value=1, max_value=100, value=20, key="n_top_genes")
            subcol1, _, _ = st.columns(3)
            submit_btn = subcol1.form_submit_button(label="Filter", use_container_width=True)

            if submit_btn:
                with st.spinner(text="Calculating highest expressed genes"):
                    ax = sc.pl.highest_expr_genes(self.adata, n_top=num_genes)
                    st.pyplot(ax)
                        
                        
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
                        sc.pp.normalize_total(self.adata, target_sum=1e4); st.session_state["script_state"].add_script("sc.pp.normalize_total(adata, target_sum=1e4)")
                        sc.pp.log1p(self.adata); st.session_state["script_state"].add_script("sc.pp.log1p(adata)")
                        sc.pp.highly_variable_genes(self.adata, min_mean=min_mean, max_mean=max_mean, min_disp=0.5); st.session_state["script_state"].add_script(f"sc.pp.highly_variable_genes(adata, min_mean={min_mean}, max_mean={max_mean}, min_disp=0.5)")
                        #make adata
                        self.make_adata(self.adata)
                        ax = sc.pl.highly_variable_genes(self.adata)
                        st.pyplot(ax)
                        #add to script state
                        st.session_state["script_state"].add_script("sc.pl.highly_variable_genes(adata)")
                        st.session_state["script_state"].add_script("plt.show()")
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
                    sc.pp.normalize_total(self.adata, target_sum=target_sum, exclude_highly_expressed=exclude_high_expr); st.session_state["script_state"].add_script(f"sc.pp.normalize_total(adata, target_sum={target_sum}, exclude_highly_expressed={exclude_high_expr})")
                    #make adata
                    self.make_adata(self.adata)
                    st.toast("Normalized data", icon='✅')
        with tab_per_cell:
            with st.form(key="form_normalize_per_cell"):
                counts_per_cell_after = st.number_input(label="Counts per cell after")

                subcol1, _, _ = st.columns(3)
                submit_btn = subcol1.form_submit_button(label="Apply", use_container_width=True)

                if submit_btn:
                    sc.pp.normalize_per_cell(self.adata, counts_per_cell_after=counts_per_cell_after); st.session_state["script_state"].add_script(f"sc.pp.normalize_per_cell(adata, counts_per_cell_after={counts_per_cell_after})")
                    #make adata
                    self.make_adata(self.adata)
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
                sc.pp.filter_cells(self.adata, max_genes=max_genes, min_genes=min_genes, max_counts=max_count, min_counts=min_count); st.session_state["script_state"].add_script(f"sc.pp.filter_cells(adata, max_genes={max_genes}, min_genes={min_genes}, max_counts={max_count}, min_counts={min_count})")
                #make adata
                self.make_adata(self.adata)
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
                sc.pp.filter_genes(self.adata, max_cells=max_cells, min_cells=min_cells, max_counts=max_count, min_counts=min_count); st.session_state["script_state"].add_script(f"sc.pp.filter_genes(adata, max_cells={max_cells}, min_cells={min_cells}, max_counts={max_count}, min_counts={min_count})")
                #make adata
                self.make_adata(self.adata)
                st.toast("Filtered genes", icon='✅')


    def recipes(self):
        with st.form(key="form_recipes"):
            st.subheader("Preprocess Recipes")
            recipe = st.selectbox(label="Recipe", key="sb_pp_recipe", options=(['Seurat', 'Weinreb17', 'Zheng17']))
            subcol1, _, _ = st.columns(3)
            submit_btn = subcol1.form_submit_button(label='Apply', use_container_width=True)
            
            if submit_btn:
                if recipe == 'Seurat':
                    sc.pp.recipe_seurat(self.adata); st.session_state["script_state"].add_script("sc.pp.recipe_seurat(adata)")
                elif recipe == 'Weinreb17':
                    sc.pp.recipe_weinreb17(self.adata); st.session_state["script_state"].add_script("sc.pp.recipe_weinreb17(adata)")
                elif recipe == 'Zheng17':
                    sc.pp.recipe_zheng17(self.adata); st.session_state["script_state"].add_script("sc.pp.recipe_zheng17(adata)")
                else:
                    st.error("Recipe not found")

                #make adata
                self.make_adata(self.adata)
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
                    

            ni_pct_counts_mt = st.number_input(label="max pct_counts_mt", key="ni_pct_counts_mt", min_value=0, value=100)

            subcol1, _, _ = st.columns(3)
            mito_annot_submit_btn = subcol1.form_submit_button(label="Apply", use_container_width=True)

            if mito_annot_submit_btn:
                self.adata = self.adata[self.adata.obs.pct_counts_mt < ni_pct_counts_mt, :]
                #write to script state
                st.session_state["script_state"].add_script("sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)")
                st.session_state["script_state"].add_script(f"adata = adata[adata.obs.pct_counts_mt < {ni_pct_counts_mt}, :]")
                #make adata
                self.make_adata(self.adata)
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

                subcol1, subcol2, subcol3 = st.columns(3, gap="small")
                with subcol1:
                    ax_scatter = sc.pl.scatter(self.adata, x='total_counts', y='pct_counts_ribo')
                    with st.expander(label="Scatter"):
                        st.pyplot(ax_scatter)

                with subcol2:
                    ax_violin = sc.pl.violin(self.adata, 'pct_counts_ribo')
                    with st.expander(label="Violin"):
                        st.pyplot(ax_violin)
                        

            max_pct_counts_ribo = st.number_input(label="max pct_counts_ribo", key="ni_pct_counts_ribo", min_value=0, value=100)

            ribo_annot_submit_btn = subcol1.form_submit_button(label="Apply", use_container_width=True)

            if ribo_annot_submit_btn:
                self.adata = self.adata[self.adata.obs.pct_counts_ribo < max_pct_counts_ribo, :]
                #add to script adata
                st.session_state["script_state"].add_script("sc.pp.calculate_qc_metrics(adata, qc_vars=['ribo'], percent_top=None, log1p=False, inplace=True)")
                st.session_state["script_state"].add_script(f"adata = adata[adata.obs.pct_counts_ribo < {max_pct_counts_ribo}, :]")
                #make adata
                self.make_adata(self.adata)
                st.toast("Filtered ribosomal genes", icon="✅")
                
                
    def annotate_hb(self):
        with st.form(key="form_annotate_hb"):
            st.subheader("Annotate haemoglobin genes", help="Filter haemoglobin gene counts. All haemoglobin genes \
                        are by default annotated and placed in the 'hb' variable.")
            
            # hemoglobin genes.
            self.adata.var['hb'] = self.adata.var_names.str.contains(("^HB[^(P)]"))
            
            sc.pp.calculate_qc_metrics(self.adata, qc_vars=['hb'], percent_top=None, log1p=False, inplace=True)
            
            st.text(f"Found {self.adata.var.hb.sum()} haemoglobin genes")

            subcol1, subcol2, _ = st.columns(3)
            
            with subcol1:
                ax_scatter = sc.pl.scatter(self.adata, x='total_counts', y='pct_counts_hb')
                with st.expander(label="Scatter"):
                    st.pyplot(ax_scatter)

            with subcol2:
                ax_violin = sc.pl.violin(self.adata, 'pct_counts_hb')
                with st.expander(label="Violin"):
                    st.pyplot(ax_violin)
                    
            max_pct_counts_hb = st.number_input(label="max pct_counts_hb", key="ni_pct_counts_hb", min_value=0, value=100)
                    
            submit_btn = subcol1.form_submit_button(label="Apply", use_container_width=True)
            
            if submit_btn:
                self.adata = self.adata[self.adata.obs.pct_counts_ribo < max_pct_counts_ribo, :]
                #add to script adata
                st.session_state["script_state"].add_script("sc.pp.calculate_qc_metrics(adata, qc_vars=['hb'], percent_top=None, log1p=False, inplace=True)")
                st.session_state["script_state"].add_script(f"adata = adata[adata.obs.pct_counts_hb < {max_pct_counts_hb}, :]")
                #make adata
                self.make_adata(self.adata)
                st.toast("Filtered ribosomal genes", icon="✅")
                
                

    def run_scrublet(self):
        with st.form(key="scrublet_form"):
            st.subheader("Doublet Prediction")
            st.write("Use Scrublet to remove cells predicted to be doublets.")
            sim_doublet_ratio = st.number_input(label="Sim doublet ratio", value=2, key="ni_sim_doublet_ratio")
            expected_doublet_rate = st.number_input(label="Expected doublet rate", value=0.05, key="ni_expected_doublet_rate")
            subcol1, _, _ = st.columns(3)
            scrublet_submit = subcol1.form_submit_button(label="Filter", use_container_width=True)

            if scrublet_submit:
                with st.spinner("Running scrublet"):
                    adata_scrublet = sc.external.pp.scrublet(self.adata, sim_doublet_ratio=sim_doublet_ratio, expected_doublet_rate=expected_doublet_rate); st.session_state["script_state"].add_script(f"adata_scrublet = sc.external.pp.scrublet(adata, sim_doublet_ratio={sim_doublet_ratio}, expected_doublet_rate={expected_doublet_rate})")
                    self.adata = adata_scrublet #TODO: only temporary, change to saving separately
                    #make adata
                    self.make_adata(self.adata)
                    st.toast("Completed doublet predictions", icon="✅")
                    
                    
    def regress_out(self):
        def regress_out_btn():
            if st.session_state.ms_regress_out_keys:
                sc.pp.regress_out(self.adata, keys=st.session_state.ms_regress_out_keys); st.session_state["script_state"].add_script(f"sc.pp.regress_out(adata, keys={st.session_state.ms_regress_out_keys})")
                st.toast("Successfully regressed out data", icon="✅")
                self.make_adata(self.adata)
            else:
                st.toast("No option selected, not regressing data.", icon="ℹ️")
        with st.form(key="regress_out_form"):
            st.subheader("Regress out", help="Regress out (mostly) unwanted sources of variation. Uses simple linear regression. This is inspired by Seurat's regressOut function in R [Satija15]. Note that this function tends to overcorrect in certain circumstances as described in :issue:526.")
            st.write("Uses linear regression to remove unwanted sources of variation.")
            regress_keys = st.multiselect(label="Keys", options=self.adata.obs_keys(), key="ms_regress_out_keys")
            subcol1, _, _ = st.columns(3)
            regress_out_btn = subcol1.form_submit_button(label="Apply", on_click=regress_out_btn, use_container_width=True)
            
            
    def scale_to_unit_variance(self):
        with st.form(key="scale_to_unit_variance_form"):
            st.subheader("Scale to unit variance")
            st.number_input(label="Max value", value=10, key="ni_scale_data_max_value")
            subcol1, _, _ = st.columns(3)
            btn_scale_data_btn = subcol1.form_submit_button(label="Apply", use_container_width=True)
            if btn_scale_data_btn:
                if st.session_state.ni_scale_data_max_value:
                    sc.pp.scale(self.adata, max_value=st.session_state.ni_scale_data_max_value); st.session_state["script_state"].add_script(f"sc.pp.scale(adata, max_value={st.session_state.ni_scale_data_max_value})")
                    self.make_adata(self.adata)
                    st.toast("Successfully scaled data", icon="✅")
                else:
                    st.toast("Max value cannot be blank", icon="❌")
            
    
    def sample_data(self):
        st.subheader("Sample data")
        downsample_tab, subsample_tab = st.tabs(['Downsample', 'Subsample'])
        
        with downsample_tab:
            with st.form(key="downsample_form"):
                st.subheader("Downsample counts")
                st.write("Downsample counts from count matrix.")
                counts_per_cell = st.number_input(label="Counts per cell")
                total_counts = st.number_input(label="Total counts")
                subcol1, _, _ = st.columns(3)
                btn_downsample = subcol1.form_submit_button(label="Apply", use_container_width=True)
                if btn_downsample:
                    sc.pp.downsample_counts(self.adata, counts_per_cell=counts_per_cell, total_counts=total_counts); st.session_state["script_state"].add_script(f"sc.pp.downsample_counts(adata, counts_per_cell={counts_per_cell}, total_counts={total_counts})")
                    self.make_adata(self.adata)
                    st.toast("Successfully downsampled data", icon="✅")
            
        with subsample_tab:
            with st.form(key="subsample_form"):
                st.subheader("Subsample counts")
                st.write("Subsample to a fraction of the number of observations.")
                n_obs = st.number_input(label="n obs")
                fraction = st.number_input(label="subsample_fraction")
                subcol1, _, _ = st.columns(3)
                btn_subsample = subcol1.form_submit_button(label="Apply", use_container_width=True)
                if btn_subsample:
                    sc.pp.subsample(self.adata, n_obs=st.session_state.ni_n_obs, fraction=st.session_state.ni_subsample_fraction); st.session_state["script_state"].add_script(f"sc.pp.subsample(adata, n_obs={st.session_state.ni_n_obs}, fraction={st.session_state.ni_subsample_fraction})")
                    self.make_adata(self.adata)
                    st.toast("Successfully subsampled data", icon="✅")
                    
    def batch_effect_removal(self):
        with st.form(key="batch_effect_removal_form"):
            st.subheader("Batch effect correction", help="Uses Combat to correct non-biological differences caused by batch effect.")
            key = st.selectbox(label="Key", options=self.adata.obs_keys(), key="sb_batch_effect_key")
            covariates = st.multiselect(placeholder="Optional", label="Covariates", options=self.adata.obs_keys())
            subcol1, _, _ = st.columns(3)
            btn_batch_effect_removal = subcol1.form_submit_button(label="Apply", use_container_width=True)
            if btn_batch_effect_removal:
                with st.spinner(text="Running Combat batch effect correction"):
                    sc.pp.combat(self.adata, key=key, covariates=covariates); st.session_state["script_state"].add_script(f"sc.pp.combat(adata, key={key}, covariates={covariates})")
                self.make_adata(self.adata)
                st.toast("Batch corrected data", icon='✅')
                
                
    def pca(self):
        with st.form(key="pca_pp_form"):
            st.subheader("PCA")
            st.write("Run PCA for guiding preprocessing.")
            
            def run_pca(adata):
                with st.spinner(text="Running PCA"):
                    sc.pp.pca(adata, random_state=42)
                    df_pca = pd.DataFrame({'pca1': adata.obsm['X_pca'][:,0], 'pca2': adata.obsm['X_pca'][:,1], 'color': adata.obs[f'{st.session_state.sb_pca_color_pp}']})  
                    pca_empty.empty()
                    pca_empty.scatter_chart(data=df_pca, x='pca1', y='pca2', color='color', size=18)
                    self.make_adata(self.adata)
                
               
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
                st.session_state["script_state"].add_script("sc.pp.pca(adata, random_state=42)")
                st.session_state["script_state"].add_script("sc.pl.pca(adata, random_state=42)")
                st.session_state["script_state"].add_script("plt.show()")
                run_pca(self.adata)
                
    
    def predict_sex(self):
        
        st.subheader("Measure gene counts")
        single_dataset, subsample = st.tabs(['Current dataset', 'Subsample'])
        
        with single_dataset:
            with st.form(key="sex_predict_single_dataset"):
                st.subheader("Collective counts across dataset")
                gene = st.selectbox(label="Gene", options=self.adata.var_names)
                subcol_btn1, _, _ = st.columns(3)
                submit_btn = subcol_btn1.form_submit_button(label="Run", use_container_width=True)
                if submit_btn:
                    with st.spinner(text="Locating sex genes"):
                        self.adata.obs["gene-counts"] = self.adata.X[:,self.adata.var_names.str.match(f'{gene}')].toarray()
                        arr = np.array([f'{st.session_state.adata_state.current.adata_name}'])
                        df = pd.DataFrame({f'{gene} count': self.adata.obs["gene-counts"], "Dataset": np.repeat(arr, self.adata.n_obs)})
                        st.bar_chart(data=df, x="Dataset", y=f'{gene} count')
        with subsample:
            with st.form(key="sex_predict_multiple_datasets"):
                st.subheader("Subsample counts in dataset")
                batch_key_sex_pred = st.selectbox(label="Obs key", options=self.adata.obs_keys(), key="sb_sex_pred_batch_key")
                gene = st.selectbox(label="Gene", options=self.adata.var_names)
                subcol_btn1, _, _ = st.columns(3)
                submit_btn = subcol_btn1.form_submit_button(label="Run", use_container_width=True)
                if submit_btn:
                    with st.spinner(text="Locating sex genes"):
                        self.adata.obs["gene-counts"] = self.adata.X[:,self.adata.var_names.str.match(f'{gene}')].toarray()
                        df = pd.DataFrame({f'{gene} count': self.adata.obs["gene-counts"], f"{batch_key_sex_pred}": self.adata.obs[f"{batch_key_sex_pred}"]})
                        st.bar_chart(data=df, x=f"{batch_key_sex_pred}", y=f'{gene} count')
                  

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
        preprocess.sample_data()
        

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



