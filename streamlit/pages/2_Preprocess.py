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
from time import sleep


warnings.filterwarnings("ignore")

sc.settings.verbosity = 3

st.set_page_config(layout="wide", page_title='Nuwa', page_icon='🧬')

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
        self.adata_raw = adata
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
                        #TODO:Add more params as input
                        sc.pp.highly_variable_genes(self.adata, min_mean=min_mean, max_mean=max_mean, min_disp=0.5)
                        ax = sc.pl.highly_variable_genes(self.adata)
                        st.pyplot(ax)

    def normalize_counts(self):
        st.subheader("Normalization")
        tab_total, tab_per_cell = st.tabs(['Total', 'Per cell'])

        with tab_total:
            with st.form(key="form_normalize_total"):
                target_sum = st.number_input(label="Target sum", value=1)
                exclude_high_expr = st.checkbox(label="Exclude highly expressed", value=False)

                subcol1, _, _ = st.columns(3)
                submit_btn = subcol1.form_submit_button(label="Apply", use_container_width=True)

                if submit_btn:
                    sc.pp.normalize_total(self.adata, target_sum=target_sum, exclude_highly_expressed=exclude_high_expr)
                    st.toast("Normalized data")
        with tab_per_cell:
            with st.form(key="form_normalize_per_cell"):
                counts_per_cell_after = st.number_input(label="Counts per cell after")

                subcol1, _, _ = st.columns(3)
                submit_btn = subcol1.form_submit_button(label="Apply", use_container_width=True)

                if submit_btn:
                    sc.pp.normalize_per_cell(self.adata, counts_per_cell_after=counts_per_cell_after)
                    st.toast("Normalized data")

    def save_adata(self, name):
        self.save_adata_to_session(name)
        self.save_adata_to_db(name)

    def save_adata_to_session(self, name):
        for i, adata in enumerate(st.session_state.adata):
            if adata.adata_name == name:
                st.session_state.adata[i] = AdataModel(work_id=st.session_state.current_workspace.id, id=i, adata_name=name, filename=f"{name}.h5ad", adata=self.adata)
                return
        st.session_state.adata.append(AdataModel(work_id=st.session_state.current_workspace.id, id=len(st.session_state.adata), adata_name=name, filename=f"{name}.h5ad", adata=self.adata))
        

    def save_adata_to_db(self, name):
        try:
            new_adata = schemas.Adata(
                work_id=int(st.session_state.current_workspace.id),
                adata_name=f"{name}",
                filename=f"/streamlit-volume/{st.session_state.current_workspace.id}/{name}.h5ad",
                notes="noteeesss"
            )
            self.conn.add(new_adata)
            self.conn.commit()
            self.conn.refresh(new_adata)
        except Exception as e:
            print(e)

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
                self.save_adata(name="adata_pp")
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
                self.save_adata(name="adata_pp")
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

                self.save_adata(name="adata_pp")
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
            submit_btn = subcol1.form_submit_button(label="Apply", use_container_width=True)

            if submit_btn:
                self.adata = self.adata[self.adata.obs.pct_counts_mt < max_pct_counts_mt, :]

            

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
            submit_btn = subcol1.form_submit_button(label="Apply", use_container_width=True)

            if submit_btn:
                self.adata = self.adata[self.adata.obs.pct_counts_ribo < max_pct_counts_ribo, :]

    def run_scrublet(self):
        with st.form(key="scrublet_form"):
            st.subheader("Doublet Prediction")
            st.write("Use Scrublet to remove cells predicted to be doublets.")
            sim_doublet_ratio = st.number_input(label="Sim doublet ratio", value=2)
            expected_doublet_rate = st.number_input(label="Expected doublet rate", value=0.05)
            scrublet_submit = st.form_submit_button(label="Filter")

            if scrublet_submit:
                with st.spinner("Running scrublet"):
                    adata_scrublet = sc.external.pp.scrublet(self.adata, sim_doublet_ratio=sim_doublet_ratio, expected_doublet_rate=expected_doublet_rate)
                    self.adata = adata_scrublet #TODO: only temporary, change to saving separately
                self.save_adata(name="adata_scrublet")

            

try:
    adata_models: [AdataModel] = st.session_state["adata"]
    show_sidebar(adata_models)

    adata = get_adata(adataList=adata_models, name=st.session_state.sb_adata_selection).adata
    st.session_state["current_adata"] = adata
    preprocess = Preprocess(adata)

    col1, col2, col3 = st.columns(3, gap="medium")

    with col1:
        preprocess.filter_highest_expr_genes()
        preprocess.filter_highly_variable_genes()
        preprocess.normalize_counts()

    with col2:
        preprocess.recipes()
        preprocess.filter_cells()
        preprocess.filter_genes()
            
    with col3:
        preprocess.annotate_mito()
        preprocess.annotate_ribo()
        preprocess.run_scrublet()

    show_preview()

except KeyError as ke:
    print("KeyError: ", ke)
    st.error("Couldn't find adata object in session, have you uploaded one?")
except Exception as e:
    print('Error: ', e)
    st.error(e)



