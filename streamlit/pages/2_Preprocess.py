import streamlit as st
import scanpy as sc
import pickle
import pandas as pd
import warnings
import os

warnings.filterwarnings("ignore")

sc.settings.verbosity = 3

st.set_page_config(layout="wide", page_title='Nuwa', page_icon='🧬')

common_style = """
            <style>
            footer {visibility: hidden;}
            .st-emotion-cache-1cypcdb {background: linear-gradient(180deg, rgb(5, 39, 103) 0%, #3a0647 70%); box-shadow: 1px 0 10px -2px #000;}
            .st-emotion-cache-86cver {color: rgba(250, 250, 250, 0.6)}
            .stButton button {
                border-radius: 0.5rem;
                background: #004dcf;
                color: #fff;
                border: 1px solid #004dcf;
                padding: 0.25rem 0.75rem;
            }

            .stButton button:hover {
                background: transparent;
                color: #004dcf;
                transition: all 0.1s ease-in-out;
                border: 1px solid #004dcf;
            }

            .stButton button:focus {
                border: 1px solid #004dcf;
                color: #004dcf;
            }

            .stButton button:active {
                border: 1px solid #004dcf;
                color: #004dcf;
            }

            .stButton button:visited {
                border: 1px solid #004dcf;
                color: #004dcf;
            }

            .st-emotion-cache-1b9yna5:focus:not(:active) {
                border: 1px solid #004dcf;
                color: #004dcf;
            }

            .st-emotion-cache-1b9yna5:focus:not(:hover) {
                border: 1px solid #004dcf;
                color: #fff;
            }

            .st-emotion-cache-oooxyj:focus:not(:active) {
                border: 1px solid #004dcf;
                color: #004dcf;
            }

            .st-emotion-cache-oooxyj:focus:not(:hover) {
                border: 1px solid #004dcf;
                color: #fff;
            }

            .st-emotion-cache-1ts31n5:focus:not(:active) {
                border: 1px solid #004dcf;
                color: #004dcf;
            }

            .st-emotion-cache-1ts31n5:focus:not(:hover) {
                border: 1px solid #004dcf;
                color: #fff;
            }
            </style>
            """
st.markdown(common_style, unsafe_allow_html=True)



st.set_option('deprecation.showPyplotGlobalUse', False)

try:
    if 'adata' not in st.session_state:
        tmp_file = open("./tmp/adata.pkl",'rb')
        cached_adata = pickle.load(tmp_file)
        st.session_state["adata"] = cached_adata
except:
    print("There was an error")


adata = st.session_state["adata"]


class Preprocess:
    def __init__(self, adata):
        self.adata = adata



        st.title("Preprocess")

    def show_preview(self):
        preview = st.container()
        preview.subheader("File preview")
        preview.code(st.session_state.adata)

    def filter_highest_expr_genes(self):
        st.subheader("Show highest expressed genes")
        st.number_input(label="Number of genes", min_value=1, max_value=100, value=20, key="n_top_genes")
        fn = 'figures/highest_expr_genes.pdf'
        sc.pl.highest_expr_genes(self.adata, n_top=st.session_state.n_top_genes, save=True)
        
        with open(fn, "rb") as img:
            ax = sc.pl.highest_expr_genes(self.adata, n_top=st.session_state.n_top_genes)
            with st.expander(label="Show figure"):
                st.pyplot(ax)

            btn = st.download_button (label="Download image", data=img, file_name=fn, mime="application/pdf", use_container_width=True)


    def filter_highly_variable_genes(self):
        self.adata_hvg = self.adata
        st.subheader("Filter highly variable genes")
        st.number_input(label="min mean", value=0.0125, key="input_highly_variable_min_mean")
        st.number_input(label="max mean", value=3.0, key="input_highly_variable_max_mean")
        sc.pp.normalize_total(self.adata_hvg, target_sum=1e4)
        sc.pp.log1p(self.adata_hvg)
        sc.pp.highly_variable_genes(self.adata_hvg, min_mean=0.0125, max_mean=3, min_disp=0.5)
        fn = 'figures/filter_genes_dispersion.pdf'
        
        sc.pl.highly_variable_genes(self.adata_hvg, save=True)

        with open(fn, "rb") as img:
            ax = sc.pl.highly_variable_genes(self.adata_hvg)
            with st.expander(label="Show figure"):
                st.pyplot(ax)
            subcol1, subcol2 = st.columns(2)
            with subcol1:
                def filter_variable_callback():
                    st.session_state.adata = self.adata_hvg
                    self.adata = self.adata_hvg
                    st.toast(body="Filtered variable genes")

                btn = st.button(label="Apply filter", key="btn_filter_highly_variable", use_container_width=True, on_click=filter_variable_callback)
            with subcol2:
                btn = st.download_button (label="Download image", data=img, file_name=fn, mime="application/pdf")


    def filter_cells(self):
        st.subheader("Filter Cells")
        subcol1, subcol2 = st.columns(2)
        self.adata_filter_cells = self.adata

        def filter_cells_callback():
            #filter min count
            if st.session_state.filter_cell_min_count:
                sc.pp.filter_cells(self.adata_filter_cells, min_counts=st.session_state.filter_cell_min_count)
            if st.session_state.filter_cell_max_count:
                sc.pp.filter_cells(self.adata_filter_cells, max_counts=st.session_state.filter_cell_max_count)
            if st.session_state.filter_cell_min_genes:
                sc.pp.filter_cells(self.adata_filter_cells, min_genes=st.session_state.filter_cell_min_genes)
            if st.session_state.filter_cell_max_genes:
                sc.pp.filter_cells(self.adata_filter_cells, max_genes=st.session_state.filter_cell_max_genes)

            self.adata = self.adata_filter_cells
            st.session_state.adata = self.adata_filter_cells
            st.toast("Filtered cells")

        with subcol1:
            st.number_input(label="Min count", min_value=1, value=None, key="filter_cell_min_count")
            st.number_input(label="min genes for cell", min_value=1, value=None, key="filter_cell_min_genes")

        with subcol2:
            st.number_input(label="Max count", min_value=1, value=None, key="filter_cell_max_count")
            st.number_input(label="max genes for cell", min_value=1, value=None, key="filter_cell_max_genes")
        
        btn = st.button(label="Apply filter", key="btn_filter_cells", on_click=filter_cells_callback)


    def filter_genes(self):
        st.subheader("Filter Genes")
        self.adata_filter_genes = self.adata

        def filter_genes_callback():
            if st.session_state.filter_gene_min_count:
                sc.pp.filter_genes(self.adata_filter_genes, min_counts=st.session_state.filter_gene_min_count)
            if st.session_state.filter_gene_max_count:
                sc.pp.filter_genes(self.adata_filter_genes, max_counts=st.session_state.filter_gene_max_count)
            if st.session_state.filter_gene_min_cells:
                sc.pp.filter_genes(self.adata_filter_genes, min_cells=st.session_state.filter_gene_min_cells)
            if st.session_state.filter_gene_max_cells:
                sc.pp.filter_genes(self.adata_filter_genes, max_cells=st.session_state.filter_gene_max_cells)

            self.adata = self.adata_filter_genes
            st.session_state.adata = self.adata_filter_genes
            st.toast("Filtered genes")

        subcol1, subcol2 = st.columns(2)
        with subcol1:
            st.number_input(label="Min count", min_value=1, value=None, key="filter_gene_min_count")
            st.number_input(label="min cells for gene", min_value=1, value=None, key="filter_gene_min_cells")

        with subcol2:
            st.number_input(label="Max count", min_value=1, value=None, key="filter_gene_max_count")
            st.number_input(label="max cells for gene", min_value=1, value=None, key="filter_gene_max_cells")
        
        btn = st.button(label="Apply filter", key="btn_filter_genes", on_click=filter_genes_callback)
    
    
    def annotate_mito(self):
        st.subheader("Annotate Mitochondrial Genes")
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

        def filter_mt_count():
            self.adata = self.adata[self.adata.obs.pct_counts_mt < st.session_state.ni_pct_counts_mt, :]
            st.session_state["adata"] = self.adata

        st.number_input(label="max pct_counts_mt", key="ni_pct_counts_mt", min_value=0)
        st.button(label="Apply filter", key="btn_annotate_mt_filter", on_click=filter_mt_count)

        st.session_state["adata"] = self.adata

    def annotate_ribo(self):
        st.subheader("Annotate Ribosomal Genes")
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

        def filter_ribo_count():
            self.adata = self.adata[self.adata.obs.pct_counts_ribo < st.session_state.ni_pct_counts_ribo, :]
            st.session_state["adata"] = self.adata

        st.number_input(label="max pct_counts_ribo", key="ni_pct_counts_ribo", min_value=0)
        st.button(label="Apply filter", key="btn_annotate_ribo_filter", on_click=filter_ribo_count)

        st.session_state["adata"] = self.adata


def add_experiment():
    print("add")


with st.sidebar:
    st.selectbox(label="Current Experiment:", options=(["raw", "adata"]))
    st.button(label="Add experiment", on_click=add_experiment, use_container_width=True)

preprocess = Preprocess(adata)

preprocess.show_preview()

col1, col2, col3 = st.columns(3, gap="large")

with col1:
    preprocess.filter_highest_expr_genes()
    st.divider()
    preprocess.filter_highly_variable_genes()
    st.divider()

with col2:
    preprocess.filter_cells()
    st.divider()
    preprocess.filter_genes()
        
with col3:
    preprocess.annotate_mito()
    st.divider()
    preprocess.annotate_ribo()
    st.divider()



