import streamlit as st
import scanpy as sc
import pickle
import pandas as pd
import warnings
import numpy as np

from models.AdataModel import AdataModel
from utils.plotting import highest_expr_genes_box_plot, plot_doubletdetection_threshold_heatmap
from components.sidebar import *
from datetime import datetime

from database.database import SessionLocal
from sqlalchemy.orm import Session

from database.schemas import schemas
from utils.AdataState import AdataState
import doubletdetection
from time import sleep
import os
import plotly.figure_factory as ff
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
    """
    Apply preprocessing on raw data for more effective analysis and detecting biological signal.

    Notes
    -----
    .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/preprocess_page.png
    
    Attributes
    ----------
    conn: Session
        Connection to postgresql database.
    adata: Anndata
        Gene expression matrix. 
    """

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
        import scanpy as sc
        sc.pl.highest_expr_genes(adata, n_top=20)
        """
        with st.form(key="form_highest_expr"):
            st.subheader("Show highest expressed genes")
            n_top_genes = st.number_input(label="Number of genes", min_value=1, max_value=100, value=20, key="ni:pp:highly_variable:n_top_genes")
            subcol1, _, _ = st.columns(3)
            submit_btn = subcol1.form_submit_button(label="Filter", use_container_width=True)

            if submit_btn:
                with st.spinner(text="Calculating highest expressed genes"):
                    st.session_state["script_state"].add_script(f"\n#Show highest expr genes\nsc.pl.highest_expr_genes(adata, n_top=20)\n")
                    fig = highest_expr_genes_box_plot(self.adata, n_top=n_top_genes)
                    st.plotly_chart(fig)
                        
                        
    def remove_genes(self):
        """
        Remove a gene from the dataset. Preserves the complete var names in raw attribute.

        Parameters
        ----------
        remove_genes: List[str]
            A list of gene names to remove.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/remove_genes.png

        Example
        -------
        import scanpy as sc

        remove_genes_list = ['malat1']
        for gene in remove_genes_list:
            remove_genes = self.adata.var_names.str.startswith(gene)
            remove = np.array(remove_genes)
            keep = np.invert(remove)
            adata = self.adata[:,keep]
        """
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
        """
        Display genes showing highly variance, optionally remove genes that don't show high variance from the dataset (preserves full genes in raw attribute).

        Parameters
        ----------
        flavour: str
            The method for computing variability in genes. Split between dispersion-based methods [Satija15] and [Zheng17] and normalized variance [Stuart19].

        n_top_genes: int
            Number of top variable genes to display. Mandatory for flavour 'seurat_v3'.

        min_mean: float
            Minimum mean. Ignored if flavour = 'seurat_v3'.

        max_mean: float
            Maximum mean. Ignored if flavour = 'seurat_v3'.

        min_disp: float
            Minimum dispersion. Ignored if flavour = 'seurat_v3'.

        max_disp: float
            Maximum dispersion. Ignored if flavour = 'seurat_v3'.

        remove: bool
            If set to True, removes non-variable genes from dataset.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/highly_variable.png

        Example
        -------
        import scanpy as sc

        sc.pp.normalize_total(adata, target_sum=1e4); 
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, min_mean=min_mean, max_mean=max_mean, min_disp=0.5)
        sc.pl.highly_variable_genes(adata)
        """
        try:

            def run_highly_variable(flavour, min_mean=None, max_mean=None, min_disp=None, max_disp=None, n_top_genes=None):
                with st.spinner(text="Calculating highly variable genes"):
                    sc.pp.normalize_total(self.adata, target_sum=1e4); 
                    sc.pp.log1p(self.adata)     
                    sc.pp.highly_variable_genes(self.adata, flavor=flavour, n_top_genes=n_top_genes, min_mean=min_mean, max_mean=max_mean, min_disp=min_disp, max_disp=max_disp)

                    #write to script state
                    st.session_state["script_state"].add_script("#Filter highly variable genes")
                    st.session_state["script_state"].add_script("sc.pp.normalize_total(adata, target_sum=1e4)")
                    st.session_state["script_state"].add_script("sc.pp.log1p(adata)")
                    st.session_state["script_state"].add_script(f"sc.pp.highly_variable_genes(adata, min_mean={min_mean}, max_mean={max_mean}, min_disp=0.5)")
                        
                    if remove:
                        self.adata = self.adata[:, self.adata.var.highly_variable]
                        
                    #make adata
                    self.save_adata()
                    #plot data
                    dispersions_tab, dispersions_norm_tab = st.tabs(['Raw', 'normalized'])
                    with dispersions_tab:
                        dispersions_tab.scatter_chart(self.adata.var, x='means', y='dispersions', color='highly_variable', size=10)
                    with dispersions_norm_tab:
                        dispersions_norm_tab.scatter_chart(self.adata.var, x='means', y='dispersions_norm', color='highly_variable', size=10)

                    #add to script state
                    st.session_state["script_state"].add_script("sc.pl.highly_variable_genes(adata)")

            st.subheader("Highly variable genes")
            seurat_tab, cell_ranger_tab, seuratv3_tab = st.tabs(['Seurat', 'Cell ranger', 'Seurat v3'])

            with seurat_tab.form(key="form_highly_variable_seurat"):
                n_top_genes = st.number_input(label="n_top_genes", min_value=1, key="ni:pp:highly_variable:seurat_n_top_genes", step=1, value=2000)
                subcol1, subcol2 = st.columns(2)
                min_mean = subcol1.number_input(label="min mean", format='%.4f', value=0.0125, key="ni:pp:highly_variable:seurat_min_mean")
                max_mean = subcol2.number_input(label="max mean", format='%.4f', value=3.0, key="ni:pp:highly_variable:seurat_max_mean")
                min_disp = subcol1.number_input(label="min_disp", format='%.4f', value=0.5, key="ni:pp:highly_variable:seurat_min_disp")
                max_disp = subcol2.number_input(label="max_disp", format='%.4f', value=3.0, key="ni:pp:highly_variable:seurat_max_disp")
                remove = st.toggle(label="Remove non-variable genes", value=False, key="toggle:pp:highly_variable:seurat_remove", help="By default, highly variable genes are only annoted. This option will remove genes without highly variable expression.")
                subcol1, _, _ = st.columns(3)
                submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)

                if submit_btn:
                    run_highly_variable(flavour="seurat", n_top_genes=n_top_genes, min_mean=min_mean, max_mean=max_mean, min_disp=min_disp, max_disp=max_disp)

            with cell_ranger_tab.form(key="form_highly_variable_cellranger"):
                n_top_genes = st.number_input(label="n_top_genes", min_value=1, key="ni:pp:highly_variable:cellranger_n_top_genes", step=1, value=2000)
                subcol1, subcol2 = st.columns(2)
                min_mean = subcol1.number_input(label="min mean", format='%.4f', value=0.0125, key="ni:pp:highly_variable:cellranger_min_mean")
                max_mean = subcol2.number_input(label="max mean", format='%.4f', value=3.0, key="ni:pp:highly_variable:cellranger_max_mean")
                min_disp = subcol1.number_input(label="min_disp", format='%.4f', value=0.5, key="ni:pp:highly_variable:cellranger_min_disp")
                max_disp = subcol2.number_input(label="max_disp", format='%.4f', value=3.0, key="ni:pp:highly_variable:cellranger_max_disp")
                remove = st.toggle(label="Remove non-variable genes", value=False, key="toggle:pp:highly_variable:cellranger_remove", help="By default, highly variable genes are only annoted. This option will remove genes without highly variable expression.")
                subcol1, _, _ = st.columns(3)
                submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)

                if submit_btn:
                    run_highly_variable(flavour="cell_ranger", n_top_genes=n_top_genes, min_mean=min_mean, max_mean=max_mean, min_disp=min_disp, max_disp=max_disp)

            with seuratv3_tab.form(key="form_highly_variable_seuratv3"):
                n_top_genes = st.number_input(label="n_top_genes", min_value=1, step=1, key="ni:pp:highly_variable:seurat_v3_n_top_genes", value=2000)
                remove = st.toggle(label="Remove non-variable genes", value=False, help="By default, highly variable genes are only annoted. This option will remove genes without highly variable expression.")
                subcol1, _, _ = st.columns(3)
                submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)

                if submit_btn:
                    run_highly_variable(flavour="seurat_v3", n_top_genes=n_top_genes)

        except Exception as e:
            st.toast(f"Failed to normalize data: {e}", icon="❌")

    def normalize_counts(self):
        """
        Normalize gene expression counts by setting the count per cell to a desired value.

        Parameters
        ----------
        target_sum: float
            The new gene counts will sum to this value. 

        exclude_high_expr: bool
            Exclude (very) highly expressed genes for the computation of the normalization factor (size factor) for each cell. 
            A gene is considered highly expressed, if it has more than max_fraction of the total counts in at least one cell. 
            The not-excluded genes will sum up to target_sum.

        log_tranform: bool
            Logarithmize data after normalization.

        max_fraction: float
            If exclude_highly_expressed=True, consider cells as highly expressed that have more counts than max_fraction of the 
            original total counts in at least one cell.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/normalize_counts.png

        Example
        -------
        import scanpy as sc

        sc.pp.normalize_total(adata, target_sum=1, exclude_highly_expressed=False, max_fraction=0.05)
        """

        with st.form(key="form_normalize_total"):
            st.subheader("Normalization")
            subcol_input1, subcol_input2 = st.columns(2, gap="medium")
            target_sum = subcol_input1.number_input(label="Target sum", value=1.0, key="ni:pp:normalize_counts:target_sum")
            max_fraction = subcol_input2.number_input(label="Max fraction", key="ni:pp:normalize_counts:max_fraction", value=0.050, min_value=0.001, max_value=1.000)
            exclude_high_expr = subcol_input1.checkbox(label="Exclude highly_expr", value=False)
            log_transform_total = subcol_input2.checkbox(label="Log transform", value=False)
            subcol1, _, _ = st.columns(3)
            submit_btn = subcol1.form_submit_button(label="Apply", use_container_width=True)

            if submit_btn:
                sc.pp.normalize_total(self.adata, target_sum=target_sum, exclude_highly_expressed=exclude_high_expr, max_fraction=max_fraction)
                if log_transform_total:
                    sc.pp.log1p(self.adata)
                #write to script state
                st.session_state["script_state"].add_script(f"#Normalize counts (total)\nsc.pp.normalize_total(adata, target_sum={target_sum}, exclude_highly_expressed={exclude_high_expr})")
                #make adata
                self.save_adata()
                st.toast("Normalized data", icon='✅')

    def filter_cells(self):
        """
        Filter cell outliers based on counts and numbers of genes expressed.

        Parameters
        ----------
        min_counts: int
            Minimum number of counts required for a cell to pass filtering.

        max_counts: int
            Maximum number of counts required for a cell to pass filtering.

        min_genes: int
            Minimum number of genes expressed required for a cell to pass filtering.

        max_genes: int
            Maximum number of genes expressed required for a cell to pass filtering.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/filter_cells.png

        Example
        -------
        import scanpy as sc
        
        sc.pp.filter_cells(adata, max_genes=None, min_genes=200, max_counts=None, min_counts=None)
        """
        with st.form(key="form_filter_cells"):
            st.subheader("Filter Cells")
            subcol1, subcol2 = st.columns(2)
            with subcol1:
                min_counts = st.number_input(label="Min counts", min_value=1, value=None, key="filter_cell_min_counts")
                min_genes = st.number_input(label="min genes for cell", min_value=1, value=None, key="filter_cell_min_genes")

            with subcol2:
                max_counts = st.number_input(label="Max counts", min_value=1, value=None, key="filter_cell_max_counts")
                max_genes = st.number_input(label="max genes for cell", min_value=1, value=None, key="filter_cell_max_genes")
      
            subcol1, _, _ = st.columns(3)
            submit_btn = subcol1.form_submit_button(label="Apply", use_container_width=True)

            if submit_btn:
                sc.pp.filter_cells(self.adata, max_genes=max_genes, min_genes=min_genes, max_counts=max_counts, min_counts=min_counts)
                #write to script state
                st.session_state["script_state"].add_script(f"#Filter cells\nsc.pp.filter_cells(adata, max_genes={max_genes}, min_genes={min_genes}, max_counts={max_counts}, min_counts={min_counts})")
                #make adata
                self.save_adata()
                st.toast("Filtered cells", icon='✅')


    def filter_genes(self):
        """
        Filter genes based on number of cells or counts.

        Parameters
        ----------
        min_counts: int
            Minimum number of counts required for a gene to pass filtering.

        max_counts: int
            Maximum number of counts required for a gene to pass filtering.

        min_cells: int
            Minimum number of cells expressed required for a gene to pass filtering.

        max_cells: int
            Maximum number of cells expressed required for a gene to pass filtering.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/filter_genes.png

        Example
        -------
        import scanpy as sc
        
        sc.pp.filter_genes(adata, max_cells=None, min_cells=None, max_counts=None, min_counts=3)
        """
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
        """
        Apply a standard preprocessing recipe to the data.

        Parameters
        ----------
        recipe: str
            Type of preprocessing recipe to apply. Available options are Seurat, Weinrev17 and Zheng17.

        log: bool
            Apply log tranform to data. If data is already logarithmized, set this to False.

        n_top_genes: int
            Number of genes to keep. Used for Zheng17 recipe.

        mean_threshold: float
            mean value for threshold. Used for Weinreb17 recipe.

        cv_threshold: float
            cv value for threshold. Used for Weinreb17 recipe.

        n_pcs: int
            Number of principle components to include. Used for Weinreb17 recipe.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/recipes.png

        Example
        -------
        import scanpy as sc

        # Seurat recipe
        sc.pp.recipe_seurat(self.adata, log=True) 
        # Weinreb17 recipe
        sc.pp.recipe_weinreb17(self.adata, log=True, mean_threshold=0.01, cv_threshold=2, n_pcs=50)
        # Zheng17 recipe
        sc.pp.recipe_zheng17(self.adata, log=True, n_top_genes=1000)
        """
        st.subheader("Preprocess Recipes")
        seurat_tab, weinreb17_tab, zheng17_tab = st.tabs(['Seurat', 'Weinreb17', 'Zheng17'])
        with seurat_tab:
            with st.form(key="form_seurat"):
                st.write("Parameters")
                log = st.checkbox(label="Log", value=True)
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


    def annotate(self):
        st.subheader("Annotate")
        mito, ribo, haem = st.tabs(['Mitochondrial', 'Ribosomal', 'Haemoglobin'])
        with mito:
            self.annotate_mito()
        with ribo:
            self.annotate_ribo()
        with haem:
            self.annotate_hb()

    
    def annotate_mito(self):
        """
        Annotate mitochondrial genes with an optional key to group by observations. Optionally remove cells that contain a high percentage of mitochondrial genes as these are likely low quality samples.

        Parameters
        ----------
        color_key: str
            Group the mitochondia annotations by this value.

        max_pct_counts_mt: int
            Maximum percentage of mitochondial genes to other genes. If a cell contains a mitochondrial gene percentage greater than this value, it will be removed.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/mito1.png
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/mito2.png

        Example
        -------
        import scanpy as sc
        """
        with st.form(key="form_annotate_mito"):
            st.subheader("Mitochondrial Genes", help="Filter mitochrondrial gene counts. All mitochrondrial genes \
                        are by default annotated and placed in the 'mt' variable.")
            
            self.adata.var['mt'] = self.adata.var_names.str.startswith(('MT-', 'mt-'))
            sc.pp.calculate_qc_metrics(self.adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
                
            st.text(f"Found {self.adata.var.mt.sum()} mitochondrial genes")

            subcol_input1, subcol_input2 = st.columns(2)
            options = np.append('None', st.session_state.adata_state.current.adata.obs_keys())
            color_key_mito = subcol_input1.selectbox(label="Color key", options=options)
            ni_pct_counts_mt = subcol_input2.number_input(label="max pct_counts_mt", key="ni_pct_counts_mt", min_value=0, value=100)

            scatter_chart, violin_chart = st.tabs(['Scatter', 'Violin'])
            mito_container_scatter = scatter_chart.empty()
            mito_container_violin = violin_chart.empty()

            def plot_charts(color=None):
                
                if color == 'None':
                    color=None
                if color != None:
                        color = self.adata.obs[color]

                with scatter_chart:
                    df_scatter = pd.DataFrame({'total_counts': self.adata.obs.total_counts, 'pct_counts_mt': self.adata.obs.pct_counts_mt, 'color': color})
                    mito_container_scatter.scatter_chart(df_scatter, x='total_counts', y='pct_counts_mt', color='color', size=10)

                with violin_chart:

                    # Plot violin plot
                    fig = go.Figure()

                    fig.add_trace(go.Violin(x=color, 
                        y=self.adata.obs.pct_counts_mt,
                        jitter=0.1, line_color='blue')
                    )

                    fig.update_traces(meanline_visible=True)
                    fig.update_layout(violingap=0, violinmode='group', xaxis_title="pct_counts_mt", yaxis_title="value") #add legend title
                    mito_container_violin.plotly_chart(fig, use_container_width=True)


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
        """
        Annotate ribosomal genes with an optional key to group by observations. Optionally remove cells that contain a high percentage of ribosomal genes as these are likely low quality samples.

        Parameters
        ----------
        color_key: str
            Group the ribosomal gene annotations by this value.

        max_pct_counts_ribo: int
            Maximum percentage of ribosomal genes to other genes. If a cell contains a ribosomal gene percentage greater than this value, it will be removed.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/ribo1.png
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/ribo2.png

        Example
        -------
        import scanpy as sc
        """
        with st.form(key="form_annotate_ribo"):
            st.subheader("Ribosomal Genes", help="Filter ribosomal gene counts. All ribosomal genes \
                        are by default annotated and placed in the 'ribo' variable.")

            with st.spinner(text="Fetching ribosomal genes"):
                ribo_url = "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"
                ribo_genes = pd.read_table(ribo_url, skiprows=2, header=None)

                self.adata.var['ribo'] = self.adata.var_names.isin(ribo_genes[0].values)
                sc.pp.calculate_qc_metrics(self.adata, qc_vars=['ribo'], percent_top=None, log1p=False, inplace=True)
                
                st.text(f"Found {self.adata.var['ribo'].sum()} ribosomal genes")

                subcol_input1, subcol_input2 = st.columns(2)
                options = np.append('None', st.session_state.adata_state.current.adata.obs_keys())
                color_key_ribo = subcol_input1.selectbox(label="Color key", options=options)
                ni_pct_counts_ribo = subcol_input2.number_input(label="max pct_counts_ribo", key="ni_pct_counts_ribo", min_value=0, value=100)

                scatter_chart, violin_chart = st.tabs(['Scatter', 'Violin'])
                ribo_container_scatter = scatter_chart.empty()
                ribo_container_violin = violin_chart.empty()

            def plot_charts(color=None):
                
                if color == 'None':
                    color=None
                if color != None:
                        color = self.adata.obs[color]

                with scatter_chart:
                    df_scatter = pd.DataFrame({'total_counts': self.adata.obs.total_counts, 'pct_counts_ribo': self.adata.obs.pct_counts_ribo, 'color': color})
                    ribo_container_scatter.scatter_chart(df_scatter, x='total_counts', y='pct_counts_ribo', color='color', size=10)

                with violin_chart:

                    # Plot violin plot
                    fig = go.Figure()

                    fig.add_trace(go.Violin(x=color, 
                        y=self.adata.obs.pct_counts_ribo,
                        jitter=0.1, line_color='blue')
                    )

                    fig.update_traces(meanline_visible=True)
                    fig.update_layout(violingap=0, violinmode='group', xaxis_title="pct_counts_ribo", yaxis_title="value") #add legend title
                    ribo_container_violin.plotly_chart(fig, use_container_width=True)


            plot_charts()

            subcol1, _, _ = st.columns(3)
            ribo_annot_submit_btn = subcol1.form_submit_button(label="Apply", use_container_width=True)

            if ribo_annot_submit_btn:
                plot_charts(color_key_ribo)
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
        """
        Annotate haemoglobin genes with an optional key to group by observations. Optionally remove cells that contain a high percentage of haemoglobin genes as these are likely low quality samples.

        Parameters
        ----------
        color_key: str
            Group the hb gene annotations by this value.

        max_pct_counts_hb: int
            Maximum percentage of hb genes to other genes. If a cell contains a hb gene percentage greater than this value, it will be removed.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/hb1.png
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/hb2.png

        Example
        -------
        import scanpy as sc


        """
        with st.form(key="form_annotate_hb"):
            st.subheader("Haemoglobin genes", help="Filter haemoglobin gene counts. All haemoglobin genes \
                        are by default annotated and placed in the 'hb' variable.")
            
            # hemoglobin genes.
            self.adata.var['hb'] = self.adata.var_names.str.contains(("^HB[^(P)]"))
            
            sc.pp.calculate_qc_metrics(self.adata, qc_vars=['hb'], percent_top=None, log1p=False, inplace=True)
            
            st.text(f"Found {self.adata.var.hb.sum()} haemoglobin genes")

            subcol_input1, subcol_input2 = st.columns(2)
            options = np.append('None', st.session_state.adata_state.current.adata.obs_keys())
            color_key_hb = subcol_input1.selectbox(label="Color key", options=options)
            ni_pct_counts_hb = subcol_input2.number_input(label="max pct_counts_hb", key="ni_pct_counts_hb", min_value=0, value=100)

            scatter_chart, violin_chart = st.tabs(['Scatter', 'Violin'])
            hb_container_scatter = scatter_chart.empty()
            hb_container_violin = violin_chart.empty()

            def plot_charts(color=None):
                
                if color == 'None':
                    color=None
                if color != None:
                        color = self.adata.obs[color]

                with scatter_chart:
                    df_scatter = pd.DataFrame({'total_counts': self.adata.obs.total_counts, 'pct_counts_hb': self.adata.obs.pct_counts_hb, 'color': color})
                    hb_container_scatter.scatter_chart(df_scatter, x='total_counts', y='pct_counts_hb', color='color', size=10)

                with violin_chart:

                    # Plot violin plot
                    fig = go.Figure()

                    fig.add_trace(go.Violin(x=color, 
                        y=self.adata.obs.pct_counts_hb,
                        jitter=0.1, line_color='blue')
                    )

                    fig.update_traces(meanline_visible=True)
                    fig.update_layout(violingap=0, violinmode='group', xaxis_title="pct_counts_hb", yaxis_title="value") #add legend title
                    hb_container_violin.plotly_chart(fig, use_container_width=True)


            plot_charts()

            subcol1, _, _ = st.columns(3)
            hb_annot_submit_btn = subcol1.form_submit_button(label="Apply", use_container_width=True)

            if hb_annot_submit_btn:
                plot_charts(color_key_hb)
                self.adata = self.adata[self.adata.obs.pct_counts_hb < ni_pct_counts_hb, :]
                #add to script adata
                st.session_state["script_state"].add_script("#Filter haemoglobin genes")
                st.session_state["script_state"].add_script("sc.pl.scatter(adata, x='total_counts', y='pct_counts_hb')")
                st.session_state["script_state"].add_script("sc.pl.violin(adata, 'pct_counts_hb')")
                st.session_state["script_state"].add_script("sc.pp.calculate_qc_metrics(adata, qc_vars=['hb'], percent_top=None, log1p=False, inplace=True)")
                st.session_state["script_state"].add_script(f"adata = adata[adata.obs.pct_counts_hb < {ni_pct_counts_hb}, :]")
                #make adata
                self.save_adata()
                st.toast("Filtered haemoglobin genes", icon="✅")


    def predict_doublets(self):
        st.subheader("Predict doublets")
        scrublet, doubletdetection = st.tabs(['scrublet', 'doubletdetection'])

        with scrublet:
            self.run_scrublet()

        with doubletdetection:
            self.run_doubletdetection()
                
                

    def run_scrublet(self):
        """
        Uses [Scrublet](https://github.com/swolock/scrublet?tab=readme-ov-file) to predict if an observation (cell) is likely to be a doublet and remove from the dataset. Doublets arise when multiple cells are mistaken as a single cell in droplet-based technologies. This affects biological signal, for example during PCA doublets may form separate clusters not reflecting biological difference.

        Parameters
        ----------
        sim_doublet_ratio: float
            Number of doublets to simulate relative to the number of observed transcriptomes.

        expected_doublet_rate: float
            Where adata_sim not suplied, the estimated doublet rate for the experiment.

        stdev_doublet_rate: float
            Where adata_sim not suplied, uncertainty in the expected doublet rate.

        batch_key: str
            Optional adata.obs column name discriminating between batches.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/scrublet.png

        Example
        -------
        !pip install scrublet # requires external library scrublet
        import scanpy as sc

        sc.external.pp.scrublet(adata, sim_doublet_ratio=2, expected_doublet_rate=0.05, stdev_doublet_rate=0.02, batch_key=None, random_state=42)

        # plot PCA with doublet predictions
        sc.pp.pca(adata)
        df = pd.DataFrame({'PCA 1': adata.obsm['X_pca'][:,0], 'PCA 2': adata.obsm['X_pca'][:,1]})
        scatter = plt.scatter(x=df['PCA 1'], y=df['PCA 2'], c=adata.obs.predicted_doublet, s=5, label=adata.obs.predicted_doublet.values)
        legend1 = plt.legend(*scatter.legend_elements(), loc="upper right", title="Predicted doublet")
        plt.ylabel('PCA 2')
        plt.xlabel('PCA 1')
        plt.title("Predicted doublets in PCA space")

        # Simulated doublets and observed transcriptomes prob density plots
        sc.external.pl.scrublet_score_distribution(adata)

        plt.show()
        """
        with st.form(key="scrublet_form"):
            st.subheader("Scrublet", help="Use Scrublet to remove cells predicted to be doublets.")
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
                    # plot PCA to see doublets
                    sc.external.pl.scrublet_score_distribution(self.adata)
                    sc.pp.pca(self.adata)
                    sc.pp.neighbors(self.adata)
                    sc.tl.umap(self.adata)
                    pca, umap, distplot, simulated_doublets = st.tabs(['PCA', 'UMAP', 'Distplot', 'Simulated doublets'])
                    with pca:
                        df = pd.DataFrame({'PCA 1': self.adata.obsm['X_pca'][:,0], 'PCA 2': self.adata.obsm['X_pca'][:,1], 'Predicted doublet': self.adata.obs.predicted_doublet})
                        st.scatter_chart(df, x='PCA 1', y='PCA 2', color='Predicted doublet', size=10)
                    with umap:
                        df = pd.DataFrame({'UMAP 1': self.adata.obsm['X_umap'][:,0], 'UMAP 2': self.adata.obsm['X_umap'][:,1], 'Predicted doublet': self.adata.obs.predicted_doublet})
                        st.scatter_chart(df, x='UMAP 1', y='UMAP 2', color='Predicted doublet', size=10)
                    with distplot:
                        fig = ff.create_distplot([self.adata.obs.doublet_score.values], group_labels=['doublet_score'], 
                            bin_size=0.02, show_rug=False, show_curve=False, colors=["#31abe8"])
                        fig.update_layout(yaxis_type="log") 
                        fig.add_vline(x=self.adata.uns['scrublet']['threshold'], line_color="red")
                        fig.update_layout(xaxis_title="Doublet score", yaxis_title="Probability density", title="Observed transcriptomes")
                        st.plotly_chart(fig, use_container_width=True)
                    with simulated_doublets:
                        st.write("To implement")
                    #write to script state
                    st.session_state["script_state"].add_script(f"#Detect doublets with scrublet\nadata_scrublet = sc.external.pp.scrublet(adata, sim_doublet_ratio={sim_doublet_ratio}, expected_doublet_rate={expected_doublet_rate})")
                    #make adata
                    self.save_adata()
                    st.toast("Run scrublet", icon="✅")


    def run_doubletdetection(self):
        """
        Uses [doubletdetection](https://zenodo.org/records/6349517) to predict if an observation (cell) is likely to be a doublet and remove from the dataset. Doublets arise when multiple cells are mistaken as a single cell in droplet-based technologies. This affects biological signal, for example during PCA doublets may form separate clusters not reflecting biological difference.
        
        Parameters
        ----------
        n_iters: int
            Number of fit operations from which to collect p-values. Defualt value is 25.

        pseudocount: float
            Pseudocount used in normalize_counts. If 1 is used, and standard_scaling=False, the classifier is much more memory efficient; however, this may result in fewer doublets detected.
        
        boost_rate: float
            Proportion of cell population size to produce as synthetic doublets.

        clustering_algorithm: str
            Clustering algorithm to use (Louvain, Leiden or phenograph).

        voter_thresh: float
            Fraction of iterations a cell must be called a doublet.

        standard_scaling: bool
            Set to True to enable standard scaling of normalized count matrix prior to PCA. Recommended when not using Phenograph.

        n_top_var_genes: int
            Number of highest variance genes to use; other genes discarded. Will use all genes when zero.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/doubletdetection.png

        Example
        -------
        # Example taken from doubletdetection tutorial for PBMC3K available at: https://doubletdetection.readthedocs.io/en/latest/tutorial.html
        !pip install doubletdetection # install doubletdetection
        import doubletdetection
        import scanpy as sc
        import matplotlib.pyplot as plt

        adata = sc.read_10x_h5(
            "pbmc_10k_v3_filtered_feature_bc_matrix.h5",
            backup_url="https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5"
        )
        adata.var_names_make_unique()

        sc.pp.filter_genes(adata, min_cells=1)

        clf = doubletdetection.BoostClassifier(
            n_iters=10,
            clustering_algorithm="louvain",
            standard_scaling=True,
            pseudocount=0.1,
            n_jobs=-1,
        )
        doublets = clf.fit(adata.X).predict(p_thresh=1e-16, voter_thresh=0.5)
        doublet_score = clf.doublet_score()

        adata.obs["doublet"] = doublets
        adata.obs["doublet_score"] = doublet_score

        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata)
        sc.tl.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)

        sc.pl.umap(adata, color=["doublet", "doublet_score"])
        sc.pl.violin(adata, "doublet_score")
        """
        with st.form(key="Doubletdetection_form"):
            st.subheader("Doubletdetection")
            col1, col2, col3 = st.columns(3)
            n_iters = col1.number_input(label="n_iters", min_value=1, step=1, value=10)
            pseudocount = col2.number_input(label="pseudocount", value=0.1, step=0.1)
            boost_rate = col3.number_input(label="boost rate", value=0.25, step=0.01)
            clustering_algorithm = col1.selectbox(label="algorithm", options=['leiden', 'louvain', 'phenograph'])
            voter_thresh = col2.number_input(label="p thresh", value=0.5, step=0.1)
            n_top_var_genes = col3.number_input(label="n_top_var_genes", value=10000, step=1)
            standard_scaling = st.toggle(label="standard scaling", value=False)
            subcol1, _, _ = st.columns(3)
            submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
            

            if "doublet_doubletdetection" in self.adata.obs:
                pca, umap, threshold  = st.tabs(['PCA', 'UMAP', 'Threshold'])
                with pca:
                    df = pd.DataFrame({'PCA 1': self.adata.obsm['X_pca'][:,0], 'PCA 2': self.adata.obsm['X_pca'][:,1], 'Predicted doublet': self.adata.obs.doublet_doubletdetection})
                    st.scatter_chart(df, x='PCA 1', y='PCA 2', color='Predicted doublet', size=10)
                with umap:
                    df = pd.DataFrame({'UMAP 1': self.adata.obsm['X_umap'][:,0], 'UMAP 2': self.adata.obsm['X_umap'][:,1], 'Predicted doublet': self.adata.obs.doublet_doubletdetection})
                    st.scatter_chart(df, x='UMAP 1', y='UMAP 2', color='Predicted doublet', size=10)
                with threshold:
                    fig_heatmap = plot_doubletdetection_threshold_heatmap(clf=st.session_state.doubletdetection_clf, height=400)
                    st.plotly_chart(fig_heatmap, use_container_width=True)

            if submit_btn:
                with st.spinner("Running doubletdetection"):
                    clf = doubletdetection.BoostClassifier(
                        n_iters=n_iters,
                        n_top_var_genes=n_top_var_genes,
                        clustering_algorithm=clustering_algorithm,
                        standard_scaling=standard_scaling,
                        pseudocount=pseudocount,
                        boost_rate=boost_rate,
                        n_jobs=-1,
                    )
                    st.session_state["doubletdetection_clf"] = clf
                    doublets = clf.fit(self.adata.X).predict(p_thresh=1e-16, voter_thresh=voter_thresh)
                    doublet_score = clf.doublet_score()

                    self.adata.obs["doublet_doubletdetection"] = doublets
                    self.adata.obs["doublet_score_doubletdetection"] = doublet_score

                    sc.pp.normalize_total(self.adata)
                    sc.pp.log1p(self.adata)
                    sc.pp.highly_variable_genes(self.adata)
                    sc.tl.pca(self.adata)
                    sc.pp.neighbors(self.adata)
                    sc.tl.umap(self.adata)

                    pca, umap, threshold  = st.tabs(['PCA', 'UMAP', 'Threshold'])
                    with pca:
                        df = pd.DataFrame({'PCA 1': self.adata.obsm['X_pca'][:,0], 'PCA 2': self.adata.obsm['X_pca'][:,1], 'Predicted doublet': self.adata.obs.doublet_doubletdetection})
                        st.scatter_chart(df, x='PCA 1', y='PCA 2', color='Predicted doublet', size=10)
                    with umap:
                        df = pd.DataFrame({'UMAP 1': self.adata.obsm['X_umap'][:,0], 'UMAP 2': self.adata.obsm['X_umap'][:,1], 'Predicted doublet': self.adata.obs.doublet_doubletdetection})
                        st.scatter_chart(df, x='UMAP 1', y='UMAP 2', color='Predicted doublet', size=10)
                    with threshold:
                        fig_heatmap = plot_doubletdetection_threshold_heatmap(clf=st.session_state.doubletdetection_clf, height=400)
                        st.plotly_chart(fig_heatmap, use_container_width=True)

                st.toast("Run doubletdetection", icon="✅")

                    
                    
    def regress_out(self):
        """
        Regress out unwanted sources of variation using linear regression. This is inspired by Seurat's regressOut function in R [Satija15]. Note that this function tends to overcorrect in certain circumstances.
        
        Parameters
        ----------
        keys: List[str]
            Keys for observation annotation on which to regress on.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/regress_out.png

        Example
        -------
        import scanpy as sc

        regress_keys = ['percent_mt', 'percent_ribo']
        sc.pp.regress_out(adata, keys=regress_keys)
        """
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
        """
        Scale data to unit variance and zero mean.

        Parameters
        ----------
        zero_center: bool
            If False, omit zero-centering variables, which allows to handle sparse input efficiently.

        max_value: float
            Clip to this value after scaling. If None, do not clip.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/scale.png

        Example
        -------
        import scanpy as sc

        sc.pp.scale(adata, zero_center=True, max_value=None)
        """
        with st.form(key="scale_to_unit_variance_form"):
            st.subheader("Scale to unit variance")
            zero_center = st.toggle(label="Zero center", value=True)
            max_value = st.number_input(label="Max value", value=10, key="ni_scale_data_max_value")
            subcol1, _, _ = st.columns(3)
            btn_scale_data_btn = subcol1.form_submit_button(label="Apply", use_container_width=True)
            if btn_scale_data_btn:
                if st.session_state.ni_scale_data_max_value:
                    sc.pp.scale(self.adata, zero_center=zero_center, max_value=max_value)
                    #write to script state
                    st.session_state["script_state"].add_script(f"#Scale to unit variance\nsc.pp.scale(adata, max_value={st.session_state.ni_scale_data_max_value})")
                    self.save_adata()
                    st.toast("Successfully scaled data", icon="✅")
                else:
                    st.toast("Max value cannot be blank", icon="❌")
            
    

    def downsample_data(self):
        """
        If counts_per_cell is specified, each cell will downsampled. If total_counts is specified, expression matrix will be downsampled to contain at most total_counts.

        Parameters
        ----------
        counts_per_cell: float
            Target total counts per cell. If a cell has more than 'counts_per_cell', it will be downsampled to this number. Resulting counts can be specified on a per cell basis by passing an array.Should be an integer or integer ndarray with same length as number of obs.
        
        total_counts: float
            Target total counts. If the count matrix has more than total_counts it will be downsampled to have this number.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/downsample.png

        Example
        -------
        import scanpy as sc

        sc.pp.downsample_counts(adata, counts_per_cell=1, total_counts=None, random_state=42)
        # counts now equal the total number of observations(cell)
        """
        st.subheader("Downsample data")
        counts_per_cell, total_counts = st.tabs(["counts_per_cell", "total_counts"])
        with counts_per_cell:
            with st.form(key="downsample_form_counts_per_cell"):
                counts_per_cell = st.number_input(label="Counts per cell", key="ni_downsample_counts_per_cell", help="Target total counts per cell. If a cell has more than 'counts_per_cell', it will be downsampled to this number. Resulting counts can be specified on a per cell basis by passing an array.Should be an integer or integer ndarray with same length as number of obs.")
                subcol1, _, _ = st.columns(3)
                btn_downsample_counts_per_cell = subcol1.form_submit_button(label="Apply", use_container_width=True)
                if btn_downsample_counts_per_cell:
                    sc.pp.downsample_counts(self.adata, counts_per_cell=counts_per_cell, random_state=42)
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
                    sc.pp.downsample_counts(self.adata, total_counts=total_counts, random_state=42)
                    #write to script state
                    st.session_state["script_state"].add_script(f"#Downsample dataset\nsc.pp.downsample_counts(adata, total_counts={total_counts})")
                    self.save_adata()
                    st.toast("Successfully downsampled data by total counts", icon="✅")

            
    def subsample_data(self):
        """
        Subsample to a fraction of the number of observations.

        Parameters
        ----------
        n_obs: int
            Subsample to this number of observations.
        
        fraction: float
            Subsample to this fraction of the number of observations.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/subsample.png

        Example
        -------
        import scanpy as sc

        # subsample data to 90% of original size (n_obs has no effect here)
        sc.pp.subsample(adata, n_obs=None, fraction=0.9, random_state=42)
        # subsample data to 1000 observations (fraction has no effect here)
        sc.pp.subsample(adata, n_obs=1000, fraction=None, random_state=42)
        """
        st.subheader("Subsample data")
        n_obs, fraction = st.tabs(["n_obs", "fraction"])
        with n_obs:
            with st.form(key="subsample_form_n_obs"):
                n_obs = st.number_input(label="n obs", key="ni_subsample_n_obs", help="Subsample to this number of observations.")
                subcol1, _, _ = st.columns(3)
                btn_subsample_n_obs = subcol1.form_submit_button(label="Apply", use_container_width=True)
                if btn_subsample_n_obs:
                    sc.pp.subsample(self.adata, n_obs=st.session_state.ni_n_obs, random_state=42)
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
                    sc.pp.subsample(self.adata, fraction=st.session_state.ni_subsample_fraction, random_state=42)
                    #write to script session
                    st.session_state["script_state"].add_script(f"#Subsample dataset\nsc.pp.subsample(adata, n_obs={st.session_state.ni_n_obs}, fraction={st.session_state.ni_subsample_fraction})")
                    self.save_adata()
                    st.toast("Successfully subsampled data", icon="✅")
                    
    def batch_effect_removal(self):
        """
        Uses ComBat for batch effect correction [Johnson07]_ [Leek12]_ [Pedersen12]_. Corrects for batch effects by fitting linear models, gains statistical power via an EB framework where information is borrowed across genes. This uses the implementation combat.py_ [Pedersen12]_.

        Parameters
        ----------
        key: str
            The batch key.

        covariates: List[str]
            Additional covariates besides the batch variable such as adjustment variables or biological condition. 
            This parameter refers to the design matrix X in Equation 2.1 in [Johnson07]_ and to the mod argument in the original combat function in the sva R package. 
            Note that not including covariates may introduce bias or lead to the removal of biological signal in unbalanced designs.
        
        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/batch_effect_removal.png

        Example
        -------
        import scanpy as sc

        sc.pp.combat(adata, key="batch", covariates=None)
        """
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
        """
        Computes PCA coordinates, loadings and variance decomposition. Uses the implementation of *scikit-learn* [Pedregosa11]_. This may be useful in the preprocessing stage, for example looking at batch effect or effect of doublets when forming clusters.

        Parameters
        ----------
        color: str
            PCA color. If the observation is categorical, this will group clusters into discreet colours. If observation is continuous, this will colour data on a scale represented by the colour bar.
        
        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/pca.png

        Example
        -------
        import scanpy as sc

        sc.pp.pca(adata, random_state=42)
        sc.pl.pca(adata, color="batch")
        """
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
                
    
    def measure_gene_counts(self):
        """
        Measure counts of genes within the whole dataset or compare across cell populations by supplying an obs key. 

        Parameters
        ----------
        gene: str
            The gene name to measure.

        obs_key: str
            Name of observation. This will measure gene counts within a sub-population.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/measure_gene_counts.png

        Example
        -------
        import scanpy as sc
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd

        fig, ax = plt.subplots()

        # subset of counts in observation
        fig, ax = plt.subplots()

        gene = 'ENSMUSG00000092341'
        obs_key = 'BATCH'
        df = pd.DataFrame({f'{gene} count': adata.to_df()[gene], f"{obs_key}": adata.obs[f"{obs_key}"]})
        ax.bar(df[f'{obs_key}'], df[f'{gene} count'])
        plt.title(f'Gene counts across {obs_key}')
        plt.xlabel(f'{obs_key}')
        plt.ylabel(f'{gene} count')
        plt.show()
        """
        
        st.subheader("Measure gene counts")
        single_dataset, subsample = st.tabs(['Whole dataset', 'Subsample'])
        
        with single_dataset:
            with st.form(key="measure_gene_counts_single_dataset"):
                st.subheader("Collective counts across dataset")
                options = self.adata.var_names
                genes = st.multiselect(label="Gene (e.g. XIST for detecting sex)", options=options)
                subcol_btn1, _, _ = st.columns(3)
                submit_btn = subcol_btn1.form_submit_button(label="Run", use_container_width=True)
                if submit_btn:
                    with st.spinner(text="Locating genes"):
                        df_whole_ds = pd.DataFrame({'genes': genes, 'counts': [self.adata.to_df()[gene].sum() for gene in genes]})
                        st.bar_chart(df_whole_ds, x='genes', y='counts', color='genes')
                        #write to script state
                        st.session_state["script_state"].add_script("#Measure gene counts in single dataset")
                        # st.session_state["script_state"].add_script(f"df_whole_ds = pd.DataFrame({{'genes': {genes}, 'counts': {[self.adata.to_df()[gene].sum() for gene in genes]}}})")
                        # st.session_state["script_state"].add_script(f"arr = np.array(['{st.session_state.adata_state.current.adata_name}'])")
                        # st.session_state["script_state"].add_script(f"df = pd.DataFrame('{gene} count': adata.obs['gene-counts'], 'Dataset': np.repeat(arr, adata.n_obs))")
        with subsample:
            with st.form(key="measure_gene_counts_multiple_datasets"):
                st.subheader("Subsample counts in dataset")
                gene_options = self.adata.var_names
                batch_key_measure_gene_counts = st.selectbox(label="Obs key", options=self.adata.obs_keys(), key="sb_sex_pred_batch_key")
                gene = st.selectbox(label="Gene (e.g. XIST for detecting sex)", options=gene_options)
                subcol_btn1, _, _ = st.columns(3)
                submit_btn = subcol_btn1.form_submit_button(label="Run", use_container_width=True)
                if submit_btn:
                    with st.spinner(text="Locating genes"):
                        df_subsample = pd.DataFrame({f'{gene} count': self.adata.to_df()[gene], f"{batch_key_measure_gene_counts}": self.adata.obs[f"{batch_key_measure_gene_counts}"]})
                        st.bar_chart(data=df_subsample, x=f"{batch_key_measure_gene_counts}", y=f'{gene} count', color=f"{batch_key_measure_gene_counts}")
                        #write to script state
                        st.session_state["script_state"].add_script("#Measure gene counts across datasets")
                        st.session_state["script_state"].add_script(f"gene = {gene}")
                        st.session_state["script_state"].add_script(f" batch_key_measure_gene_counts = {batch_key_measure_gene_counts}")
                        st.session_state["script_state"].add_script(f"df_subsample = pd.DataFrame({{f'{{gene}} count': self.adata.to_df()[gene], f'{{batch_key_measure_gene_counts}}': self.adata.obs[f'{{batch_key_measure_gene_counts}}']}})")
                  
    def cell_cycle_scoring(self):
        """
        Use marker genes to assign an S and G2/M score to approximate cell cycle distribution in dataset.

        Parameters
        ----------
        cell_cycle_file: UploadedFile
            A csv or tsv file containing marker genes for S and G2/M phase.

        gene_column: str
            The dataframe column name for the marker genes.

        phase_column: str
            The dataframe column name for the phase.

        group_by: Optional[str]
            The name of the obs value to group by when plotting.

        bandwidth: float
            The bandwidth of the violins (a higher value increases the width of each violin).

        jitter: float
            The offset of points on the violin plot.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/cell_cycle_scoring.png

        Example
        -------
        import scanpy as sc
        import pandas as pd
        import plotly.graph_objects as go
        import numpy as np
        import re

        df = pd.read_csv("./mus_musculus_cell_cycle.csv", delimiter=',', index_col=0)

        phase_column_index = 0
        gene_column_index = 1
                
        s_genes = df.iloc[:, phase_column_index].str.contains("s", flags=re.IGNORECASE, regex=True)
        g2m_genes = df.iloc[:, phase_column_index].str.contains("g2m", flags=re.IGNORECASE, regex=True)

        s_genes = df[s_genes].iloc[:, gene_column_index].values
        g2m_genes = df[g2m_genes].iloc[:, gene_column_index].values

        # Apply processing
        adata.raw = adata
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
        sc.pp.log1p(adata)
        sc.pp.scale(adata)

        # In this case group by batches
        group_by = 'BATCH'

        sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

        # Plot violin plot
        fig = go.Figure()

        s_score_df = pd.DataFrame({'phase': np.repeat('S_score', len(adata.obs['S_score'])), 'score': adata.obs['S_score']})
        g2m_score_df = pd.DataFrame({'phase': np.repeat('G2M_score', len(adata.obs['G2M_score'])), 'score': adata.obs['G2M_score']})

        violin_df = pd.concat([s_score_df, g2m_score_df])

        violin_df["group"] = adata.obs[group_by]

        fig.add_trace(go.Violin(x=violin_df['group'][violin_df['phase'] == 'S_score'], 
            y=violin_df['score'][violin_df['phase'] == 'S_score'],
            legendgroup='S', scalegroup='S', name='S',
            bandwidth=0.4, jitter=0.1, line_color='blue')
        )

        fig.add_trace(go.Violin(x=violin_df['group'][violin_df['phase'] == 'G2M_score'], 
            y=violin_df['score'][violin_df['phase'] == 'G2M_score'],
            legendgroup='G2M', scalegroup='G2M', name='G2M',
            bandwidth=0.4, jitter=0.1, line_color='orange')
        )

        fig.update_traces(meanline_visible=True)
        fig.update_layout(violingap=0, violinmode='group', xaxis_title=group_by, yaxis_title="Score", legend_title="Phase") #add legend title
        fig.show()
        """
        st.subheader("Cell cycle score")
        col1, col2 = st.columns(2)
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
        preprocess.filter_highly_variable_genes()
        preprocess.normalize_counts()
        preprocess.scale_to_unit_variance()
        preprocess.subsample_data()
        preprocess.downsample_data()
        

    with col2:
        preprocess.filter_cells()
        preprocess.filter_genes()
        preprocess.predict_doublets()
        preprocess.recipes()
        preprocess.batch_effect_removal()
        preprocess.measure_gene_counts()
        
            
    with col3:
        preprocess.remove_genes()
        preprocess.annotate()
        preprocess.pca()
        preprocess.regress_out()
        
        
        
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
