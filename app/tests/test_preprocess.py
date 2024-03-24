from streamlit.testing.v1 import AppTest
import os
from anndata import AnnData
from database.database import SessionLocal
from sqlalchemy.orm import Session
from database.schemas import schemas
import scanpy as sc
import math
import numpy as np
import pandas as pd
from utils.plotting import highest_expr_genes_box_plot, plot_doubletdetection_threshold_heatmap
from pandas.testing import assert_frame_equal, assert_series_equal

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

class Test_Preprocess:
    def __init__(self, session_state = None, pipeline: int = 1):
        print(f"{bcolors.OKBLUE}Initialising page... {bcolors.ENDC}", end="")
        self.at = AppTest.from_file("pages/2_Preprocess.py")
        self.conn: Session = SessionLocal()
        if session_state is not None:
            self.at.session_state = session_state
            
        self.at.run(timeout=500)
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")

        if pipeline == 1:
            self.pipeline1()
        elif pipeline == 2:
            self.pipeline2()
        elif pipeline == 3:
            self.pipeline3()
        else:
            raise Exception("Unknown run")
        
        
    def pipeline1(self):
        self.test_add_and_delete_adata()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_notes()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_filter_highest_expressed()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")

        self.test_scrublet()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")

        self.test_scale_data()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_normalize_data()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_highly_variable_genes()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_filter_cells()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_filter_genes()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_batch_effect_removal_and_pca()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")

        self.test_cell_cycle_scoring()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_downsampling_data_total_counts()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")

        self.test_downsampling_data_counts_per_cell()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")

        self.test_subsampling_data_fraction()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_regress_out()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        # self.test_save_adata()
        # assert not self.at.exception
        # print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
    def pipeline2(self):

        self.test_geneID_conversion()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")

        self.test_remove_genes()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")

        self.test_annot_mito()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_annot_ribo()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_annot_hb()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")

        self.test_measure_gene_counts()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")

        self.test_doubletdetection()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")

        self.test_subsampling_data_n_obs()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")


    def pipeline3(self):
        self.test_pp_recipe()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")


        
    def test_add_and_delete_adata(self):
        print(f"{bcolors.OKBLUE}test_add_and_delete_data {bcolors.ENDC}", end="")
        self.at.text_input(key="ti_new_adata_name").input("adata_test")
        self.at.text_input(key="ti_new_adata_notes").input("test_note")
        self.at.button(key="btn_add_adata").click().run(timeout=100)
        #test added to db
        assert self.conn.query(schemas.Adata).filter(schemas.Adata.adata_name == "adata_test") \
        .filter(schemas.Adata.notes == "test_note") \
        .filter(schemas.Adata.work_id == self.at.session_state.current_workspace.id).count() == 1
        #switch adata
        self.at.selectbox(key="sb_adata_selection").select("adata_test").run(timeout=100)
        #test notes
        assert self.at.text_area(key="sidebar_notes").value == "test_note"
        #test delete adata
        self.at.button(key="btn_delete_adata").click().run(timeout=100)
        assert self.conn.query(schemas.Adata).filter(schemas.Adata.adata_name == "adata_test") \
        .filter(schemas.Adata.notes == "test_note") \
        .filter(schemas.Adata.work_id == self.at.session_state.current_workspace.id).count() == 0
        #Select original adata
        self.at.selectbox(key="sb_adata_selection").select("mouse_mammary_epithelial").run(timeout=100)
        
    def test_filter_highest_expressed(self):
        print(f"{bcolors.OKBLUE}test_filter_highest_expressed {bcolors.ENDC}", end="")
        # Check inputs are correct
        adata: AnnData = self.at.session_state.adata_state.current.adata
        adata_file = self.at.session_state.adata_state.current.filename
        adata_from_file = sc.read_h5ad(adata_file)
        assert adata.n_obs == adata_from_file.n_obs == 9288
        assert adata.n_vars == adata_from_file.n_vars == 1222
        fig1 = highest_expr_genes_box_plot(adata, n_top=21)
        # Run highest expr box plots
        self.at.number_input(key="ni:pp:highly_variable:n_top_genes").set_value(21) #increase to 21
        self.at.button(key="FormSubmitter:form_highest_expr-Filter").click().run(timeout=100)
        fig2 = self.at.session_state["highest_expr_box_plot"]
        assert fig1 == fig2

        
    def test_highly_variable_genes(self):
        print(f"{bcolors.OKBLUE}test_highly_variable_genes {bcolors.ENDC}", end="")
        # Check inputs are correct
        adata1: AnnData = self.at.session_state.adata_state.current.adata.copy()
        adata_file = self.at.session_state.adata_state.current.filename
        adata_from_file = sc.read_h5ad(adata_file)
        assert adata1.n_obs == adata_from_file.n_obs == 9288
        assert adata1.n_vars == adata_from_file.n_vars == 1222

        sc.pp.highly_variable_genes(adata1, n_top_genes=2000, min_mean=0.0125, max_mean=3.0, min_disp=0.5)

        # Run highly variable
        self.at.slider(key="sl:pp:highly_variable:mean").set_range(lower=0.0125, upper=3.0)
        self.at.number_input(key="ni:pp:highly_variable:seurat_n_top_genes").set_value(2000)
        self.at.button(key="FormSubmitter:form_highly_variable_seurat-Run").click().run(timeout=100)

        adata2: AnnData = self.at.session_state.adata_state.current.adata.copy()

        assert_series_equal(adata1.var.highly_variable, adata2.var.highly_variable)
        assert_series_equal(adata1.var.means, adata2.var.means)
        assert_series_equal(adata1.var.dispersions_norm, adata2.var.dispersions_norm)
        assert_series_equal(adata1.var.dispersions, adata2.var.dispersions)


    def test_filter_cells(self):
        print(f"{bcolors.OKBLUE}test_filter_cells {bcolors.ENDC}", end="")
        # Check inputs are correct
        adata: AnnData = self.at.session_state.adata_state.current.adata.copy()
        adata_file = self.at.session_state.adata_state.current.filename
        adata_from_file = sc.read_h5ad(adata_file)
        assert adata.n_obs == adata_from_file.n_obs == 9288
        assert adata.n_vars == adata_from_file.n_vars == 1222
        # Run filter cells
        self.at.number_input(key="ni:pp:filter_cells:min_genes").set_value(250)
        #TODO: add other inputs
        self.at.button(key="FormSubmitter:form_filter_cells-Apply").click().run(timeout=100)
        # test filter from state
        sc.pp.filter_cells(adata, min_genes=250)
        assert self.at.session_state.adata_state.current.adata.n_obs == adata.n_obs == 9283
        assert self.at.session_state.adata_state.current.adata.n_vars == adata.n_vars == 1222
        # test adata from file
        from_file_adata = sc.read(self.at.session_state.adata_state.current.filename)
        assert len(from_file_adata.var) == len(adata.var) == 1222
        assert len(from_file_adata.obs) == len(adata.obs) == 9283
        

    def test_filter_genes(self):
        print(f"{bcolors.OKBLUE}test_filter_genes {bcolors.ENDC}", end="")
        # Check inputs are correct
        adata: AnnData = self.at.session_state.adata_state.current.adata.copy()
        adata_file = self.at.session_state.adata_state.current.filename
        adata_from_file = sc.read_h5ad(adata_file)
        assert adata.n_obs == adata_from_file.n_obs == 9283
        assert adata.n_vars == adata_from_file.n_vars == 1222
        # Run filter genes
        self.at.number_input(key="ni:pp:filter_genes:min_cells").set_value(5000)
        self.at.button(key="FormSubmitter:form_filter_genes-Apply").click().run(timeout=100)
        # Test filter in state
        sc.pp.filter_genes(adata, min_cells=5000)
        assert self.at.session_state.adata_state.current.adata.n_obs == 9283
        assert self.at.session_state.adata_state.current.adata.n_vars == 1216
        # Test filter in file
        from_file_adata = sc.read(self.at.session_state.adata_state.current.filename)
        assert from_file_adata.n_vars == 1216
        assert from_file_adata.n_obs == 9283

        
    def test_pp_recipe(self):
        print(f"{bcolors.OKBLUE}test_pp_recipe {bcolors.ENDC}", end="")
        adata: AnnData = sc.datasets.pbmc3k()
        self.at.session_state.adata_state.current.adata = adata
        sc.write(adata=adata, filename=self.at.session_state.adata_state.current.filename)
        sc.pp.recipe_seurat(adata, log=True)
        self.at.checkbox(key="cb:pp:recipe:seurat:log").set_value(True)
        self.at.button(key="FormSubmitter:form_seurat-Apply").click().run(timeout=100)
        assert_frame_equal(self.at.session_state.adata_state.current.adata.to_df(), adata.to_df())
        assert_frame_equal(self.at.session_state.adata_state.current.adata.obs, adata.obs)
        assert_frame_equal(self.at.session_state.adata_state.current.adata.var, adata.var)


    def test_scrublet(self):
        print(f"{bcolors.OKBLUE}test_scrublet {bcolors.ENDC}", end="")
        # Check inputs
        assert self.at.session_state.adata_state.current.adata.n_obs == 9288
        assert self.at.session_state.adata_state.current.adata.n_vars == 1222
        from_file_adata = sc.read(self.at.session_state.adata_state.current.filename)
        assert from_file_adata.n_vars == 1222
        assert from_file_adata.n_obs == 9288
        adata = self.at.session_state.adata_state.current.adata.copy()
        self.at.number_input(key="ni:pp:scrublet:sim_doublet_ratio").set_value(2.10)
        self.at.number_input(key="ni:pp:scrublet:expected_doublet_rate").set_value(0.06)
        self.at.number_input(key="ni:pp:scrublet:stdev_doublet_rate").set_value(0.02)
        self.at.selectbox(key="sb:pp:scrublet:batch_key").select("BATCH")
        self.at.button(key="FormSubmitter:scrublet_form-Run").click().run(timeout=100)

        sc.external.pp.scrublet(adata, sim_doublet_ratio=2.10, expected_doublet_rate=0.06, stdev_doublet_rate=0.02, verbose=False, batch_key="BATCH", random_state=42)

        assert_series_equal(adata.obs.doublet_score, self.at.session_state.adata_state.current.adata.obs.doublet_score)
        assert_series_equal(adata.obs.predicted_doublet, self.at.session_state.adata_state.current.adata.obs.predicted_doublet)

    def test_doubletdetection(self):
        print(f"{bcolors.OKBLUE}test_doubletdetection {bcolors.ENDC}", end="")
        self.at.button(key="FormSubmitter:doubletdetection_form-Run").click().run(timeout=100)
        #TODO: implement tests for charts

        
    def test_annot_mito(self):
        print(f"{bcolors.OKBLUE}test_annot_mito {bcolors.ENDC}", end="")
        #self.at.number_input(key="ni:pp:pct_counts_mt").set_value(5)
        #self.at.button(key="FormSubmitter:form_annotate_mito-Apply").click().run(timeout=100)
        
    def test_annot_ribo(self):
        print(f"{bcolors.OKBLUE}test_annot_ribo {bcolors.ENDC}", end="")
        #self.at.number_input(key="ni_pct_counts_ribo").set_value(5)
        #self.at.button(key="FormSubmitter:form_annotate_ribo-Apply").click().run(timeout=100)
        
    def test_annot_hb(self):
        print(f"{bcolors.OKBLUE}test_annot_hb {bcolors.ENDC}", end="")
        
    def test_measure_gene_counts(self):
        print(f"{bcolors.OKBLUE}test_measure_gene_counts {bcolors.ENDC}", end="")
        # Single dataset
        self.at.multiselect(key="ms:pp:measure_genes:genes").set_value(['XIST'])
        self.at.button(key="FormSubmitter:measure_gene_counts_single_dataset-Run").click().run(timeout=100)
        # Multiple datasets
        self.at.selectbox(key="sb:pp:measure_genes:batch").select(['BATCH'])
        self.at.multiselect(key="ms:pp:measure_genes_batch:genes").set_value(['XIST'])
        self.at.button(key="FormSubmitter:measure_gene_counts_multiple_datasets-Run").click().run(timeout=100)
        
    def test_remove_genes(self):
        print(f"{bcolors.OKBLUE}test_remove_genes {bcolors.ENDC}", end="")
        self.at.multiselect(key="ms:pp:remove_genes:genes").set_value(['MALAT1'])
        self.at.button(key="FormSubmitter:remove_genes_form-Run").click().run(timeout=100)
        assert 'MALAT1' not in self.at.session_state.adata_state.current.adata.var_names
        
        
    def test_regress_out(self):
        print(f"{bcolors.OKBLUE}test_regress_out {bcolors.ENDC}", end="")
        
    def test_normalize_data(self):
        print(f"{bcolors.OKBLUE}test_normalize_data {bcolors.ENDC}", end="")
        # Check inputs
        adata: AnnData = self.at.session_state.adata_state.current.adata
        adata_file = self.at.session_state.adata_state.current.filename
        adata_from_file = sc.read_h5ad(adata_file)
        assert adata.n_obs == adata_from_file.n_obs == 9288
        assert adata.n_vars == adata_from_file.n_vars == 1222
        # Run Normalize
        self.at.number_input(key="ni:pp:normalize_counts:target_sum").set_value(1.00)
        self.at.button(key="FormSubmitter:form_normalize_total-Apply").click().run(timeout=100)
        # Assert counts from adata
        sum_counts = round(self.at.session_state.adata_state.current.adata.to_df().iloc[:, :].values.sum())
        num_obs = self.at.session_state.adata_state.current.adata.n_obs
        assert sum_counts == num_obs
        # Assert counts from file
        from_file_adata = sc.read(self.at.session_state.adata_state.current.filename)
        sum_counts = round(from_file_adata.to_df().iloc[:, :].values.sum())
        num_obs = int(len(from_file_adata.obs))
        assert sum_counts == num_obs

    def test_notes(self):
        print(f"{bcolors.OKBLUE}test_notes {bcolors.ENDC}", end="")
        self.at.selectbox(key="sb_adata_selection").select("mouse_mammary_epithelial").run(timeout=150)
        self.at.text_area(key="sidebar_notes").input("Important notes").run(timeout=150)
        adata = self.conn.query(schemas.Adata).filter(schemas.Adata.adata_name == "mouse_mammary_epithelial") \
        .filter(schemas.Adata.work_id == self.at.session_state.current_workspace.id).first()
        assert adata.notes == "Important notes"
        self.at.selectbox(key="sb_adata_selection").select("mouse_mammary_epithelial").run(timeout=150)
        assert self.at.text_area(key="sidebar_notes").value == "Important notes"
        
    def test_save_adata(self):
        print(f"{bcolors.OKBLUE}test_save_adata {bcolors.ENDC}", end="")
        self.at.selectbox(key="sb_adata_selection").select("mouse_mammary_epithelial").run(timeout=150)
        self.at.button(key="btn_save_adata").click().run(timeout=100)
        assert os.path.isfile(os.path.join(os.getenv('WORKDIR'), 'downloads', 'mouse_mammary_epithelial', 'mouse_mammary_epithelial.h5ad'))
        downloaded_adata = sc.read_h5ad(os.path.join(os.getenv('WORKDIR'), 'downloads', 'mouse_mammary_epithelial', 'mouse_mammary_epithelial.h5ad'))
        current_adata = sc.read_h5ad(os.path.join(os.getenv('WORKDIR'), 'adata', 'mouse_mammary_epithelial.h5ad'))
        #test seurat format
        self.at.checkbox(key="cb_seurat_format").check().run(timeout=100)
        self.at.button(key="btn_save_adata").click().run(timeout=500)
        assert os.path.isfile(os.path.join(os.getenv('WORKDIR'), 'downloads', 'mouse_mammary_epithelial', 'seurat', 'barcodes.tsv'))
        assert os.path.isfile(os.path.join(os.getenv('WORKDIR'), 'downloads', 'mouse_mammary_epithelial', 'seurat', 'features.tsv'))
        assert os.path.isfile(os.path.join(os.getenv('WORKDIR'), 'downloads', 'mouse_mammary_epithelial', 'seurat', 'matrix.mtx'))
        assert os.path.isfile(os.path.join(os.getenv('WORKDIR'), 'downloads', 'mouse_mammary_epithelial', 'seurat', 'metadata.csv'))
        
        
    def test_scale_data(self):
        print(f"{bcolors.OKBLUE}test_scale_adata {bcolors.ENDC}", end="")
        # Check inputs are correct
        adata: AnnData = self.at.session_state.adata_state.current.adata
        adata_file = self.at.session_state.adata_state.current.filename
        adata_from_file = sc.read_h5ad(adata_file)
        assert adata.n_obs == adata_from_file.n_obs == 9288
        assert adata.n_vars == adata_from_file.n_vars == 1222
        # Run scaling
        self.at.number_input(key="ni:pp:scale_data:max_value").set_value(10)
        self.at.button(key="FormSubmitter:scale_to_unit_variance_form-Apply").click().run(timeout=100)
        #assert each row is equal to scaled value, then check the mean is 0 and std 1 and within defined range
        assert round(self.at.session_state.adata_state.current.adata.to_df().mean().mean(), ndigits=1) == 0.0 # test mean
        assert round(self.at.session_state.adata_state.current.adata.to_df().std().mean(), ndigits=1) == 1.0 # test std
        assert round(self.at.session_state.adata_state.current.adata.to_df().max().max(), ndigits=1) <= 10.0 # test clipped max value
        

    def test_batch_effect_removal_and_pca(self):
        print(f"{bcolors.OKBLUE}test_batch_effect_removal_and_adata {bcolors.ENDC}", end="")

        adata: AnnData = self.at.session_state.adata_state.current.adata.copy()
        sc.pp.combat(adata, key='BATCH', inplace=True)
        sc.pp.pca(adata, random_state=42)

        # Run Combat
        self.at.selectbox(key="sb:pp:combat:batch_key").select("BATCH")
        self.at.button(key="FormSubmitter:batch_effect_removal_form-Apply").click().run(timeout=500)
        
        # Now run PCA
        self.at.selectbox(key="sb:pp:pca:color").select("BATCH")
        self.at.button(key="FormSubmitter:pca_pp_form-Apply").click().run(timeout=100)

        assert np.array_equal(adata.obsm['X_pca'], self.at.session_state.adata_state.current.adata.obsm['X_pca'])

        # Test plot generated
        pp_pca_df = pd.DataFrame({'pca1': adata.obsm['X_pca'][:,0], 'pca2': adata.obsm['X_pca'][:,1], 'color': adata.obs['BATCH']})  
        plots_dict: dict = self.at.session_state["preprocess_plots"]["pca"]
        assert_frame_equal(plots_dict.get('df'), pp_pca_df)

        
    def test_downsampling_data_total_counts(self):
        print(f"{bcolors.OKBLUE}test_downsampling_data_total_counts {bcolors.ENDC}", end="")
        # Total counts
        self.at.number_input(key="ni:pp:downsample:total_counts").set_value(1000.0)
        self.at.button(key="FormSubmitter:downsample_form_total_counts-Apply").click().run(timeout=100)
        assert self.at.session_state.adata_state.current.adata.to_df().sum().sum() == 1000.0


    def test_downsampling_data_counts_per_cell(self):
        print(f"{bcolors.OKBLUE}test_downsampling_data_counts_per_cell {bcolors.ENDC}", end="")
        # Counts per cell
        # Reset to raw to avoid negative values which skew counts
        # TODO: Possibly reset to raw via another method?
        self.at.number_input(key="ni:pp:downsample:counts_per_cell").set_value(2)
        self.at.button(key="FormSubmitter:downsample_form_counts_per_cell-Apply").click().run(timeout=100)
        #assert self.at.session_state.adata_state.current.adata.to_df().sum().sum() == self.at.session_state.adata_state.current.adata.n_obs * 2
        
        
    def test_subsampling_data_fraction(self):
        print(f"{bcolors.OKBLUE}test_subsampling_data_fraction {bcolors.ENDC}", end="")
        original_n_obs = self.at.session_state.adata_state.current.adata.n_obs
        fraction = 0.95
        self.at.number_input(key="ni:pp:subsample:fraction").set_value(fraction)
        self.at.button(key="FormSubmitter:subsample_form_fraction-Apply").click().run(timeout=100)
        subsampled_n_obs = self.at.session_state.adata_state.current.adata.n_obs
        assert subsampled_n_obs == math.floor(original_n_obs * fraction)


    def test_subsampling_data_n_obs(self):
        print(f"{bcolors.OKBLUE}test_subsampling_data_n_obs {bcolors.ENDC}", end="")
        fraction = 0.95
        original_n_obs = self.at.session_state.adata_state.current.adata.n_obs
        n_obs = math.floor(original_n_obs * fraction)
        print(n_obs)
        self.at.number_input(key="ni:pp:subsample:n_obs").set_value(n_obs)
        self.at.button(key="FormSubmitter:subsample_form_n_obs-Apply").click().run(timeout=100)
        print(self.at.session_state.adata_state.current.adata.n_obs)
        assert self.at.session_state.adata_state.current.adata.n_obs == n_obs

    def test_geneID_conversion(self):
        print(f"{bcolors.OKBLUE}test_geneID_conversion {bcolors.ENDC}", end="")
        self.at.toggle(key="toggle_gene_format").set_value(False).run()


    def test_cell_cycle_scoring(self):
        print(f"{bcolors.OKBLUE}test_cell_cycle_scoring {bcolors.ENDC}", end="")
        #TODO: Figure out how to run this test, can't simulate loading a file
        # self.at.selectbox(key="sb_gene_col_cell_cycle").select("genes").run()
        # self.at.selectbox(key="sb_phase_col_cell_cycle").select("phase").run()
        # self.at.selectbox(key="sb_group_cell_cycle").select("BATCH").run()
        # self.at.button(key="FormSubmitter:cell_cycle_scoring_form-Run").click().run()

    def get_final_session_state(self):
        return self.at.session_state
    
    
        