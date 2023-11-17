from streamlit.testing.v1 import AppTest
import streamlit as st
import time

from utils.AdataState import AdataState
from models.AdataModel import AdataModel
from utils.AdataState import AdataState
from models.AdataModel import AdataModel
from models.WorkspaceModel import WorkspaceModel
from database.database import SessionLocal
from sqlalchemy.orm import Session
from database.schemas import schemas
from utils.AdataState import AdataState

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
    def __init__(self, session_state = None):
        print(f"{bcolors.OKBLUE}Initialising page... {bcolors.ENDC}", end="")
        self.at = AppTest.from_file("pages/2_Preprocess.py")
        self.conn: Session = SessionLocal()
        if session_state is not None:
            self.at.session_state = session_state
            
        self.at.run(timeout=100)
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_add_and_delete_adata()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_notes()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_filter_highest_expressed()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_highest_variable_genes()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_filter_cells()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_filter_genes()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_pp_recipe()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_doublet_prediction()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_annot_mito()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_annot_ribo()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_normalize_data()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        
    def test_add_and_delete_adata(self):
        print(f"{bcolors.OKBLUE}test_add_and_delete_data {bcolors.ENDC}", end="")
        self.at.button(key="btn_add_adata").click().run(timeout=100)
        self.at.text_input(key="ti_new_adata_name").input("adata_test")
        self.at.text_input(key="ti_new_adata_notes").input("test_note")
        self.at.button(key="FormSubmitter:new_adata_form-Save").click().run(timeout=100)
        #test added to db
        assert self.conn.query(schemas.Adata).filter(schemas.Adata.adata_name == "adata_test") \
        .filter(schemas.Adata.notes == "test_note") \
        .filter(schemas.Adata.work_id == self.at.session_state.current_workspace.id).count() == 1
        #switch adata
        self.at.selectbox(key="sb_adata_selection").select("adata_test").run(timeout=100)
        #test notes
        assert self.at.text_area(key="sidebar_notes").value == "test_note"
        #test delete adata
        self.at.button(key="btn_delete_adata").click().run()
        assert self.conn.query(schemas.Adata).filter(schemas.Adata.adata_name == "adata_test") \
        .filter(schemas.Adata.notes == "test_note") \
        .filter(schemas.Adata.work_id == self.at.session_state.current_workspace.id).count() == 0
        #Select original raw
        self.at.selectbox(key="sb_adata_selection").select("adata_raw").run(timeout=100)
        
    def test_filter_highest_expressed(self):
        print(f"{bcolors.OKBLUE}test_filter_highest_expressed {bcolors.ENDC}", end="")
        self.at.number_input(key="n_top_genes").increment() #increase to 21
        self.at.button(key="FormSubmitter:form_highest_expr-Filter").click().run(timeout=100)
        #TODO: Find a way to test if top genes are correct
        
    def test_highest_variable_genes(self):
        print(f"{bcolors.OKBLUE}test_highest_variable_genes {bcolors.ENDC}", end="")
        self.at.number_input(key="input_highly_variable_min_mean").set_value(0.0125)
        self.at.number_input(key="input_highly_variable_max_mean").set_value(3.00)
        self.at.button(key="FormSubmitter:form_highly_variable-Filter").click().run(timeout=100)
        #TODO: add test to check adata
        
    def test_filter_cells(self):
        print(f"{bcolors.OKBLUE}test_filter_cells {bcolors.ENDC}", end="")
        self.at.number_input(key="filter_cell_min_genes").set_value(3)
        #TODO: add other inputs
        self.at.button(key="FormSubmitter:form_filter_cells-Apply").click().run()
        
    def test_filter_genes(self):
        print(f"{bcolors.OKBLUE}test_filter_genes {bcolors.ENDC}", end="")
        self.at.number_input(key="filter_gene_min_cells").set_value(100)
        self.at.button(key="FormSubmitter:form_filter_genes-Apply").click().run()
        
    def test_pp_recipe(self):
        print(f"{bcolors.OKBLUE}test_pp_recipe {bcolors.ENDC}", end="")
        self.at.selectbox(key="sb_pp_recipe").select("Seurat")
        self.at.button(key="FormSubmitter:form_recipes-Apply").click().run()
        
    def test_doublet_prediction(self):
        print(f"{bcolors.OKBLUE}test_doublet_prediction {bcolors.ENDC}", end="")
        self.at.number_input(key="ni_sim_doublet_ratio").set_value(2)
        self.at.number_input(key="ni_expected_doublet_rate").set_value(0.05)
        self.at.button(key="FormSubmitter:scrublet_form-Filter").click().run(timeout=100)
        
    def test_annot_mito(self):
        print(f"{bcolors.OKBLUE}test_annot_mito {bcolors.ENDC}", end="")
        self.at.number_input(key="ni_pct_counts_mt").set_value(5)
        self.at.button(key="FormSubmitter:form_annotate_mito-Apply").click().run()
        
    def test_annot_ribo(self):
        print(f"{bcolors.OKBLUE}test_annot_ribo {bcolors.ENDC}", end="")
        self.at.number_input(key="ni_pct_counts_ribo").set_value(5)
        self.at.button(key="FormSubmitter:form_annotate_ribo-Apply").click().run()
        
    def test_normalize_data(self):
        print(f"{bcolors.OKBLUE}test_normalize_data {bcolors.ENDC}", end="")
        self.at.number_input(key="ni_target_sum").set_value(1)
        self.at.button(key="FormSubmitter:form_normalize_total-Apply").click().run()
        
    def test_notes(self):
        print(f"{bcolors.OKBLUE}test_notes {bcolors.ENDC}", end="")
        self.at.selectbox(key="sb_adata_selection").select("adata_raw").run()
        self.at.text_area(key="sidebar_notes").input("Important notes").run()
        adata = self.conn.query(schemas.Adata).filter(schemas.Adata.adata_name == "adata_raw") \
        .filter(schemas.Adata.work_id == self.at.session_state.current_workspace.id).first()
        assert adata.notes == "Important notes"
        assert self.at.text_area(key="sidebar_notes").value == "Important notes"

    def get_final_session_state(self):
        return self.at.session_state





        