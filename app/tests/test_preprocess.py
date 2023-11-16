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
        
        
    def test_add_and_delete_adata(self):
        print(f"{bcolors.OKBLUE}test_add_and_delete_data {bcolors.ENDC}", end="")
        self.at.button(key="btn_add_adata").click().run()
        self.at.text_input(key="ti_new_adata_name").input("adata_test")
        self.at.text_input(key="ti_new_adata_notes").input("test_note")
        self.at.button(key="FormSubmitter:new_adata_form-Save").click().run()
        #test added to db
        assert self.conn.query(schemas.Adata).filter(schemas.Adata.adata_name == "adata_test") \
        .filter(schemas.Adata.notes == "test_note") \
        .filter(schemas.Adata.work_id == self.at.session_state.current_workspace.id).count() == 1
        #test delete adata
        self.at.selectbox(key="sb_adata_selection").select("adata_test").run()
        self.at.button(key="btn_delete_adata").click().run()
        assert self.conn.query(schemas.Adata).filter(schemas.Adata.adata_name == "adata_test") \
        .filter(schemas.Adata.notes == "test_note") \
        .filter(schemas.Adata.work_id == self.at.session_state.current_workspace.id).count() == 0
        #select original raw
        self.at.selectbox(key="sb_adata_selection").select("adata_raw").run()
        
    def test_filter_highest_expressed(self):
        print(f"{bcolors.OKBLUE}test_filter_highest_expressed {bcolors.ENDC}", end="")
        self.at.number_input(key="n_top_genes").increment() #increase to 21
        self.at.button(key="FormSubmitter:form_highest_expr-Filter").click().run()
        #TODO: Find a way to test if top genes are correct
        
    def test_highest_variable_genes(self):
        print(f"{bcolors.OKBLUE}test_highest_variable_genes {bcolors.ENDC}", end="")
        self.at.number_input(key="input_highly_variable_min_mean").increment()
        self.at.number_input(key="input_highly_variable_min_mean").decrement()
        self.at.number_input(key="input_highly_variable_max_mean").increment()
        self.at.number_input(key="input_highly_variable_max_mean").decrement()
        self.at.button(key="FormSubmitter:form_highly_variable-Filter").click().run()
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
        self.at.button(key="FormSubmitter:scrublet_form-Filter").click().run()
        

    def get_final_session_state(self):
        return self.at.session_state





        