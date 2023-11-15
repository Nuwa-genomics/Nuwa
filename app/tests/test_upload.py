from streamlit.testing.v1 import AppTest
import streamlit as st
import scanpy as sc
import os
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

class Test_Upload:
    def __init__(self, session_state = None):
        print(f"{bcolors.OKBLUE}Initialising page...{bcolors.ENDC}")
        self.at = AppTest.from_file("pages/1_Upload.py")
        if session_state is not None:
            self.at.session_state = session_state

        self.adata = sc.datasets.pbmc3k()
        self.conn: Session = SessionLocal()
        
        self.at.run(timeout=100)
        self.test_load_adata()
        self.test_upload_data()
        self.test_insert_record()
        assert not self.at.exception

    def test_insert_record(self):
        print(f"{bcolors.OKBLUE}test_insert_record...{bcolors.ENDC}")
        #data should already by inserted when adata state is instantiated. check it is there
        query = self.conn.query(schemas.Adata).filter(schemas.Adata.adata_name == "adata_raw").filter(schemas.Adata.work_id == self.at.session_state.current_workspace.id)
        assert query.count() == 1

    def test_load_adata(self):
        print(f"{bcolors.OKBLUE}test_load_adata...{bcolors.ENDC}")
        #adata needs to be written first
        sc.write(filename=f"{os.getenv('WORKDIR')}uploads/adata_raw.h5ad", adata=self.adata)
        sc.write(filename=f"{os.getenv('WORKDIR')}adata/adata_raw.h5ad", adata=self.adata)
        self.at.session_state["adata_state"] = AdataState(workspace_id=self.at.session_state.current_workspace.id)
        
    def test_upload_data(self):
        print(f"{bcolors.OKBLUE}test_upload_data...{bcolors.ENDC}")
        path = os.path.join(os.getenv('WORKDIR'), 'adata/', "adata_raw.h5ad")
        file_read = sc.read_h5ad(path)
        assert file_read.obs.to_dict() == self.adata.obs.to_dict()
        assert file_read.var.to_dict() == self.adata.var.to_dict()

    def get_final_session_state(self):
        return self.at.session_state




        