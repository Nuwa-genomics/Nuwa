from streamlit.testing.v1 import AppTest
import streamlit as st
import os

from utils.AdataState import AdataState
from models.AdataModel import AdataModel
from utils.AdataState import AdataState
from models.AdataModel import AdataModel
from models.WorkspaceModel import WorkspaceModel
from database.database import SessionLocal
from sqlalchemy.orm import Session
from database.schemas import schemas
from utils.AdataState import AdataState

import scanpy as sc


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

class Test_Integrate:
    def __init__(self, session_state = None):
        print(f"{bcolors.OKBLUE}Initialising page... {bcolors.ENDC}", end="")
        self.at = AppTest.from_file("pages/3_Integrate.py")
        self.conn: Session = SessionLocal()
        if session_state is not None:
            self.at.session_state = session_state
            
        self.at.run(timeout=1000)
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_ingest()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_concat_dataset()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_quick_map()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_scanorama_integrate()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_umap()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
    def test_ingest(self):
        print(f"{bcolors.OKBLUE}test_ingest{bcolors.ENDC}", end="")
        
    def test_concat_dataset(self):
        print(f"{bcolors.OKBLUE}test_concat_dataset{bcolors.ENDC}", end="")
        
    def test_quick_map(self):
        print(f"{bcolors.OKBLUE}test_quick_map{bcolors.ENDC}", end="")
        
    def test_scanorama_integrate(self):
        print(f"{bcolors.OKBLUE}test_scanorama_integrate{bcolors.ENDC}", end="")
        
    def test_umap(self):
        print(f"{bcolors.OKBLUE}test_umap{bcolors.ENDC}", end="")
        
    def get_final_session_state(self):
        return self.at.session_state