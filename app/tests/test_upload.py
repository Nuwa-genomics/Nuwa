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
        print(f"{bcolors.OKBLUE}Initialising page... {bcolors.ENDC}", end="")
        self.at = AppTest.from_file("pages/1_Upload.py")
        if session_state is not None:
            self.at.session_state = session_state

        self.conn: Session = SessionLocal()
        
        self.at.run(timeout=100)
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        print(f"{bcolors.OKBLUE}Testing select dataset... {bcolors.ENDC}", end="")
        self.at.selectbox(key="sb_sc_datasets").select("pbmc3k").run()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")

    def get_final_session_state(self):
        return self.at.session_state




        