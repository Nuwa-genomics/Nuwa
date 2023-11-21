from streamlit.testing.v1 import AppTest
from utils.AdataState import AdataState
import time

from models.WorkspaceModel import WorkspaceModel
from database.database import SessionLocal
from sqlalchemy.orm import Session
from database.schemas import schemas

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

class Test_Dashboard:
    def __init__(self, workspace_name, session_state = None):
        print(f"{bcolors.OKBLUE}Initialising page... {bcolors.ENDC}")
        self.workspace_name = workspace_name
        self.conn: Session = SessionLocal()
        
        self.at = AppTest.from_file("Dashboard.py")
        if session_state is not None:
            self.at.session_state = session_state
            
        self.at.run(timeout=100)
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_create_new_workspace()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_delete_workspace()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_select_created_workspace()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        

    def test_create_new_workspace(self):
        print(f"{bcolors.OKBLUE}test_create_new_workspace... {bcolors.ENDC}")
        self.at.button("btn_new_workspace").click().run(timeout=100)
        self.at.text_input(key="ti_new_workspace_name").input(self.workspace_name)
        self.at.text_input(key="ti_new_workspace_desc").input("test description")
        self.at.button(key="FormSubmitter:new_workspace_form-Save").click().run(timeout=100)
        assert self.conn.query(schemas.Workspaces).filter(schemas.Workspaces.workspace_name == self.workspace_name).count() == 1
        

    def test_select_created_workspace(self):
        print(f"{bcolors.OKBLUE}test_select_created_workspace... {bcolors.ENDC}")
        for button in self.at.button:
            if button.label == self.workspace_name:
                button.click().run(timeout=100)
                
    def test_delete_workspace(self):
        print(f"{bcolors.OKBLUE}test_delete_workspace... {bcolors.ENDC}")
        #create new workspace
        self.at.button("btn_new_workspace").click().run(timeout=100)
        self.at.text_input(key="ti_new_workspace_name").input("test_delete")
        self.at.text_input(key="ti_new_workspace_desc").input("test_delete")
        self.at.button(key="FormSubmitter:new_workspace_form-Save").click().run(timeout=100)
        assert self.conn.query(schemas.Workspaces).filter(schemas.Workspaces.workspace_name == "test_delete").count() == 1
        #select workspace
        for button in self.at.button:
            if button.label == "test_delete":
                button.click().run(timeout=100)
        #delete workspace
        self.at.button(key="btn_delete_workspace").click().run(timeout=100)
        assert self.conn.query(schemas.Workspaces).filter(schemas.Workspaces.workspace_name == "test_delete").count() == 0

    def get_final_session_state(self):
        return self.at.session_state

        