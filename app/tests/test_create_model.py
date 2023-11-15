from streamlit.testing.v1 import AppTest
from utils.AdataState import AdataState
import time

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

class Test_Create_Model:
    def __init__(self, session_state = None, model=None):
        print(f"{bcolors.OKBLUE}Initialising page...{bcolors.ENDC}")
        self.at = AppTest.from_file("pages/3_Create_model.py")
        if session_state is not None:
            self.at.session_state = session_state

        self.at.run(timeout=100)
        
        if model != None:
            self.at.selectbox(key="sb_model_selection").select(model).run()
        assert not self.at.exception
        

    def get_final_session_state(self):
        return self.at.session_state

        