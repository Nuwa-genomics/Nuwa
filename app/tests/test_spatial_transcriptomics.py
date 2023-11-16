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

class Test_Spatial_Transcriptomics:
    def __init__(self, session_state = None):
        print(f"{bcolors.OKBLUE}Initialising page... {bcolors.ENDC}", end="")
        self.at = AppTest.from_file("pages/7_Spatial_Transcriptomics.py")
        if session_state is not None:
            self.at.session_state = session_state
        self.at.run(timeout=1000)
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        

    def get_final_session_state(self):
        return self.at.session_state

        