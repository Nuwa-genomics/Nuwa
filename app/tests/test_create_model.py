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
        print(f"{bcolors.OKBLUE}Initialising page... {bcolors.ENDC}", end="")
        self.at = AppTest.from_file("pages/3_Create_model.py")
        if session_state is not None:
            self.at.session_state = session_state

        self.at.run(timeout=100)
        
        if model != None:
            self.at.selectbox(key="sb_model_selection").select(model).run(timeout=100)
            #reduce epochs for quick testing
            if model == "Citeseq (dimensionality reduction)":
                self.at.number_input(key="ni_citeseq_epochs").set_value(10).run(timeout=100)
            elif model == "Solo (doublet removal)":
                self.at.number_input(key="ni_vae_epochs").set_value(10).run(timeout=100)
            elif model == "DeepST (identify spatial domains)":
                self.at.number_input(key="ni_deepst_epochs").set_value(10).run(timeout=100)
                self.at.number_input(key="ni_deepst_preepochs").set_value(10).run(timeout=100)
            else:
                print("Error: Unknown model")
                assert self.at.exception
            
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        

    def get_final_session_state(self):
        return self.at.session_state

        