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
        self.at = AppTest.from_file("pages/4_Create_model.py")
        if session_state is not None:
            self.at.session_state = session_state

        self.at.run(timeout=100)
        
        if model != None:
            self.at.selectbox(key="sb_model_selection").select(model).run(timeout=100)
            #reduce epochs for quick testing
            if model == "Citeseq (dimensionality reduction)":
                self.at.number_input(key="ni_citeseq_epochs").set_value(10).run(timeout=100)
                self.at.number_input(key="ni_citeseq_lr").set_value(0.0015).run(timeout=100)
                self.at.selectbox(key="sb_citeseq_optim").set_value("SGD").run(timeout=100)
                self.at.selectbox(key="sb_citeseq_optim").set_value("Adam").run(timeout=100)
                self.at.slider(key="citeseq_train_test_split").set_value(89).run(timeout=100)
                
                assert self.at.session_state.model_obj["lr"] == 0.0015
                assert self.at.session_state.model_obj["n_epochs"] == 10
                assert self.at.session_state.model_obj["test_split"] == 0.11
                assert self.at.session_state.model_obj["optim"] == "Adam"
                assert self.at.session_state.model_obj["n_features"] == self.at.session_state.adata_state.current.adata.to_df().shape[1]
                
            elif model == "Solo (doublet removal)":
                self.at.number_input(key="ni_vae_epochs").set_value(10).run(timeout=100)
                self.at.number_input(key="ni_solo_lr").set_value(0.0015).run(timeout=100)
                self.at.slider(key="input_train_size_solo_vae").set_value(89).run(timeout=100)
                
                assert self.at.session_state.model_obj["lr"] == 0.0015
                assert self.at.session_state.model_obj["n_epochs"] == 10
                assert self.at.session_state.model_obj["train_size"] == 0.11

            else:
                print("Error: Unknown model")
                assert self.at.exception
            
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        

    def get_final_session_state(self):
        return self.at.session_state

        