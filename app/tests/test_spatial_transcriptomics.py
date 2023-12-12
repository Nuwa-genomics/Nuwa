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
        self.at = AppTest.from_file("pages/8_Spatial_Transcriptomics.py")
        if session_state is not None:
            self.at.session_state = session_state
            
        self.at.run(timeout=1000)
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_spatial_scatter()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_neighbourhood_enrichment()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_ripley_score()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_cooccurance_score()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_interaction_matrix()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_centrality_score()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        self.test_ligand_receptor_interation()
        assert not self.at.exception
        print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
    def test_spatial_scatter(self):
        print(f"{bcolors.OKBLUE}test_spatial_scatter {bcolors.ENDC}", end="")
        
    def test_neighbourhood_enrichment(self):
        print(f"{bcolors.OKBLUE}test_neighbourhood_enrichment {bcolors.ENDC}", end="")

    def test_ripley_score(self):
        print(f"{bcolors.OKBLUE}test_ripley_score {bcolors.ENDC}", end="")
    
    def test_cooccurance_score(self):
        print(f"{bcolors.OKBLUE}test_cooccurance_score {bcolors.ENDC}", end="")
    
    def test_interaction_matrix(self):
        print(f"{bcolors.OKBLUE}test_interaction_matrix {bcolors.ENDC}", end="")
    
    def test_centrality_score(self):
        print(f"{bcolors.OKBLUE}test_centrality_score {bcolors.ENDC}", end="")
        
    def test_ligand_receptor_interation(self):
        print(f"{bcolors.OKBLUE}test_ligand_receptor_interaction {bcolors.ENDC}", end="")
    
    def get_final_session_state(self):
        return self.at.session_state

        