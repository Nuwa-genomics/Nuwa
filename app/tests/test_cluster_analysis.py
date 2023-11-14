from streamlit.testing.v1 import AppTest
from utils.AdataState import AdataState
import time

class Test_Cluster_Analysis:
    def __init__(self, session_state = None):
        self.at = AppTest.from_file("pages/5_Cluster_Analysis.py")
        if session_state is not None:
            self.at.session_state = session_state
        self.at.run(timeout=1000)
        assert not self.at.exception
        

    def get_final_session_state(self):
        return self.at.session_state

        