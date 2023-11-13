from streamlit.testing.v1 import AppTest
from utils.AdataState import AdataState
import time

class Test_Spatial_Transcriptomics:
    def __init__(self, workspace_name, session_state = None):
        self.workspace_name = workspace_name
        self.at = AppTest.from_file("pages/7_Spatial_Transcriptomics.py")
        if session_state is not None:
            self.at.session_state = session_state
        self.at.run()
        assert not self.at.exception
        

    def get_final_session_state(self):
        return self.at.session_state

        