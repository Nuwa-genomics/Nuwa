from streamlit.testing.v1 import AppTest
import streamlit as st
import time

from utils.AdataState import AdataState

class Test_Preprocess:
    def __init__(self, session_state = None):
        self.at = AppTest.from_file("pages/2_Preprocess.py")
        if session_state is not None:
            self.at.session_state = session_state
            
        self.at.run(timeout=100)
        
        assert not self.at.exception

    def get_final_session_state(self):
        return self.at.session_state





        