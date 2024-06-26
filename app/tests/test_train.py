from streamlit.testing.v1 import AppTest
from state.AdataState import AdataState
import time
import streamlit as st

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

class Test_Train:
    def __init__(self, session_state = None):

        #TODO: Figure out why this fails on exit ??
        try:
            print(f"{bcolors.OKBLUE}Initialising page... {bcolors.ENDC}", end="")
            self.at = AppTest.from_file("pages/5_Train.py")
            if session_state is not None:
                self.at.session_state = session_state

            self.at.run(timeout=1000)
            
        except Exception as e:
            print("Error: ", e)
        finally:
            assert not self.at.exception
            print(f"{bcolors.OKGREEN}OK{bcolors.ENDC}")
        
        
    def get_final_session_state(self):
        return self.at.session_state

        