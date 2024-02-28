import streamlit as st
from state.ScriptState import ScriptState

class Script:

    def __init__(self):
        self.script_state: ScriptState = st.session_state.script_state

    def add_script(self):
        return