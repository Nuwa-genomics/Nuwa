import streamlit as st
from state.ScriptState import ScriptState
from models.ScriptModel import Language

class Script:
    """
    Base class for scripts added to script state.
    """

    def __init__(self, language: Language):
        self.script_state: ScriptState = st.session_state.script_state

        if isinstance(language, Language):
            self.language = language
        else:
            print("Error: Unknown language")
            return

    def add_script(self):
        return