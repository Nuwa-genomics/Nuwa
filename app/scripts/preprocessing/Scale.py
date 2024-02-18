from models.ScriptModel import Language
import streamlit as st
from state.ScriptState import ScriptState
import numpy as np

class Scale:
    """
    Exports a python or R script for scaling data.
    """

    @staticmethod
    def add_script(language: Language | str, max_value: float, zero_center: bool = True, object: str = None):

        script_state: ScriptState = st.session_state.script_state

        if language == Language.R or language == Language.R.value or language == Language.ALL_SUPPORTED:
            if object == None:
                object = "pbmc"
            script = f""" \
            \n# TODO: find equivalent R script \
            """
            script_state.add_script(script, language=Language.R)

        if language == Language.python or language == Language.python.value or language == Language.ALL_SUPPORTED:
            if object == None:
                object = "adata"
            script = f""" \
            \n# Scale data \
            \nsc.pp.scale({object}, zero_center={zero_center}, max_value={max_value})
            """
            script_state.add_script(script, language=Language.python)

        if not isinstance(language, Language):
            print("Error: Unknown language, not adding to script state")
            return
