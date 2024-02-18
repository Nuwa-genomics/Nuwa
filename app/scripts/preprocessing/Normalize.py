from models.ScriptModel import Language
import streamlit as st
from state.ScriptState import ScriptState
import numpy as np

class Normalize:
    """
    Exports a python or R script for normalizing counts.
    """

    @staticmethod
    def add_script(language: Language | str, scale_factor: float, object: str = None, log_norm = True):

        script_state: ScriptState = st.session_state.script_state

        if language == Language.R or language == Language.R.value or language == Language.ALL_SUPPORTED:
            if object == None:
                object = "pbmc"
            script = f""" \
            \n# Normalize counts \
            \n{object} <- NormalizeData({object}, normalization.method = 'LogNormalize', scale.factor = {scale_factor}) \
            """
            script_state.add_script(script, language=Language.R)

        if language == Language.python or language == Language.python.value or language == Language.ALL_SUPPORTED:
            if object == None:
                object = "adata"
            script = f""" \
            \n# Normalize counts \
            \n log_norm = {log_norm} \
            \nsc.pp.normalize_total({object}, target_sum={scale_factor}) \
            \nif log_norm: \
            \n\tsc.pp.log1p({object}) \
            """
            script_state.add_script(script, language=Language.python)

        if not isinstance(language, Language):
            print("Error: Unknown language, not adding to script state")
            return
