from models.ScriptModel import Language
import streamlit as st
from state.ScriptState import ScriptState
import numpy as np

class PCA:
    """
    Exports a python or R script for computing pca.
    """

    @staticmethod
    def add_script(language: Language | str, object: str = None, color: str = None):

        script_state: ScriptState = st.session_state.script_state

        if language == Language.R or language == Language.R.value or language == Language.ALL_SUPPORTED:
            if object == None:
                object = "pbmc"
            script = f""" \
            \n#Compute pca \
            \npbmc <- RunPCA({object}, features = VariableFeatures(object = {object})) \
            \nVizDimLoadings({object}, dims = 1:2, reduction = "pca") \
            \nDimPlot({object}, reduction = "pca") + NoLegend() \
            """
            script_state.add_script(script, language=Language.R)

        if language == Language.python or language == Language.python.value or language == Language.ALL_SUPPORTED:
            if object == None:
                object = "adata"
            script = f""" \
            \n# Compute PCA \
            \nsc.pp.pca({object}, svd_solver="arpack", random_state=42)
            \nsc.pl.pca({object}, color={color})
            """
            script_state.add_script(script, language=Language.python)

        if not isinstance(language, Language):
            print("Error: Unknown language, not adding to script state")
            return
