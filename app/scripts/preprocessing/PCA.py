from enums.Language import Language
import streamlit as st
from state.ScriptState import ScriptState
import numpy as np
from scripts.Script import Script

class PCA(Script):
    """
    Exports a python or R script for computing pca.
    """

    def __init__(self, language: Language | str, color: str = None, object: str = None):
        super().__init__(language=language)
        
        self.color = color
        self.object = object


    def add_script(self):

        if self.language == Language.R or self.language == Language.R.value or self.language == Language.ALL_SUPPORTED:
            if self.object == None:
                self.object = "pbmc"
            script = f""" \
            \n#Compute pca \
            \npbmc <- RunPCA({self.object}, features = VariableFeatures(object = {self.object})) \
            \nVizDimLoadings({self.object}, dims = 1:2, reduction = "pca") \
            \nDimPlot({self.object}, reduction = "pca") + NoLegend() \
            """
            self.script_state.add_script(script, language=Language.R)

        if self.language == Language.python or self.language == Language.python.value or self.language == Language.ALL_SUPPORTED:
            if self.object == None:
                self.object = "adata"
            script = f""" \
            \n# Compute PCA \
            \nsc.pp.pca({self.object}, svd_solver="arpack", random_state=42)
            \nsc.pl.pca({self.object}, color={self.color})
            """
            self.script_state.add_script(script, language=Language.python)