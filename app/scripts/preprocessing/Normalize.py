from enums.Language import Language
import streamlit as st
from state.ScriptState import ScriptState
import numpy as np
from scripts.Script import Script

class Normalize(Script):
    """
    Exports a python or R script for normalizing counts.
    """

    def __init__(self, language: Language | str, scale_factor: float, object: str = None, log_norm = True):
        super().__init__(language=language)
        
        self.scale_factor = scale_factor
        self.log_norm = log_norm
        self.object = object


    def add_script(self):

        if self.language == Language.R or self.language == Language.R.value or self.language == Language.ALL_SUPPORTED:
            if self.object == None:
                self.object = "pbmc"
            script = f""" \
            \n# Normalize counts \
            \n{self.object} <- NormalizeData({self.object}, normalization.method = 'LogNormalize', scale.factor = {self.scale_factor}) \
            """
            self.script_state.add_script(script, language=Language.R)

        if self.language == Language.python or self.language == Language.python.value or self.language == Language.ALL_SUPPORTED:
            if self.object == None:
                self.object = "adata"
            script = f""" \
            \n# Normalize counts \
            \n log_norm = {self.log_norm} \
            \nsc.pp.normalize_total({self.object}, target_sum={self.scale_factor}) \
            \nif log_norm: \
            \n\tsc.pp.log1p({self.object}) \
            """
            self.script_state.add_script(script, language=Language.python)

