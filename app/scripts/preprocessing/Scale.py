from models.ScriptModel import Language
import streamlit as st
from state.ScriptState import ScriptState
import numpy as np
from scripts.Script import Script

class Scale(Script):
    """
    Exports a python or R script for scaling data.
    """

    def __init__(self, language: Language | str, max_value: float, zero_center: bool, object: str = None):
        super().__init__(language=language)
        
        self.max_value = max_value
        self.zero_center = zero_center
        self.object = object


    def add_script(self):

        if self.language == Language.R or self.language == Language.R.value or self.language == Language.ALL_SUPPORTED:
            if self.object == None:
                self.object = "pbmc"
            script = f""" \
            \n# TODO: find equivalent R script \
            """
            self.script_state.add_script(script, language=Language.R)

        if self.language == Language.python or self.language == Language.python.value or self.language == Language.ALL_SUPPORTED:
            if self.object == None:
                self.object = "adata"
            script = f""" \
            \n# Scale data \
            \nsc.pp.scale({self.object}, zero_center={self.zero_center}, max_value={self.max_value})
            """
            self.script_state.add_script(script, language=Language.python)
