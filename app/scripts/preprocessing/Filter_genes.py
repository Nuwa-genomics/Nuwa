from enums.Language import Language
import streamlit as st
from state.ScriptState import ScriptState
from scripts.Script import Script

class Filter_genes(Script):
    """
    Exports an R or python script for filtering genes.
    """

    def __init__(self, language: Language | str, min_cells: int, object: str = None):
        super().__init__(language=language)
        
        self.min_cells = min_cells
        self.object = object


    def add_script(self):

        script_state: ScriptState = st.session_state.script_state

        if self.language == Language.R or self.language == Language.R.value or self.language == Language.ALL_SUPPORTED:
            if self.object == None:
                self.object = "pbmc.data"
            script = f""" \
            \n# Load/reload Seurat object with min genes/cells parameters \
            \npbmc <- CreateSeuratObject(counts = {self.object}, project = "pbmc3k", min.cells = {self.min_cells})
            """
            script_state.add_script(script, language=Language.R)

        if self.language == Language.python or self.language == Language.python.value or self.language == Language.ALL_SUPPORTED:
            if self.object == None:
                self.object = "adata"
            script = f"""
                \n# Filter cells \
                \nsc.pp.filter_genes(self.{self.object}, min_cells={self.min_cells})
            """
            script_state.add_script(script, language=Language.python)

