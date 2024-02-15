from models.ScriptModel import Language
import streamlit as st
from state.ScriptState import ScriptState

class Filter_cells:
    """
    Exports an R or python script for filtering cell counts.
    """

    @staticmethod
    def add_script(language: Language | str, min_genes: int, object: str = None):

        script_state: ScriptState = st.session_state.script_state

        if language == Language.R or language == Language.R.value or language == Language.ALL_SUPPORTED:
            if object == None:
                object = "pbmc.data"
            script = f""" \
            \n# Load/reload Seurat object with min genes/cells parameters \
            \npbmc <- CreateSeuratObject(counts = {object}, project = "pbmc3k", min.features = {min_genes})
            """
            script_state.add_script(script, language=Language.R)

        if language == Language.python or language == Language.python.value or language == Language.ALL_SUPPORTED:
            if object == None:
                object = "adata"
            script = f"""
                \n# Filter cells \
                \nsc.pp.filter_cells(self.{object}, min_genes={min_genes})
            """
            script_state.add_script(script, language=Language.python)

        if not isinstance(language, Language):
            print("Error: Unknown language, not adding to script state")
            return

        
