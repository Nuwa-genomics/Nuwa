from models.ScriptModel import Language
import streamlit as st
from state.ScriptState import ScriptState

class Highest_expr_genes:
    """
    Exports an R or python script for plotting highest expressed genes.
    """

    @staticmethod
    def add_script(language: Language | str, n_top_genes: int = 20, object: str = None):

        script_state: ScriptState = st.session_state.script_state

        if language == Language.R or language == Language.R.value or language == Language.ALL_SUPPORTED:
            if object == None:
                object = "pbmc.data"
            script = f""" \
            \n# This uses the scater library \
            \nplotHighestExprs( \
                \n\t{object}, \
                \n\tn = {n_top_genes}, \
                \n\tcolour_cells_by = NULL, \
                \n\tdrop_features = NULL, \
                \n\texprs_values = "counts", \
                \n\tby_exprs_values = exprs_values, \
                \n\tfeature_names_to_plot = NULL, \
                \n\tas_percentage = TRUE, \
                \n\tswap_rownames = NULL \
            \n)
            """
            script_state.add_script(script, language=Language.R)

        if language == Language.python or language == Language.python.value or language == Language.ALL_SUPPORTED:
            if object == None:
                object = "adata"
            script = f"""
                \n# Plot highest expr genes \
                \nsc.pl.highest_expr_genes({object}, n_top={n_top_genes})
            """
            script_state.add_script(script, language=Language.python)

        if not isinstance(language, Language):
            print("Error: Unknown language, not adding to script state")
            return

