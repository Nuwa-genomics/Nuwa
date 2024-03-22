from enums.Language import Language
import streamlit as st
from state.ScriptState import ScriptState
from scripts.Script import Script

class Highest_expr_genes(Script):
    """
    Exports an R or python script for plotting highest expressed genes.
    """

    def __init__(self, language: Language | str, n_top_genes: int = 20, object: str = None):
        super().__init__(language=language)
        
        self.n_top_genes = n_top_genes
        self.object = object


    def add_script(self):

        if self.language == Language.R or self.language == Language.R.value or self.language == Language.ALL_SUPPORTED:
            if self.object == None:
                self.object = "pbmc.data"
            script = f""" \
            \n# This uses the scater library \
            \nplotHighestExprs( \
                \n\t{self.object}, \
                \n\tn = {self.n_top_genes}, \
                \n\tcolour_cells_by = NULL, \
                \n\tdrop_features = NULL, \
                \n\texprs_values = "counts", \
                \n\tby_exprs_values = exprs_values, \
                \n\tfeature_names_to_plot = NULL, \
                \n\tas_percentage = TRUE, \
                \n\tswap_rownames = NULL \
            \n)
            """
            self.script_state.add_script(script, language=Language.R)

        if self.language == Language.python or self.language == Language.python.value or self.language == Language.ALL_SUPPORTED:
            if self.object == None:
                self.object = "adata"
            script = f"""
                \n# Plot highest expr genes \
                \nsc.pl.highest_expr_genes({self.object}, n_top={self.n_top_genes})
            """
            self.script_state.add_script(script, language=Language.python)

