from models.ScriptModel import Language
import streamlit as st
from state.ScriptState import ScriptState
from scripts.Script import Script

class Filter_mito(Script):
    """
    Exports an R or python script for annotating mitochondrial genes.
    """

    def __init__(self, language: Language | str, mito_pct: int, object: str = None):
        super().__init__(language=language)
        
        self.mito_pct = mito_pct
        self.object = object

    def add_script(self):

        if self.language == Language.R or self.language == Language.R.value or self.language == Language.ALL_SUPPORTED:
            if self.object == None:
                self.object = "pbmc"
            script = f""" \
            \n# Annotate mitochondrial genes \
            \n{self.object}[["percent.mt"]] <- PercentageFeatureSet({self.object}, pattern = "^MT-") \
            \nVlnPlot({self.object}, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
            """
            self.script_state.add_script(script, language=Language.R)

        if self.language == Language.python or self.language == Language.python.value or self.language == Language.ALL_SUPPORTED:
            if self.object == None:
                self.object = "adata"
            script = f"""
                \n# Annotate mitochondrial genes \
                \n{self.object}.var['mt'] = {self.object}.var_names.str.startswith(('MT-', 'mt-')) \
                \nsc.pp.calculate_qc_metrics({self.object}, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True) \
                \nsc.pl.violin({self.object}, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True) \
                \nsc.pl.scatter({self.object}, x='total_counts', y='pct_counts_mt') \
                \nsc.pl.scatter({self.object}, x='total_counts', y='n_genes_by_counts') \
                \n{self.object} = {self.object}[{self.object}.obs.pct_counts_mt < {self.mito_pct}, :]
            """
            self.script_state.add_script(script, language=Language.python)

        
