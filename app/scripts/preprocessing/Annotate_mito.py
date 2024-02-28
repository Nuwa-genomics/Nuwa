from models.ScriptModel import Language
import streamlit as st
from state.ScriptState import ScriptState

class Annotate_mito:
    """
    Exports an R or python script for annotating mitochondrial genes.
    """

    def add_script(language: Language | str, object: str = None):

        script_state: ScriptState = st.session_state.script_state

        if language == Language.R or language == Language.R.value or language == Language.ALL_SUPPORTED:
            if object == None:
                object = "pbmc"
            script = f""" \
            \n# Annotate mitochondrial genes \
            \n{object}[["percent.mt"]] <- PercentageFeatureSet({object}, pattern = "^MT-") \
            \nVlnPlot({object}, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
            """
            script_state.add_script(script, language=Language.R)

        if language == Language.python or language == Language.python.value or language == Language.ALL_SUPPORTED:
            if object == None:
                object = "adata"
            script = f"""
                \n# Annotate mitochondrial genes \
                \n{object}.var['mt'] = {object}.var_names.str.startswith(('MT-', 'mt-')) \
                \nsc.pp.calculate_qc_metrics({object}, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True) \
                \nsc.pl.violin({object}, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True) \
                \nsc.pl.scatter({object}, x='total_counts', y='pct_counts_mt') \
                \nsc.pl.scatter({object}, x='total_counts', y='n_genes_by_counts')
            """
            script_state.add_script(script, language=Language.python)

        if not isinstance(language, Language):
            print("Error: Unknown language, not adding to script state")
            return
        

class Filter_mito:
    """
    Exports an R or python script for filtering out mitochondrial genes.
    """

    def add_script(language: Language | str, mito_pct: int, object: str = None):

        script_state: ScriptState = st.session_state.script_state

        if language == Language.R or language == Language.R.value or language == Language.ALL_SUPPORTED:
            if object == None:
                object = "pbmc"
            script = f""" \
            \n# TODO: Add R script for filtering out genes
            """
            script_state.add_script(script, language=Language.R)

        if language == Language.python or language == Language.python.value or language == Language.ALL_SUPPORTED:
            if object == None:
                object = "adata"
            script = f"""
                \n# Filter out counts with high mitochondrial gene count \
                \n{object} = {object}[{object}.obs.pct_counts_mt < {mito_pct}, :]
            """
            script_state.add_script(script, language=Language.python)

        if not isinstance(language, Language):
            print("Error: Unknown language, not adding to script state")
            return
        
