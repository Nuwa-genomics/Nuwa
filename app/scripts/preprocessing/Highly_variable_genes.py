from models.ScriptModel import Language
import streamlit as st
from state.ScriptState import ScriptState
import numpy as np

class Highly_variable_genes:
    """
    Exports a python or R script for computing highly variable genes.
    """

    @staticmethod
    def add_script(language: Language | str, min_mean: float, max_mean: float, min_disp: float = 0.5, max_disp=np.Inf, n_top_genes: int = 2000, span=0.3, object: str = None):

        script_state: ScriptState = st.session_state.script_state

        if language == Language.R or language == Language.R.value or language == Language.ALL_SUPPORTED:
            if object == None:
                object = "pbmc.data"
            script = f""" \
            \n#Compute highly variable genes. This uses the FindVariableFeatures method in seurat \
            \n#If you have been using bioconductor's sce object you will need to reload it in as a seurat object. \
            \n{object} <- Read10X(data.dir = 'path_to_mtx_files') \
            \n{object} <- CreateSeuratObject(counts = {object}, project = 'pbmc3k') \
            \n#Log normalize \
            \n{object} <- NormalizeData({object}, normalization.method = 'LogNormalize', scale.factor = 10000) \
            \n{object} <- FindVariableFeatures({object}, selection.method = 'dispersion', nfeatures = {n_top_genes}, mean.cutoff = c({min_mean}, {max_mean}), dispersion.cutoff = c({min_disp}, {max_disp}), loess.span={span}) \
            \n# create plot \
            \nplot_highly_variable <- VariableFeaturePlot({object}) \
            \nplot_highly_variable
            """
            script_state.add_script(script, language=Language.R)

        if language == Language.python or language == Language.python.value or language == Language.ALL_SUPPORTED:
            if object == None:
                object = "adata"
            script = f""" \
            \n#Filter highly variable genes \
            \nsc.pp.normalize_total({object}, target_sum=1e4) \
            \nsc.pp.log1p({object}) \
            \nsc.pp.highly_variable_genes({object}, min_mean={min_mean}, max_mean={max_mean}, min_disp={min_disp}, max_disp={max_disp}, n_top_genes={n_top_genes}, span={span}) \
            \nsc.pl.highly_variable_genes({object})
            """
            script_state.add_script(script, language=Language.python)

        if not isinstance(language, Language):
            print("Error: Unknown language, not adding to script state")
            return
