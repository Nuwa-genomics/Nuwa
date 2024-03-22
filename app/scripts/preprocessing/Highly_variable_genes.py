from enums.Language import Language
import streamlit as st
from state.ScriptState import ScriptState
import numpy as np
from scripts.Script import Script

class Highly_variable_genes(Script):
    """
    Exports a python or R script for computing highly variable genes.
    """

    def __init__(self, language: Language | str, 
        min_mean: float, 
        max_mean: float, 
        min_disp: float = 0.5, 
        max_disp=np.Inf, 
        n_top_genes: int = 2000, 
        span=0.3, 
        object: str = None):

        super().__init__(language=language)
        
        self.min_mean = min_mean
        self.max_mean = max_mean
        self.min_disp = min_disp
        self.max_disp = max_disp
        self.n_top_genes = n_top_genes
        self.span = span
        self.object = object


    def add_script(self):

        if self.language == Language.R or self.language == Language.R.value or self.language == Language.ALL_SUPPORTED:
            if self.object == None:
                self.object = "pbmc.data"
            script = f""" \
            \n#Compute highly variable genes. This uses the FindVariableFeatures method in seurat \
            \n#If you have been using bioconductor's sce object you will need to reload it in as a seurat object. \
            \n{self.object} <- Read10X(data.dir = 'path_to_mtx_files') \
            \n{self.object} <- CreateSeuratObject(counts = {self.object}, project = 'pbmc3k') \
            \n#Log normalize \
            \n{self.object} <- NormalizeData({self.object}, normalization.method = 'LogNormalize', scale.factor = 10000) \
            \n{self.object} <- FindVariableFeatures({self.object}, selection.method = 'dispersion', nfeatures = {self.n_top_genes}, mean.cutoff = c({self.min_mean}, {self.max_mean}), dispersion.cutoff = c({self.min_disp}, {self.max_disp}), loess.span={self.span}) \
            \n# create plot \
            \nplot_highly_variable <- VariableFeaturePlot({self.object}) \
            \nplot_highly_variable
            """
            self.script_state.add_script(script, language=Language.R)

        if self.language == Language.python or self.language == Language.python.value or self.language == Language.ALL_SUPPORTED:
            if self.object == None:
                self.object = "adata"
            script = f""" \
            \n#Filter highly variable genes \
            \nsc.pp.normalize_total({self.object}, target_sum=1e4) \
            \nsc.pp.log1p({self.object}) \
            \nsc.pp.highly_variable_genes({self.object}, min_mean={self.min_mean}, max_mean={self.max_mean}, min_disp={self.min_disp}, max_disp={self.max_disp}, n_top_genes={self.n_top_genes}, span={self.span}) \
            \nsc.pl.highly_variable_genes({self.object})
            """
            self.script_state.add_script(script, language=Language.python)
