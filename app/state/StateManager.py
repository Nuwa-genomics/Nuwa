from state.AdataState import AdataState
from state.ScriptState import ScriptState
from models.AdataModel import AdataModel
from models.ScriptModel import ScriptModel
from models.WorkspaceModel import WorkspaceModel
from scripts.Script import Script
from utils.session_cache import cache_data_to_session, load_data_from_cache
import scanpy as sc
import os
import streamlit as st
from anndata import AnnData

class StateManager:
    """
    Makes changes made to the database and filesystems synchronously using data from session state in an atomic way. Also responsible for loading data into session state and initialising files when loading new dataset.
    """

    def add_script(self, script: Script):
        # add script if present
        if script is not None:
            if isinstance(script, Script):
                self.script = script
        return self

    def add_adata(self, adata: AnnData):
        if isinstance(adata, AnnData):
            self.add_adata = adata
        return self

    def load_session():
        raise NotImplementedError
    

    def save_session(self):
        """
        Takes a copy of current session state and saves in pickle format. Extracts adata from session state and updates filesystem and database. Also
        takes an optional script to update database.

        """ 
        # write adata h5ad object to file
        if hasattr(self, 'adata'):
            sc.write(filename=os.path.join(os.getenv('WORKDIR'), 'adata', st.session_state.adata_state.current.adata_name), adata=self.adata)
        
        # add script if present
        if hasattr(self, 'script'):
            self.script.add_script()

        # cache data to pickle file
        cache_data_to_session()


    def init_session():
        raise NotImplementedError