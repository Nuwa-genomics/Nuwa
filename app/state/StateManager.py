from typing import List
from state.AdataState import AdataState
from state.ScriptState import ScriptState
from models.AdataModel import AdataModel
from models.ScriptModel import ScriptModel
from models.WorkspaceModel import WorkspaceModel
from scripts.Script import Script
from utils.session_cache import cache_data_to_session, load_data_from_cache
import scanpy as sc
from database.schemas import schemas
import os
import streamlit as st
from anndata import AnnData
from database.database import SessionLocal
from enums.ErrorMessage import ErrorMessage

class StateManager:
    """
    Makes changes made to the database and filesystems synchronously using data from session state in an atomic way. Also responsible for loading data into session state and initialising files when loading new dataset.
    """

    ########## Factory methods ##########

    def add_script(self, script: Script):
        if isinstance(script, Script):
            self.script = script
        return self

    def add_adata(self, adata: AnnData):
        if isinstance(adata, AnnData):
            self.adata = adata
        return self
    
    def add_description(self, description: str):
        self.description = description
        return self
    
    ########## session ##########

    def load_session(self):

        # get current adata id either from session state or environment
        if "adata_state" in st.session_state:
            current_adata_id = st.session_state.adata_state.current.id
        else:
            current_adata_id = os.getenv('CURRENT_ADATA_ID')

        # fetch from database
        conn = SessionLocal()
        cache_files = conn.query(schemas.Session) \
            .filter(schemas.Session.adata_id == int(current_adata_id)) \
            .all() \
            .sort(schemas.Session.created)
        
        st.write(cache_files)
        
        #load_data_from_cache(cache_file)
    

    def save_session(self):
        """
        Takes a copy of current session state and saves in pickle format. Extracts adata from session state and updates filesystem and database. Also
        takes an optional script to update database.

        """ 
        # write adata h5ad object to file
        if hasattr(self, 'adata'):
            sc.write(filename=os.path.join(os.getenv('WORKDIR'), 'adata', self.adata_state().current.adata_name), adata=self.adata)
            st.session_state.adata_state.current.adata = self.adata
        
        # add script if present
        if hasattr(self, 'script'):
            self.script.add_script()

        if not hasattr(self, 'description'):
            self.description = ""


        # cache data to pickle file
        cache_data_to_session(description=self.description)

    
    ########## Adata state ##########
    
    def get_current_adata(self) -> AnnData:
        """Return a copy of current adata."""
        return st.session_state.adata_state.current.adata.copy()
    
    def adata_state(self) -> AdataState:
        """An accessor for the adata state."""
        return st.session_state.adata_state
    
    ########## Script state ##########

    def script_state(self) -> ScriptState:
        """An accessor for the script state."""
        return st.session_state.script_state
    
    ########## Workspace ##########

    def get_current_workspace(self) -> WorkspaceModel:
        """Get current workspace."""
        return st.session_state.current_workspace
    
    def workspaces(self) -> List[WorkspaceModel]:
        """An accessor for workspaces."""
        return st.session_state.workspaces
    
        