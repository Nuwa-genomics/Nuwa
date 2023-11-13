from streamlit.testing.v1 import AppTest
import streamlit as st
import scanpy as sc
import os

from utils.AdataState import AdataState

class Test_Upload:
    def __init__(self, session_state = None):
        self.at = AppTest.from_file("pages/1_Upload.py")
        if session_state is not None:
            self.at.session_state = session_state

        self.adata = sc.datasets.pbmc3k()
        
        self.at.run()
        self.test_load_adata()
        self.test_upload_data()
        assert not self.at.exception

    def test_load_adata(self):
        self.at.session_state["adata_state"] = AdataState(workspace_id=self.at.session_state.current_workspace.id)
        
    def test_upload_data(self):
        path = os.path.join(os.getenv('WORKDIR'), 'uploads/', "adata_raw.h5ad")
        sc.write(filename=path, adata=self.adata)
        file_read = sc.read_h5ad(path)
        assert file_read.obs.to_dict() == self.adata.obs.to_dict()
        assert file_read.var.to_dict() == self.adata.var.to_dict()

    def get_final_session_state(self):
        return self.at.session_state




        