from pydantic import ValidationError
import streamlit as st
import scanpy as sc
import pickle
import os
from models.AdataModel import AdataModel
from models.WorkspaceModel import WorkspaceModel
from database.database import SessionLocal
from sqlalchemy.orm import Session
from database.schemas import schemas
from utils.AdataState import AdataState

st.set_page_config(page_title='Nuwa', page_icon='🧬')

with open('css/common.css') as f:
    common_style = f"""
                <style>
                {f.read()}
                </style>
                """
    st.markdown(common_style, unsafe_allow_html=True)
    
    
class Upload:
    def __init__(self):
        self.upload_file()
        self.scanpy_dataset()
        self.external_sources()
        
    def upload_file(self):
        st.title("Upload a file")
        
        uploaded_f = st.file_uploader("Choose a file to upload", type=["csv", "h5ad", "h5", "loom", "mtx", "tsv"], accept_multiple_files=True)
        
        if not uploaded_f:
            st.text("No File uploaded 🧬")
            
        try:
            workspace_model: WorkspaceModel = st.session_state["current_workspace"]

            if uploaded_f is not None:

                with st.spinner(text="Loading Data"):
                    #create workspace dirs
                    upload_path = os.path.join(workspace_model.data_dir, "uploads")
                    download_path = os.path.join(workspace_model.data_dir, "downloads")
                    adata_path = os.path.join(workspace_model.data_dir, "adata")
                    if not os.path.exists(upload_path):
                        os.mkdir(upload_path)
                    if not os.path.exists(download_path):
                        os.mkdir(download_path)
                    if not os.path.exists(adata_path):
                        os.mkdir(adata_path)

                    for f in uploaded_f:
                        bytes_data = f.read()
                        #add to uploads dir
                        path = upload_path + f.name
                        with open(path, 'wb') as file:
                            file.write(bytes_data)

                        file_type = f.name.split(".")[-1]

                        if file_type == "mtx":
                            adata = sc.read_10x_mtx(upload_path, var_names='gene_symbols', cache=True)
                            self.show_anndata(adata, f)
                        if file_type == "h5ad":
                            adata = sc.read_h5ad(f)
                            self.show_anndata(adata, f)

        except KeyError as ke:
            print("KeyError: ", ke)
            st.error("Couldn't find workspace in session, have you selected one?")
            
        
        
    def scanpy_dataset(self):
        st.subheader("Or select a scanpy dataset")
        
        dataset_options = ['none', 'krumsiek11', 'moignard15', 'pbmc3k', 'pbmc3k_processed', 'pbmc68k_reduced', 'paul15', 'visium_sge']   
        scanpy_ds = st.selectbox(label="Dataset", options=dataset_options, key="sb_sc_datasets")
        if scanpy_ds != 'none':
            with st.spinner(text="Loading dataset"):
                if(st.session_state.sb_sc_datasets == 'krumsiek11'):
                    self.adata = sc.datasets.krumsiek11()
                if(st.session_state.sb_sc_datasets == 'moignard15'):
                    self.adata = sc.datasets.moignard15()
                if(st.session_state.sb_sc_datasets == 'pbmc3k'):
                    self.adata = sc.datasets.pbmc3k()
                if(st.session_state.sb_sc_datasets == 'pbmc3k_processed'):
                    self.adata = sc.datasets.pbmc3k_processed()
                if(st.session_state.sb_sc_datasets == 'pbmc68k_reduced'):
                    self.adata = sc.datasets.pbmc68k_reduced()
                if(st.session_state.sb_sc_datasets == 'paul15'):
                    self.adata = sc.datasets.paul15()
                if(st.session_state.sb_sc_datasets == 'visium_sge'):
                    self.adata = sc.datasets.visium_sge()
            
            self.show_anndata(self.adata)
        

    def external_sources(self):

        st.subheader("External sources")
        with st.form(key="ebi_ds_form"):
            st.subheader("EBI Expression Atlas")
            accession_str = st.text_input(label="EBI accession string")
            ebi_form_btn = st.form_submit_button("Find dataset")
            if ebi_form_btn:
                with st.spinner(text="Fetching dataset"):
                    dataset = sc.datasets.ebi_expression_atlas(accession=accession_str)
                    self.show_anndata(dataset)


    def show_anndata(self, adata, f = None):
        try:
            #upload raw adata
            sc.write(filename=f"{os.getenv('WORKDIR')}uploads/adata_raw.h5ad", adata=adata)
            sc.write(filename=f"{os.getenv('WORKDIR')}adata/adata_raw.h5ad", adata=adata)

            st.session_state["adata_state"] = AdataState(workspace_id=st.session_state.current_workspace.id)

        except ValidationError as e:
            st.error(e)

        with st.sidebar:
            st.subheader("File info")
            st.write(f"AnnData object with n_obs x n_vars = {adata.n_obs} x {adata.n_vars}")
            if f is not None:
                st.write(f"Size: {f.size} bytes")

            st.subheader("Obs")
            if not adata.obs.empty:
                st.dataframe(adata.obs.head())
                st.write(f"Showing 5 rows of {len(adata.obs.columns)} columns")
            else:
                st.text("No obs to show")
            
            st.subheader("Var")
            if not adata.var.empty:
                st.dataframe(adata.var.head())
                st.write(f"Showing 5 rows of {len(adata.var.columns)} columns")
            else:
                st.text("No var to show")

        st.toast("Successfully uploaded file", icon='✅')
        
        
upload_page = Upload()



