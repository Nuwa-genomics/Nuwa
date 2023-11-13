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


st.title("Upload a file")


uploaded_f = st.file_uploader("Choose a file to upload", type=["csv", "h5ad", "h5", "loom", "mtx", "tsv"], accept_multiple_files=True)


if not uploaded_f:
    st.text("No File uploaded 🧬")
    

def show_anndata(adata):
    try:
        st.session_state["adata_state"] = AdataState(workspace_id=st.session_state.current_workspace.id)

        #upload raw adata
        sc.write(filename=os.path.join(os.getenv('WORKDIR'), 'uploads/', "adata_raw.h5ad"), adata=adata)

        st.session_state.adata_state.add_adata(
            AdataModel(work_id=st.session_state.current_workspace.id, adata_name="adata_raw", 
                filename=os.path.join(os.getenv('WORKDIR'), "adata/", "adata_raw.h5ad"), adata=adata, notes=""))
    except ValidationError as e:
        st.error(e)

    st.subheader("File info")
    st.text(f"AnnData object with n_obs x n_vars = {adata.n_obs} x {adata.n_vars}")
    st.text(f"Size: {f.size} bytes")

    st.subheader("Obs")
    if not adata.obs.empty:
        st.dataframe(adata.obs.head())
        st.text(f"Showing 5 rows of {len(adata.obs.columns)} columns")
    else:
        st.text("No obs to show")
    
    st.subheader("Var")
    if not adata.var.empty:
        st.dataframe(adata.var.head())
        st.text(f"Showing 5 rows of {len(adata.var.columns)} columns")
    else:
        st.text("No var to show")

    st.toast("Successfully uploaded file", icon='✅')


try:
    workspace_model: WorkspaceModel = st.session_state["current_workspace"]

    if uploaded_f is not None:
        

        with st.spinner(text="Loading Data"):
            #create workspace dirs
            upload_path = f'{workspace_model.data_dir}/uploads/'
            download_path = f'{workspace_model.data_dir}/downloads/'
            adata_path = f'{workspace_model.data_dir}/adata/'
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
                    show_anndata(adata)
                if file_type == "h5ad":
                    adata = sc.read_h5ad(f)
                    show_anndata(adata)

except KeyError as ke:
    print("KeyError: ", ke)
    st.error("Couldn't find workspace in session, have you selected one?")