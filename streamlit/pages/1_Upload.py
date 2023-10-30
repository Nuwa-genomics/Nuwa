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
else:
    st.subheader("File info")
    st.toast("Successfully uploaded file", icon='✅')


def show_anndata(adata):

    #store adata in session storage
    try:
        adata_model = AdataModel(
            work_id=st.session_state.current_workspace.id, 
            id=0, 
            adata_name="adata_raw", 
            filename="adata_raw.h5ad",
            adata=adata
        )
        
        #fetch from db
        conn: Session = SessionLocal()
        adatas = conn.query(schemas.Adata).all()
        st.session_state["adata"] = adatas
        st.session_state["adata"].append(adata_model)
    except ValidationError as e:
        st.error(e)

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


def show_mtx_file_info(p):
    with st.spinner(text="Loading data files"):
        adata = sc.read_10x_mtx(p, var_names='gene_symbols', cache=True)
        show_anndata(adata)
        
def show_h5ad_file_info(f):
    with st.spinner(text="Loading data files"):
        adata = sc.read_h5ad(f)
        show_anndata(adata)



try:

    workspace_model: WorkspaceModel = st.session_state["current_workspace"]

    if uploaded_f is not None:
        #create uploads dir
        upload_path = f'{workspace_model.data_dir}/uploads/'
        if not os.path.exists(upload_path):
            os.mkdir(upload_path)

        for f in uploaded_f:
            bytes_data = f.read()
            path = upload_path + f.name
            with open(path, 'wb') as file:
                file.write(bytes_data)

            file_type = f.name.split(".")[-1]

            if file_type == "mtx":
                show_mtx_file_info(upload_path)
            if file_type == "h5ad":
                show_h5ad_file_info(f)

except KeyError as ke:
    print("KeyError: ", ke)
    st.error("Couldn't find workspace in session, have you selected one?")