from pydantic import ValidationError
import streamlit as st
import scanpy as sc
import pickle
import os
from models.AdataModel import AdataModel

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
        st.session_state["adata"] = [adata_model]
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

def show_csv_file_info():
    st.text("This is a csv file")

if uploaded_f is not None:
    #create uploads dir
    if not os.path.exists('./uploads'):
        os.mkdir('./uploads')

    for f in uploaded_f:
        bytes_data = f.read()
        path = './uploads/' + f.name
        with open(path, 'wb') as file:
            file.write(bytes_data)

        file_type = f.name.split(".")[-1]

        if file_type == "mtx":
            show_mtx_file_info('./uploads')
        if file_type == "h5ad":
            show_h5ad_file_info(f)