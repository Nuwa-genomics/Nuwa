import streamlit as st
import scanpy as sc
import math
import time
import pickle
import os

common_style = """
            <style>
            footer {visibility: hidden;}
            .st-emotion-cache-1cypcdb {background: linear-gradient(180deg, rgb(5, 39, 103) 0%, #3a0647 70%); box-shadow: 1px 0 10px -2px #000;}
            .st-emotion-cache-86cver {rgba(250, 250, 250, 0.6)}
            </style>
            """
st.markdown(common_style, unsafe_allow_html=True)


st.title("Upload a file")


uploaded_f = st.file_uploader("Choose a file to upload", type=["csv", "h5ad", "loom", "mtx", "tsv"], accept_multiple_files=True)


if not uploaded_f:
    st.text("No File uploaded 🧬")
else:
    st.subheader("File info")
    st.toast("Successfully uploaded file", icon='✅')


def show_anndata(adata):
    st.text(f"AnnData object with n_obs x n_vars = {adata.n_obs} x {adata.n_vars}")
    st.text(f"Size: {f.size} bytes")
    st.subheader("Obs")
    st.dataframe(adata.obs.head())
    st.text(f"Showing 5 rows of {len(adata.obs.columns)} columns")
    st.subheader("Var")
    st.dataframe(adata.var.head())
    st.text(f"Showing 5 rows of {len(adata.var.columns)} columns")
    st.session_state["adata"] = adata
    with open('./tmp/adata.pkl', 'wb') as tmp:
        pickle.dump(adata, tmp, pickle.HIGHEST_PROTOCOL)


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