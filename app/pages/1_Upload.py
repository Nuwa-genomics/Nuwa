from pydantic import ValidationError
import streamlit as st
import scanpy as sc
import squidpy as sq
import pickle
import os
from models.AdataModel import AdataModel
from models.WorkspaceModel import WorkspaceModel
from database.database import SessionLocal
from sqlalchemy.orm import Session
from database.schemas import schemas
from utils.AdataState import AdataState

st.set_page_config(page_title='Nuwa', page_icon='🧬')

os.chdir('/app')

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
        st.title("Upload a dataset")
        
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
        st.subheader("Or select a built-in dataset")
        
        #scanpy datasets
        dataset_options = ['none', 'pbmc3k', 'pbmc3k_processed', 'pbmc68k_reduced', 'paul15', 'four_i', 'imc', 'seqfish', 'merfish', 'mibitof', 'slideseqv2', 'sc_mouse_cortex', 'visium', 'visium_hne_adata', 'visium_hne_adata_crop', 'visium_fluo_adata', 'visium_fluo_adata_crop', 'visium_hne_image', 'visium_hne_crop', 'visium_fluo_image_crop']   
        scanpy_ds = st.selectbox(label="Dataset", options=dataset_options, key="sb_sc_datasets")
        ds_empty = st.empty()
        
        if scanpy_ds != 'none':
            with st.spinner(text="Loading dataset"):
                if(st.session_state.sb_sc_datasets == 'pbmc3k'):
                    self.adata = sc.datasets.pbmc3k()
                    ds_empty.empty()
                    ds_empty.info("3k PBMCs from 10x Genomics. The data consists in 3k PBMCs from a Healthy Donor and is available from 10x Genomics [here](%s)" % "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz")
                if(st.session_state.sb_sc_datasets == 'pbmc3k_processed'):
                    self.adata = sc.datasets.pbmc3k_processed()
                    ds_empty.empty()
                    ds_empty.info("Processed 3k PBMCs from 10x Genomics.")
                if(st.session_state.sb_sc_datasets == 'pbmc68k_reduced'):
                    self.adata = sc.datasets.pbmc68k_reduced()
                    ds_empty.empty()
                    ds_empty.info("Subsampled and processed 68k PBMCs. 10x PBMC 68k dataset from [here](%s). The original PBMC 68k dataset was preprocessed using scanpy and was saved keeping only 724 cells and 221 highly variable genes. The saved file contains the annotation of cell types (key: 'bulk_labels'), UMAP coordinates, louvain clustering and gene rankings based on the bulk_labels." % "https://support.10xgenomics.com/single-cell-gene-expression/datasets")
                if(st.session_state.sb_sc_datasets == 'paul15'):
                    self.adata = sc.datasets.paul15()
                    ds_empty.empty()
                    ds_empty.info("Development of Myeloid Progenitors. Non-logarithmized raw data.")
                if(st.session_state.sb_sc_datasets == 'four_i'):
                    self.adata = sq.datasets.four_i()
                    ds_empty.empty()
                    ds_empty.info("Pre-processed subset 4i dataset from [Gut et al](%s)." % "https://doi.org/10.1126/science.aar7042")
                if(st.session_state.sb_sc_datasets == 'imc'):
                    self.adata = sq.datasets.imc()
                    ds_empty.empty()
                    ds_empty.info("Pre-processed subset IMC dataset from [Jackson et al](%s)." % "https://www.nature.com/articles/s41586-019-1876-x")
                if(st.session_state.sb_sc_datasets == 'seqfish'):
                    self.adata = sq.datasets.seqfish()
                    ds_empty.empty()
                    ds_empty.info("Pre-processed subset seqFISH dataset from [Lohoff et al](%s)." % "https://www.biorxiv.org/content/10.1101/2020.11.20.391896v1")
                if(st.session_state.sb_sc_datasets == 'merfish'):
                    self.adata = sq.datasets.merfish()
                    ds_empty.empty()
                    ds_empty.info("Pre-processed MERFISH dataset from [Moffitt et al](%s)." % "https://doi.org/10.1126/science.aau5324")
                if(st.session_state.sb_sc_datasets == 'mibitof'):
                    self.adata = sq.datasets.mibitof()
                    ds_empty.empty()
                    ds_empty.info("Pre-processed MIBI-TOF dataset from [Hartmann et al](%s)." % "https://doi.org/10.1101/2020.01.17.909796")
                if(st.session_state.sb_sc_datasets == 'slideseqv2'):
                    self.adata = sq.datasets.slideseqv2()
                    ds_empty.empty()
                    ds_empty.info("Pre-processed SlideseqV2 dataset from [Stickles et al](%s)." % "https://doi.org/10.1038/s41587-020-0739-1")
                if(st.session_state.sb_sc_datasets == 'sc_mouse_cortex'):
                    self.adata = sq.datasets.sc_mouse_cortex()
                    ds_empty.empty()
                    ds_empty.info("Pre-processed [scRNA-seq mouse cortex](%s)." % "https://doi.org/10.1038/s41586-018-0654-5")
                if(st.session_state.sb_sc_datasets == 'visium'):
                    self.adata = sq.datasets.visium()
                    ds_empty.empty()
                    ds_empty.info("Download Visium [datasets](%s) from 10x Genomics." % "https://support.10xgenomics.com/spatial-gene-expression/datasets")
                if(st.session_state.sb_sc_datasets == 'visium_hne_adata'):
                    self.adata = sq.datasets.visium_hne_adata()
                    ds_empty.empty()
                    ds_empty.info("Pre-processed 10x Genomics Visium H&E [dataset](%s)." % "https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Adult_Mouse_Brain")
                if(st.session_state.sb_sc_datasets == 'visium_hne_adata_crop'):
                    self.adata = sq.datasets.visium_hne_adata_crop()
                    ds_empty.empty()
                    ds_empty.info("Pre-processed subset 10x Genomics Visium H&E [dataset](%s)." % "https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Adult_Mouse_Brain")
                if(st.session_state.sb_sc_datasets == 'visium_fluo_adata'):
                    self.adata = sq.datasets.visium_fluo_adata()
                    ds_empty.empty()
                    ds_empty.info("Pre-processed 10x Genomics Visium Fluorecent [dataset](%s)." % "https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Adult_Mouse_Brain_Coronal_Section_2")
                if(st.session_state.sb_sc_datasets == 'visium_fluo_adata_crop'):
                    self.adata = sq.datasets.visium_fluo_adata_crop()
                    ds_empty.empty()
                    ds_empty.info("Pre-processed subset 10x Genomics Visium Fluorescent [dataset](%s)." % "https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Adult_Mouse_Brain_Coronal_Section_2")

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
            
            #make var unique
            
            if f is not None:
                st.write(f"Size: {f.size} bytes")
                
            

            st.subheader("Obs")
            
            if not adata.obs.empty:
                make_obs_unique = st.checkbox(label="Obs names make unique", key="btn_make_obs_unique", value=False)
                if make_obs_unique:
                    st.write("to implement")
                st.dataframe(adata.obs.head())
                st.write(f"Showing 5 rows of {len(adata.obs.columns)} columns")
                
            else:
                st.text("No obs to show")
            
            st.subheader("Var")
            
            if not adata.var.empty:
                make_var_unique = st.checkbox(label="Var names make unique", key="btn_make_var_unique", value=False)
                if make_var_unique:
                    st.write("To implement")
                st.dataframe(adata.var.head())
                st.write(f"Showing 5 rows of {len(adata.var.columns)} columns")
                
            else:
                st.text("No var to show")

        st.toast("Successfully uploaded file", icon='✅')
        
        
upload_page = Upload()



