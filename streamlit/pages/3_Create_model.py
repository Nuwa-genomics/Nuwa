import streamlit as st
from ml.citeseq.model import CiteAutoencoder
from ml.citeseq.dataset import TabularDataset
from sklearn.model_selection import train_test_split
from torch.utils.data import DataLoader
import numpy as np
import pickle
import scvi
from ml.solo_scvi.solo_model import *
from ml.DeepST.deepst.main import *
import torch
from components.sidebar import *
from models.AdataModel import AdataModel


st.set_page_config(layout="wide", page_title='Nuwa', page_icon='🧬')

with open('css/common.css') as f:
    common_style = f"""
                <style>
                {f.read()}
                </style>
                """
    st.markdown(common_style, unsafe_allow_html=True)

try:
    adata_model = st.session_state["adata"]
    show_sidebar(adata_model)
except KeyError as ke:
    print('Key Not Found in Employee Dictionary:', ke)


adata = get_adata(adataList=adata_model, name=st.session_state.sb_adata_selection).adata
st.session_state["current_adata"] = adata


class CreateCiteSeqModel:
    def __init__(self, adata):
        self.adata = adata
        

    def draw_page(self):
        col1, _, _, _ = st.columns(4)

        with col1:
            st.selectbox(label="Device", options=(self.device_options), on_change=self.set_device, key="sb_select_device_citeseq")

        col1, col2 = st.columns(2, gap="large")

        with col1:
            self.set_hyperparams()
            self.set_train_test_split()
        with col2:
            self.set_autoencoder()
        
    def init_device(self):
        #init devices if not already exists
        if torch.cuda.is_available():
            self.device_options = ["Cuda (GPU)", "CPU"]
        else:
            self.device_options = ["CPU"]
        
        if not 'sb_select_device_citeseq' in st.session_state:
            self.device = "cuda" if torch.cuda.is_available() else "cpu"
        else:
            self.device = "cpu" if st.session_state.sb_select_device_citeseq == "CPU" else "cuda"

        st.session_state.device = self.device

    def set_device(self):
        
        self.device = "cpu" if st.session_state["sb_select_device_citeseq"] == "CPU" else "cuda"
        st.session_state["device"] = self.device

    def init_model(self):
        #initialize model object
        self.model_dict = {
            "model": None,
            "lr": 1e-2,
            "n_epochs": 100,
            "n_features": self.adata.to_df().shape[1],
            "optim": 'Adam',
            "test_split": 0.1,
            "train_dl": None,
            "valid_dl": None
        }

        st.session_state['model_obj'] = self.model_dict

    def create_datasets(self, df, test_size=None):

        if test_size == None:
            #get value from input
            input_test_value = st.session_state["input_train_test_split"]
            test_size = input_test_value / 100

            if test_size > 0.9:
                st.warning("This Train/Test split may be too high leading to overfitting")
            elif test_size < 0.5:
                st.warning("The train/test split may be too low to provide enough training data")
        
        train_data, valid_data = train_test_split(df.to_numpy(dtype=np.float32), test_size=test_size)

        train_ds = TabularDataset(train_data)
        valid_ds = TabularDataset(valid_data)

        train_dl = DataLoader(train_ds, batch_size=64, shuffle=True)
        valid_dl = DataLoader(valid_ds, batch_size=64, shuffle=False)

        st.session_state.model_obj["train_dl"] = train_dl
        st.session_state.model_obj["valid_dl"] = valid_dl
        st.session_state.model_obj["test_split"] = test_size

        return train_dl, valid_dl
    
    def change_hyperparams(self):
        st.session_state.model_obj["lr"] = st.session_state.input_lr
        st.session_state.model_obj["n_epochs"] = st.session_state.input_nepochs
        st.session_state.model_obj["optim"] = st.session_state.input_optim

    def set_hyperparams(self):
        st.subheader("Set model hyperparameters")
        st.number_input(label="Epochs", min_value=1, key="input_nepochs", value=100, on_change=self.change_hyperparams)
        st.number_input(label="Learning rate", min_value=1e-4, max_value=1.0, value=1e-3, key="input_lr", step=1e-3, format='%.4f', on_change=self.change_hyperparams)
        st.selectbox(label="Optimizer", options=(["Adam", "SGD", "RMSProp"]), key="input_optim", on_change=self.change_hyperparams)

    def set_autoencoder(self):
        st.subheader("Autoencoder")
        st.number_input(label="Latent Dimensions", key="input_latent", min_value=2, value=30)
        st.number_input(label="Hidden Layers", key="input_hidden", min_value=1, value=100)

    def set_train_test_split(self):
        st.subheader("Train/test split")
        st.slider(label=f"Train data %", min_value=1, max_value=99, value=90, key="input_train_test_split", on_change=self.create_datasets)

    def build_model(self, df):

        n_features = df.shape[1]

        model = CiteAutoencoder(
            nfeatures_rna=n_features, 
            nfeatures_pro=0, 
            hidden_rna=100, 
            hidden_pro=0, 
            z_dim=50
        )

        # Create train and test datasets
        train_data, valid_data = train_test_split(df.to_numpy(dtype=np.float32), test_size=0.1)

        train_ds = TabularDataset(train_data)
        valid_ds = TabularDataset(valid_data)

        train_dl = DataLoader(train_ds, batch_size=64, shuffle=True)
        valid_dl = DataLoader(valid_ds, batch_size=64, shuffle=False)

        st.session_state.model_obj["model"] = model
        st.session_state.model_obj["n_features"] = n_features

class CreateSoloModel:
    def __init__(self, adata):
        self.adata = adata

    def draw_page(self):

        col1, _, _, _ = st.columns(4)

        with col1:
            st.selectbox(label="Device", options=(self.device_options), on_change=self.set_device, key="sb_select_device_solo")

        col1, col2 = st.columns(2)
        with col1:
            st.subheader("Model parameters")
            st.number_input(label="Epochs", min_value=1, value=400, key="ni_vae_epochs")
            st.number_input(label="Learning rate", min_value=1e-4, max_value=10.0, value=1e-3, format='%.4f', key="ni_solo_lr")
            st.subheader("Train/test split")
            st.slider(label=f"Train data %", min_value=1, max_value=99, value=90, key="input_train_test_split_solo_vae")

    def init_device(self):
        #init devices if not already exists
        if torch.cuda.is_available():
            self.device_options = ["Cuda (GPU)", "CPU"]
        else:
            self.device_options = ["CPU"]
        
        if not 'sb_select_device_solo' in st.session_state:
            self.device = "cuda" if torch.cuda.is_available() else "cpu"
        else:
            self.device = "cpu" if st.session_state.sb_select_device_solo == "CPU" else "cuda"

        st.session_state.device = self.device

    def set_device(self):
        
        self.device = "cpu" if st.session_state["sb_select_device_solo"] == "CPU" else "cuda"
        st.session_state["device"] = self.device

    def init_model(self):
        #initialize model object
        self.model_dict = {
            "model": None,
            "n_epochs": 400,
            "test_split": 0.1,
            "lr": 1e-3,
        }

        st.session_state['model_obj'] = self.model_dict

    def build_model(self, adata):
        model = solo_model(adata)
        st.session_state.model_obj["model"] = model

    def change_hyperparams(self):
        st.session_state.model_obj["lr"] = 1e-3
        st.session_state.model_obj["n_epochs"] = 400

class CreateDeepSTModel:
    def __init__(self, adata):
        self.adata = adata

    def draw_page(self):

        col1, col2, _, _ = st.columns(4, gap="large")
        col1.selectbox(label="Device", options=(self.device_options), on_change=self.set_device, key="sb_select_device_deepst")

        col1, col2 = st.columns(2, gap="large")

        with col1:
            st.subheader("Set model hyperparameters")
            st.number_input(label="Pre-epochs", min_value=1, key="ni_deepst_preepochs", value=800)
            st.number_input(label="Epochs", min_value=1, key="ni_deepst_epochs", value=1000)
            st.slider(label="No. Principal components", min_value=1, max_value=500, value=200, key="sl_n_pca")
            
        
        with col2:
            st.subheader("Spatial parameters")
            st.number_input(label="Number of spatial domains", min_value=1, key="ni_sp_domains", value=10)
            st.selectbox(label="Task", options=(['Identify_Domain', 'Integration']), key="sb_deepst_task")
            st.selectbox(label="Platform", options=['Visium'], key="sb_deepst_platform")
            st.selectbox(label="Disttype", options=['BallTree', 'KDTree', 'kneighbors_graph', 'radius'], key="sb_deepst_disttype")

    def init_device(self):
        #init devices if not already exists
        if torch.cuda.is_available():
            self.device_options = ["Cuda (GPU)", "CPU"]
        else:
            self.device_options = ["CPU"]
        
        if not 'sb_select_device_deepst' in st.session_state:
            self.device = "cuda" if torch.cuda.is_available() else "cpu"
        else:
            self.device = "cpu" if st.session_state.sb_select_device_deepst == "CPU" else "cuda"

        st.session_state.device = self.device

    def set_device(self):
        
        self.device = "cpu" if st.session_state["sb_select_device_deepst"] == "CPU" else "cuda"
        st.session_state["device"] = self.device

    def init_model(self):
        #initialize model object
        self.model_dict = {
            "model": None,
            "n_preepochs": 800,
            "n_epochs": 1000,
            "n_domains": 10,
            "platform": "visium",
            "task": "Identify_Domain",
            "dist type": "BallTree",
            "n_pca_comp": 200,
        }

        st.session_state['model_obj'] = self.model_dict

    def build_model(self, adata):
        model = DeepSTModel(
            task=st.session_state.sb_deepst_task,
            n_domains=st.session_state.ni_sp_domains,
            pre_epochs=st.session_state.ni_deepst_preepochs,
            epochs=st.session_state.ni_deepst_epochs,
            use_gpu=(self.device == "cuda"),
            platform=st.session_state.sb_deepst_platform,
            dist_type=st.session_state.sb_deepst_disttype,
            n_pca_comps=st.session_state.sl_n_pca
        )
        st.session_state.model_obj["model"] = model
        


def create_citeseq(adata):
    create_model = CreateCiteSeqModel(adata)

    create_model.init_device()

    create_model.draw_page()

    create_model.init_model()

    create_model.build_model(adata.to_df())

    create_model.create_datasets(adata.to_df(), test_size=st.session_state.model_obj["test_split"])

    st.subheader("Model summary")
    st.json(st.session_state.model_obj, expanded=False)

def create_solo(adata):
    create_model = CreateSoloModel(adata)

    create_model.init_device()

    create_model.draw_page()

    create_model.init_model()

    create_model.build_model(adata)

    st.subheader("Model summary")
    st.json(st.session_state.model_obj, expanded=False)

def create_deepst(adata):
    create_model = CreateDeepSTModel(adata)

    create_model.init_device()

    create_model.draw_page()

    create_model.init_model()

    create_model.build_model(adata)

    st.subheader("Model summary")
    st.json(st.session_state.model_obj, expanded=False)


def change_model():
    adata = st.session_state["current_adata"]
    if st.session_state.sb_model_selection == 'Citeseq (dimensionality reduction)':
        create_citeseq(adata)
    elif st.session_state.sb_model_selection == 'Solo (doublet removal)':
        create_solo(adata)
    elif st.session_state.sb_model_selection == 'DeepST (identify spatial domains)':
        create_deepst(adata)

st.title("Create Model")


col1, _, _, _ = st.columns(4)
col1.selectbox(label="model", options=([
    "Citeseq (dimensionality reduction)", 
    "Solo (doublet removal)",
    "DeepST (identify spatial domains)"
    ]), key='sb_model_selection')

change_model()

show_preview()