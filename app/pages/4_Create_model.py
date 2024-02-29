import os
os.chdir('/app')

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
from state.AdataState import AdataState
from state.StateManager import StateManager


st.set_page_config(layout="wide", page_title='Nuwa', page_icon='🧬')

with open('css/common.css') as f:
    common_style = f"""
                <style>
                {f.read()}
                </style>
                """
    st.markdown(common_style, unsafe_allow_html=True)


class Create_CiteSeq_model:
    """
    Train a deep autoencoder model for dimensionality reduction on gene expression data. This is an alternative to other DR algorithms such as tSNE. Original Github repo can be found [here](https://github.com/naity/citeseq_autoencoder)
    
    Notes
    -----
    .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/create_model_citeseq_page.png
    """
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
        """
        Init Cite-seq model with given hyperparams.

        Parameters
        ----------
        model: CiteAutoencoder
            CiteAutoencoder model object.

        device: str
            Device to run the training on (CPU or Cuda if GPU is available).

        learning_rate: float
            Adjustable learning rate for improving gradient descent optimizer.

        epochs: int
            Number of epochs to do the training.

        optimizer: str
            Optimizer for neural network.

        train_split: int
            Percntage of dataset to be used for training (not testing/validation).

        latent_dimensions: int
            Dimensionality of latent space in autoencoder.

        hidden_layers: int
            Number of hidden layers in autoencoder model.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/citeseq_model.png

        Example
        -------
        """
        self.model_dict = {
            "model": None,
            "lr": st.session_state.ni_citeseq_lr,
            "n_epochs": st.session_state.ni_citeseq_epochs,
            "n_features": self.adata.to_df().shape[1],
            "optim": st.session_state.sb_citeseq_optim,
            "test_split": round(1 - (st.session_state.citeseq_train_test_split / 100), 2),
            "train_dl": None,
            "valid_dl": None
        }

        st.session_state['model_obj'] = self.model_dict

    def create_datasets(self, test_size=None):
 
        df = self.adata.to_df()

        if test_size == None:
            #get value from input
            input_test_value = st.session_state["citeseq_train_test_split"]
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
        st.session_state.model_obj["lr"] = st.session_state.ni_citeseq_lr
        st.session_state.model_obj["n_epochs"] = st.session_state.ni_citeseq_epochs
        st.session_state.model_obj["optim"] = st.session_state.sb_citeseq_optim

    def set_hyperparams(self):
        st.subheader("Set model hyperparameters")
        st.number_input(label="Epochs", min_value=1, key="ni_citeseq_epochs", value=100, on_change=self.change_hyperparams)
        st.number_input(label="Learning rate", min_value=1e-4, max_value=1.0, value=1e-3, key="ni_citeseq_lr", step=1e-3, format='%.4f', on_change=self.change_hyperparams)
        st.selectbox(label="Optimizer", options=(["Adam", "SGD", "RMSProp"]), key="sb_citeseq_optim", on_change=self.change_hyperparams)

    def set_autoencoder(self):
        st.subheader("Autoencoder")
        st.number_input(label="Latent Dimensions", key="input_latent", min_value=2, value=30)
        st.number_input(label="Hidden Layers", key="input_hidden", min_value=1, value=100)

    def set_train_test_split(self):
        st.subheader("Train/test split")
        st.slider(label=f"Train data %", min_value=1, max_value=99, value=90, key="citeseq_train_test_split", on_change=self.create_datasets)

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

class Create_Solo_model:
    """
    Create a Solo model for doublet detection and removal. A package in [scvi tools](https://docs.scvi-tools.org/en/stable/user_guide/models/solo.html)

    Notes
    -----
    .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/create_model_solo_page.png
    """
    def __init__(self, adata):
        self.adata = adata

    def draw_page(self):

        col1, _, _, _ = st.columns(4)

        with col1:
            st.selectbox(label="Device", options=(self.device_options), on_change=self.set_device, key="sb_select_device_solo")

        col1, col2 = st.columns(2)
        with col1:
            st.subheader("Model parameters")
            st.number_input(label="Epochs", min_value=1, value=400, key="ni_vae_epochs", on_change=self.change_hyperparams)
            st.number_input(label="Learning rate", min_value=1e-4, max_value=10.0, value=1e-3, format='%.4f', key="ni_solo_lr", on_change=self.change_hyperparams)
            st.subheader("Train size")
            st.slider(label=f"Train data %", min_value=1, max_value=99, value=90, key="input_train_size_solo_vae", on_change=self.change_hyperparams)

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
        """
        Init solo model object with given hyperparameters.

        Parameters
        ----------
        device: str
            Device to run the training on (CPU or Cuda if GPU is available).

        epochs: str
            Number of epochs to train models for.

        learning_rate: float
            Adjustable learning rate for improving gradient descent opttimizer.

        train_size: int
            The percentage of the dataset to be used for training (not testing/validation).

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/solo_model.png

        Example
        -------
        """
        self.model_dict = {
            "model": None,
            "n_epochs": st.session_state.ni_vae_epochs,
            "train_size": round((100 - st.session_state.input_train_size_solo_vae) / 100, 2),
            "lr": st.session_state.ni_solo_lr,
            "device": self.device
        }

        st.session_state['model_obj'] = self.model_dict

    def build_model(self):
        lr = st.session_state.ni_solo_lr
        epochs = st.session_state.ni_vae_epochs
        train_size = round((100 - st.session_state.input_train_size_solo_vae) / 100, 2)
        model = solo_model(adata=self.adata, epochs=epochs, lr=lr, train_size=train_size, use_gpu=(self.device == "cuda"))
        st.session_state.model_obj["model"] = model

    def change_hyperparams(self):
        self.init_model()
        self.build_model()
        

def create_citeseq(adata):
    create_model = Create_CiteSeq_model(adata)

    create_model.init_device()

    create_model.draw_page()

    create_model.init_model()

    create_model.build_model(adata.to_df())

    create_model.create_datasets(test_size=st.session_state.model_obj["test_split"])

    st.subheader("Model summary")
    st.json(st.session_state.model_obj, expanded=False)

def create_solo(adata):
    create_model = Create_Solo_model(adata)

    create_model.init_device()

    create_model.draw_page()

    create_model.init_model()

    create_model.build_model()

    st.subheader("Model summary")
    st.json(st.session_state.model_obj, expanded=False)


def change_model():
    adata = st.session_state.adata_state.current.adata
    if st.session_state.sb_model_selection == 'Citeseq (dimensionality reduction)':
        create_citeseq(adata)
    elif st.session_state.sb_model_selection == 'Solo (doublet removal)':
        create_solo(adata)


try:
    sidebar = Sidebar()
    sidebar.show()

    st.title("Create Model")

    col1, _, _, _ = st.columns(4)
    col1.selectbox(label="model", options=([
        "Citeseq (dimensionality reduction)", 
        "Solo (doublet removal)"
        ]), key='sb_model_selection')

    change_model()

    sidebar.show_preview()
    sidebar.export_script()
    sidebar.delete_experiment_btn()
    sidebar.show_version()
    
except Exception as e:
    if(st.session_state == {}):
        StateManager().load_session()
        st.rerun()
    else:
        st.toast(e, icon="❌")