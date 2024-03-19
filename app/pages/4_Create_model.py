import os
os.chdir('/app')

from anndata import AnnData
import streamlit as st
from ml.citeseq.model import CiteAutoencoder
from ml.citeseq.dataset import TabularDataset
from sklearn.model_selection import train_test_split
from torch.utils.data import DataLoader
import numpy as np
from ml.solo_scvi.solo_model import *
from ml.DeepST.deepst.main import *
from ml.linear_vae.linear_vae import *
import torch
from components.sidebar import *
from state.StateManager import StateManager
from enum import Enum


st.set_page_config(layout="wide", page_title='Nuwa', page_icon='🧬')

with open('css/common.css') as f:
    common_style = f"""
                <style>
                {f.read()}
                </style>
                """
    st.markdown(common_style, unsafe_allow_html=True)


class NNModels(Enum):
    CITESEQ = "Citeseq (dimensionality reduction)"
    SOLO = "Solo (doublet removal)"
    LDVAE = "LDVAE (dimensionality reduction)"

class Devices(Enum):
    CPU = "CPU"
    CUDA = "Cuda (GPU)"


class CreateModel:

    def __init__(self, adata: AnnData):
        self.adata = adata
        self.device = "cuda" if st.session_state["create_model:sb_device"] == Devices.CUDA.value else "cpu"
        st.session_state["device"] = self.device

    def draw_page(self):
        pass

    def init_model(self):
        pass

    def build_model(self):
        pass

    def summary(self):
        st.subheader("Model summary")
        st.json(st.session_state.model_obj, expanded=False)


class Create_CiteSeq_model(CreateModel):
    """
    Train a deep autoencoder model for dimensionality reduction on gene expression data. This is an alternative to other DR algorithms such as tSNE. Original Github repo can be found [here](https://github.com/naity/citeseq_autoencoder)
    
    Notes
    -----
    .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/create_model_citeseq_page.png
    """
    def __init__(self, adata):
        super(self.__class__, self).__init__(adata)
        

    def draw_page(self):
        col1, col2 = st.columns(2, gap="large")

        with col1:
            self.set_hyperparams()
            self.set_train_test_split()
        with col2:
            self.set_autoencoder()
        

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

    def build_model(self):

        df = self.adata.to_df()
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

        self.create_datasets(test_size=st.session_state.model_obj["test_split"])

class Create_Solo_model(CreateModel):
    """
    Create a Solo model for doublet detection and removal. A package in [scvi tools](https://docs.scvi-tools.org/en/stable/user_guide/models/solo.html)

    Notes
    -----
    .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/create_model_solo_page.png
    """
    def __init__(self, adata):
        super(self.__class__, self).__init__(adata)

    def draw_page(self):
        col1, col2 = st.columns(2)
        with col1:
            st.subheader("Model parameters")
            st.number_input(label="Epochs", min_value=1, value=400, key="ni_vae_epochs", on_change=self.change_hyperparams)
            st.number_input(label="Learning rate", min_value=1e-4, max_value=10.0, value=1e-3, format='%.4f', key="ni_solo_lr", on_change=self.change_hyperparams)
            st.subheader("Train size")
            st.slider(label=f"Train data %", min_value=1, max_value=99, value=90, key="input_train_size_solo_vae", on_change=self.change_hyperparams)
            

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


class Create_linear_vae(CreateModel):
    """Create a Solo model for doublet detection and removal. A package in [scvi tools](https://docs.scvi-tools.org/en/1.0.2/tutorials/notebooks/linear_decoder.html)

    Notes
    -----
    .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/create_model_solo_page.png
    """

    def __init__(self, adata):
        super(self.__class__, self).__init__(adata)

    def draw_page(self):
        col1, col2, _, _ = st.columns(4)
        col1.number_input(label="max_epochs", value=250, step=1, min_value=1, format="%i", key="ni:create_model:ldvae:max_epochs")
        col2.number_input(label="lr", value=5e-3, step=1e-3, format='%.4f', key="ni:create_model:ldvae:lr")

    def init_model(self):
        """
        Init Linearly decoded VAE model object with given hyperparameters.

        Parameters
        ----------
        device: str
            Device to run the training on (CPU or Cuda if GPU is available).

        epochs: str
            Number of epochs to train models for.

        learning_rate: float
            Adjustable learning rate for improving gradient descent opttimizer.


        Notes
        -----
        .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/ldvae_model.png

        Example
        -------
        """
        self.model_dict = {
            "model": None,
            "n_epochs": st.session_state["ni:create_model:ldvae:max_epochs"],
            "lr": st.session_state["ni:create_model:ldvae:lr"],
            "device": self.device
        }

        st.session_state['model_obj'] = self.model_dict

    def build_model(self):
        max_epochs = st.session_state["ni:create_model:ldvae:max_epochs"]
        lr = st.session_state["ni:create_model:ldvae:lr"]
        model = LDVAE(adata=self.adata, max_epochs=max_epochs, lr=lr, use_gpu=(self.device == "cuda"))
        st.session_state.model_obj["model"] = model
        

def create_model(adata):
    if st.session_state["create_model:sb_model"] == NNModels.CITESEQ.value:
        model = Create_CiteSeq_model(adata)
    elif st.session_state["create_model:sb_model"] == NNModels.SOLO.value:
        model = Create_Solo_model(adata)
    elif st.session_state["create_model:sb_model"] == NNModels.LDVAE.value:
        model = Create_linear_vae(adata)

    model.draw_page()

    model.init_model()

    model.build_model()

    model.summary()


try:
    sidebar = Sidebar()
    sidebar.show()

    st.title("Create Model")

    col1, col2, _, _ = st.columns(4)
    col1.selectbox(label="model", options=([NNModels.CITESEQ.value, NNModels.SOLO.value, NNModels.LDVAE.value]), key='create_model:sb_model')
    col2.selectbox(label="Device", options=([Devices.CUDA.value, Devices.CPU.value] if torch.cuda.is_available() else [Devices.CPU.value]), key="create_model:sb_device")

    st.divider()

    adata = st.session_state.adata_state.current.adata.copy()

    create_model(adata)

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