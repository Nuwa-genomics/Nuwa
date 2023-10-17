import streamlit as st
from ml.model import CiteAutoencoder
from ml.dataset import TabularDataset
from sklearn.model_selection import train_test_split
from torch.utils.data import DataLoader
import numpy as np
import pickle
import sys
import torch

st.set_page_config(layout="wide", page_title='Nuwa', page_icon='🧬')

common_style = """
            <style>
            footer {visibility: hidden;}
            .st-emotion-cache-1cypcdb {background: linear-gradient(180deg, rgb(5, 39, 103) 0%, #3a0647 70%); box-shadow: 1px 0 10px -2px #000;}
            .st-emotion-cache-86cver {rgba(250, 250, 250, 0.6)}
            </style>
            """
st.markdown(common_style, unsafe_allow_html=True)

if 'adata' not in st.session_state:
    tmp_file = open("./tmp/adata.pkl",'rb')
    cached_adata = pickle.load(tmp_file)
    st.session_state["adata"] = cached_adata

try:
    adata = st.session_state["adata"]
except KeyError as ke:
    print('Key Not Found in Employee Dictionary:', ke)


class CreateModel:
    def __init__(self, adata):
        st.title("Create Model")

        self.adata = adata
        print("running init")


    def draw_page(self):
        col1, _, _, _ = st.columns(4)

        with col1:
            st.selectbox(label="model", options=(["Citeseq"]))
            st.selectbox(label="Device", options=(self.device_options), on_change=self.set_device, key="select_device")

        col1, col2 = st.columns(2, gap="large")

        with col1:
            self.set_hyper_params()
            self.set_train_test_split()
        with col2:
            self.set_autoencoder()
        
    def init_device(self):
        #init devices if not already exists
        if torch.cuda.is_available():
            self.device_options = ["Cuda (GPU)", "CPU"]
        else:
            self.device_options = ["CPU"]
        
        if not 'select_device' in st.session_state:
            self.device = "cuda" if torch.cuda.is_available() else "cpu"
        else:
            self.device = "cpu" if st.session_state.select_device == "CPU" else "cuda"

        st.session_state.device = self.device

    def set_device(self):
        
        self.device = "cpu" if st.session_state["select_device"] == "CPU" else "cuda"
        st.session_state["device"] = self.device

    def init_model(self):
        #initialize model object
        if 'model_obj' not in st.session_state:
            self.model_dict = {
                "model": None,
                "lr": 1e-2,
                "n_epochs": 100,
                "n_features": adata.to_df().shape[1],
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
    
    def change_hyperparam(self):
        st.session_state.model_obj["lr"] = st.session_state.input_lr
        st.session_state.model_obj["n_epochs"] = st.session_state.input_nepochs
        st.session_state.model_obj["optim"] = st.session_state.input_optim

    def set_hyper_params(self):
        st.subheader("Set model hyperparameters")
        st.number_input(label="Epochs", min_value=1, key="input_nepochs", value=100, on_change=self.change_hyperparam)
        st.number_input(label="Learning rate", min_value=1e-4, max_value=1.0, value=1e-3, key="input_lr", step=1e-3, format='%.4f', on_change=self.change_hyperparam)
        st.selectbox(label="Optimizer", options=(["Adam", "SGD", "RMSProp"]), key="input_optim", on_change=self.change_hyperparam)

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


create_model = CreateModel(adata)

create_model.init_device()

create_model.draw_page()

create_model.init_model()

create_model.build_model(adata.to_df())

create_model.create_datasets(adata.to_df(), test_size=st.session_state.model_obj["test_split"])

st.subheader("Model summary")
st.json(st.session_state.model_obj, expanded=False)