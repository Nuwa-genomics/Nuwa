import streamlit as st
import time
import json
import streamlit.components.v1 as com
from streamlit_lottie import st_lottie
from ml.train import train_model
import pickle
import matplotlib.pyplot as plt
import numpy as np
from ml.train import train_model

common_style = """
            <style>
            footer {visibility: hidden;}
            .st-emotion-cache-1cypcdb {background: linear-gradient(180deg, rgb(5, 39, 103) 0%, #3a0647 70%); box-shadow: 1px 0 10px -2px #000;}
            .st-emotion-cache-86cver {rgba(250, 250, 250, 0.6)}
            .st-e2 {background-color: #004dcf;}
            </style>
            """
st.markdown(common_style, unsafe_allow_html=True)


if 'adata' not in st.session_state:
    tmp_file = open("./tmp/adata.pkl",'rb')
    cached_adata = pickle.load(tmp_file)
    st.session_state["adata"] = cached_adata

if 'train_loss' not in st.session_state:
    st.session_state["train_loss"] = 0.00


class Train:
    def __init__(self, adata):
        self.model_params = st.session_state["model_obj"]

        self.placeholder = st.empty()

        self.pb = st.progress(value=0, text=f"Epoch 0/{self.model_params['n_epochs']}")

        if 'device' not in st.session_state:
            self.device = "cpu"
        else:
            self.device = st.session_state.device

        self.animation_fn = "animations/cpu.json" if self.device == "cpu" else "animations/gpu.json"

    def draw_animation(self):
        with open(self.animation_fn) as source:
            animation = json.load(source)

        with self.placeholder.container():
            st.markdown("<h3 style='text-align: center; margin-bottom: 2rem'>Training Model</h3>", unsafe_allow_html=True)
            st_lottie(animation, width=700, height=200)


    def train_pgb(self, epoch, train_loss, valid_loss):
        n_epochs = self.model_params['n_epochs']
        self.pb.progress(round((epoch/n_epochs)*100), text=f"Epoch {epoch}/{n_epochs}")
        if epoch == self.model_params['n_epochs']:
            with self.placeholder.container():
                st.markdown("<h3 style='text-align: center; margin-bottom: 2rem'>Training Complete</h3>", unsafe_allow_html=True)
                with open("animations/tick.json") as source:
                    complete_animation = json.load(source)
                    st_lottie(complete_animation, width=700, height=200, loop=False, speed=1.2)

                    subcol1, subcol2 = st.columns(2, gap="medium")
                    
                    with subcol1:
                        st.metric(label="Train loss", value=round(train_loss, 4))
                    with subcol2:
                        st.metric(label="Validation loss", value=round(valid_loss, 4))


    def train(self):
        trained_model, losses = train_model(
        model=self.model_params["model"], 
        train_dl=self.model_params["train_dl"], 
        valid_dl=self.model_params["valid_dl"], 
        lr=self.model_params["lr"], 
        epochs=self.model_params["n_epochs"], 
        callback_on_epoch=self.train_pgb, 
        verbose=True,
        device=self.device
        )

        st.session_state["losses"] = losses #save losses
        st.session_state["trained_model"] = trained_model #save model to local session


adata = st.session_state["adata"]

train = Train(adata)

train.draw_animation()

train.train()

