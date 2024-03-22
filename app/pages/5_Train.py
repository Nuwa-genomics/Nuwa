from enum import Enum
import streamlit as st
import json
from streamlit_lottie import st_lottie
from ml.citeseq.model import CiteAutoencoder
from lightning.pytorch.callbacks import Callback
from ml.citeseq.train import train_model
from ml.citeseq.train import train_model
from ml.solo_scvi.solo_model import solo_model
from ml.linear_vae.linear_vae import LDVAE
from ml.DeepST.deepst.main import *
from models.AdataModel import AdataModel
from components.sidebar import *
from state.AdataState import AdataState
import os
from state.StateManager import StateManager

st.set_page_config(page_title='Nuwa', page_icon='🧬', layout="centered")

os.chdir('/app')

with open('css/common.css') as f:
    common_style = f"""
                <style>
                {f.read()}
                </style>
                """
    st.markdown(common_style, unsafe_allow_html=True)


if 'train_loss' not in st.session_state:
    st.session_state["train_loss"] = 0.00


class Devices(Enum):
    CPU = "CPU"
    CUDA = "Cuda (GPU)"


class ProgressCallback(Callback):

    def __init__(self, pg_prefix="", n_epochs=None, play_complete_animation=True):
        self.pg_prefix = pg_prefix
        self.play_complete_animation = play_complete_animation
        if n_epochs == None:
            self.n_epochs = st.session_state["model_obj"]['n_epochs']
        else:
            self.n_epochs = n_epochs

    def on_train_start(self, trainer, pl_module):
        print("Training is starting")
      

    def on_train_end(self, trainer, pl_module):
        train.pb.progress(100, text=f"{self.pg_prefix}Epoch {trainer.current_epoch}/{self.n_epochs}")
        if self.play_complete_animation:
            with train.placeholder.container():
                st.markdown("<h3 style='text-align: center; margin-bottom: 2rem'>Training Complete</h3>", unsafe_allow_html=True)
                with open("animations/tick.json") as source:
                        
                    complete_animation = json.load(source)
                    st_lottie(complete_animation, width=700, height=200, loop=False, speed=1.2)

                    train_loss = 600
                    valid_loss = 600

                    subcol1, subcol2, subcol3, subcol4, subcol5 = st.columns(5, gap="medium")
                        
                    with subcol2:
                        st.metric(label="Train loss", value=round(train_loss, 4))
                    with subcol4:
                        st.metric(label="Validation loss", value=round(valid_loss, 4))

    def on_train_epoch_end(self, trainer, pl_module):
        current_epoch = trainer.current_epoch
        train.pb.progress(round((current_epoch / self.n_epochs) * 100), text=f"{self.pg_prefix}Epoch {current_epoch}/{self.n_epochs}")

        print(trainer.callback_metrics)
        if "elbo_validation" in trainer.callback_metrics and "elbo_train" in trainer.callback_metrics:
            print(trainer.callback_metrics["elbo_validation"].item())
            print(trainer.callback_metrics["elbo_train"].item())

        #print(trainer.callback_metrics)


class Train:
    def __init__(self, adata):
        self.model = st.session_state["model_obj"]

        self.placeholder = st.empty()

        self.pb = st.progress(value=0, text=f"Epoch 0/{self.model['n_epochs']}")

        if "device" not in st.session_state:
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
        n_epochs = self.model['n_epochs']
        self.pb.progress(round((epoch/n_epochs)*100), text=f"Epoch {epoch}/{n_epochs}")
        if epoch == self.model['n_epochs']:
            with self.placeholder.container():
                st.markdown("<h3 style='text-align: center; margin-bottom: 2rem'>Training Complete</h3>", unsafe_allow_html=True)
                with open("animations/tick.json") as source:
                    
                    complete_animation = json.load(source)
                    st_lottie(complete_animation, width=700, height=200, loop=False, speed=1.2)

                    subcol1, subcol2, subcol3, subcol4, subcol5 = st.columns(5, gap="medium")
                    
                    with subcol2:
                        st.metric(label="Train loss", value=round(train_loss, 4))
                    with subcol4:
                        st.metric(label="Validation loss", value=round(valid_loss, 4))
                    
                    

    def train_pgb_non_specific(self, percent, text):
        n_epochs = self.model['n_epochs']
        self.pb.progress(percent, text=text)
        if percent >= 100:
            with self.placeholder.container():
                st.markdown("<h3 style='text-align: center; margin-bottom: 2rem'>Training Complete</h3>", unsafe_allow_html=True)
                with open("animations/tick.json") as source:
                    complete_animation = json.load(source)
                    st_lottie(complete_animation, width=700, height=200, loop=False, speed=1.2)

                    subcol1, subcol2 = st.columns(2, gap="medium")
                    
                    #with subcol1:
                        #st.metric(label="Train loss", value=round(train_loss, 4))
                    #with subcol2:
                        #st.metric(label="Validation loss", value=round(valid_loss, 4))


    def train(self):
        if(isinstance(self.model['model'], CiteAutoencoder)):
            trained_model, losses = train_model(
            model=self.model["model"], 
            train_dl=self.model["train_dl"], 
            valid_dl=self.model["valid_dl"], 
            lr=self.model["lr"], 
            epochs=self.model["n_epochs"], 
            callback_on_epoch=self.train_pgb, 
            verbose=True,
            device=self.device
            )

            st.session_state["losses"] = losses #save losses
            st.session_state["trained_model"] = trained_model #save model to local session
       
        elif(isinstance(self.model['model'], solo_model)):
            self.model['model'].train_vae(callbacks=[ProgressCallback(pg_prefix="Vae model: ", play_complete_animation=False)])
            self.model['model'].train_solo(callbacks=[ProgressCallback(pg_prefix="Solo model: ")])
            self.model['model'].predict_solo()
            st.session_state["trained_model"] = self.model['model']

        elif(isinstance(self.model['model'], LDVAE)):
            st.session_state["trained_model"] = self.model['model'].train(callbacks=[ProgressCallback()])

try:
    adata = st.session_state.adata_state.current.adata

    sidebar = Sidebar()

    sidebar.show()
    
    train = Train(adata)

    train.draw_animation()

    train.train()

    
except Exception as e:
    if(st.session_state == {}):
        StateManager().load_session()
        st.rerun()
    else:
        st.toast(e, icon="❌")




