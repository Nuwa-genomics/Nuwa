import streamlit as st
import json
from streamlit_lottie import st_lottie
from ml.citeseq.model import CiteAutoencoder
from ml.citeseq.train import train_model
import pickle
from ml.citeseq.train import train_model
from ml.solo_scvi.solo_model import solo_model
from ml.DeepST.deepst.main import *
from models.AdataModel import AdataModel
from components.sidebar import *

st.set_page_config(page_title='Nuwa', page_icon='🧬', layout="centered")

with open('css/common.css') as f:
    common_style = f"""
                <style>
                {f.read()}
                </style>
                """
    st.markdown(common_style, unsafe_allow_html=True)


if 'train_loss' not in st.session_state:
    st.session_state["train_loss"] = 0.00


class Train:
    def __init__(self, adata):
        self.model = st.session_state["model_obj"]

        self.placeholder = st.empty()

        self.pb = st.progress(value=0, text=f"Epoch 0/{self.model['n_epochs']}")

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
        n_epochs = self.model['n_epochs']
        self.pb.progress(round((epoch/n_epochs)*100), text=f"Epoch {epoch}/{n_epochs}")
        if epoch == self.model['n_epochs']:
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
            st.session_state["trained_model"] = self.model['model'].train(callback=self.train_pgb_non_specific)
        elif(isinstance(self.model['model'], DeepSTModel)):
            st.session_state["trained_model"] = self.model['model'].run(callback=self.train_pgb_non_specific)

try:
    adata_model: AdataModel = st.session_state["adata"]
    show_sidebar(adata_model)

    adata = get_adata(adataList=adata_model, name=st.session_state.sb_adata_selection).adata
    st.session_state["current_adata"] = adata

    show_preview()

    train = Train(adata)

    train.draw_animation()

    train.train()

except KeyError as ke:
    print("KeyError: ", ke)
    st.error("Couldn't find adata object in session, have you uploaded one?")
    
except Exception as e:
    st.error(e)



