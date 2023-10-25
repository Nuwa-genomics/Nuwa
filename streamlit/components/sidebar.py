import streamlit as st
import scanpy as sc
import pickle

from models.AdataModel import AdataModel

def show_preview():
        with st.sidebar:
            with st.expander(label="Show Preview"):
                st.subheader("Anndata preview")
                with st.container():
                    st.markdown(f"<p style='font-size: 14px; color: rgba(255, 255, 255, 0.75)'>{st.session_state.current_adata}</p>", unsafe_allow_html=True)
                    #st.code(st.session_state.adata)

def get_adata(adataList, name) -> AdataModel:
    for adata in adataList:
        if adata.name == name:
            return adata
    return 0

def show_sidebar(adataList):
    with st.sidebar:
        def add_experiment():
            print("hi")
        def set_adata():
            adata_bytes = get_adata(adataList, name=st.session_state.sb_adata_selection).adata
            st.session_state["current_adata"] = pickle.loads(adata_bytes)
        def save_file():
            selected_adata = get_adata(adataList, name=st.session_state.sb_adata_selection)
            if not selected_adata:
                st.toast("Couldn't find selected adata to save")
            else:
                sc.write(filename=f"downloads/{selected_adata.name}.h5ad", adata=pickle.loads(selected_adata.adata))
                st.toast("Downloaded file", icon='✅')
        options=[item.name for item in adataList]
        st.selectbox(label="Current Experiment:", options=options, key="sb_adata_selection", on_change=set_adata)
        st.button(label="Download adata file", on_click=save_file, use_container_width=True, key="btn_save_adata")
        st.button(label="Add experiment", on_click=add_experiment, use_container_width=True, key="btn_add_adata")