import os
import streamlit as st
import scanpy as sc
from database.schemas import schemas

from models.AdataModel import AdataModel
from database.database import SessionLocal

sidebar_conn: SessionLocal = SessionLocal()

def show_preview():
        with st.sidebar:
            with st.expander(label="Show Preview"):
                st.subheader("Anndata preview")
                with st.container():
                    st.markdown(f"<p style='font-size: 14px; color: rgba(255, 255, 255, 0.75)'>{st.session_state.current_adata}</p>", unsafe_allow_html=True)
                    #st.code(st.session_state.adata)

def get_adata(adataList, name) -> AdataModel:
    for adata in adataList:
        if adata.adata_name == name:
            return adata
    return 0


def write_adata_to_db():
    name = st.session_state.ti_new_adata_name
    new_adata = schemas.Adata(
        work_id = st.session_state.current_workspace.id,
        adata_name = name,
        filename = os.path.join(os.getenv('WORKDIR'), name),
        notes = st.session_state.ti_new_adata_notes
    )
    sidebar_conn.add(new_adata)
    sidebar_conn.commit()
    sidebar_conn.refresh(new_adata)
    st.toast("Created new adata", icon="✅")


def add_experiment():
    with st.sidebar:
        try:
            with st.form(key="new_workspace_form"):
                st.subheader("Create New Adata")
                st.text_input(label="Name", key="ti_new_adata_name")
                st.text_input(label="Notes", key="ti_new_adata_notes")
                st.form_submit_button(label="Save", on_click=write_adata_to_db)
        except Exception as e:
            print("Error: ", e)
            st.error(e)

def show_sidebar(adataList):
    with st.sidebar:
        def set_adata():
            adata = get_adata(adataList, name=st.session_state.sb_adata_selection).adata
            st.session_state["current_adata"] = adata
        def save_file():
            selected_adata = get_adata(adataList, name=st.session_state.sb_adata_selection)
            if not selected_adata:
                st.toast("Couldn't find selected adata to save")
            else:
                sc.write(filename=f"{os.getenv('WORKDIR')}/downloads/{selected_adata.adata_name}.h5ad", adata=st.session_state.current_adata)
                st.toast("Downloaded file", icon='✅')
        options=[item.adata_name for item in adataList]
        st.selectbox(label="Current Experiment:", options=reversed(options), key="sb_adata_selection", on_change=set_adata)
        st.button(label="Download adata file", on_click=save_file, use_container_width=True, key="btn_save_adata")
        st.button(label="Add experiment", on_click=add_experiment, use_container_width=True, key="btn_add_adata")