import os
import streamlit as st
import scanpy as sc
from database.schemas import schemas

from models.AdataModel import AdataModel
from database.database import SessionLocal
from sqlalchemy.orm import Session
from utils.AdataState import AdataState

class Sidebar:
    def __init__(self):
        self.conn: Session = SessionLocal()

    def show_preview(self):
            with st.sidebar:
                with st.expander(label="Show Preview"):
                    st.subheader("Anndata preview")
                    st.markdown(f"""<p style='font-size: 14px; color: rgba(255, 255, 255, 0.75)'>{st.session_state.current_adata}</p>""", unsafe_allow_html=True)

    def delete_experiment_btn(self):
        with st.sidebar:
            delete_btn = st.button(label="🗑️ Delete Experiment", use_container_width=True)
            if delete_btn:
                st.session_state.adata_state.delete_record()



    def show_notes(self):
        with st.sidebar:
            notes = st.session_state.adata_state.current.notes
            display_notes = notes if notes != None else ""
            notes_ta = st.text_area(label="Notes", placeholder="Notes", value=display_notes)
            try:
                st.session_state.adata_state.current.notes = notes_ta
                st.session_state.adata_state.update_record()
            except Exception as e:
                st.toast("Notes failed to save", icon="❌")
                print("Error: ", e)

    def get_adata(self, adataList, name) -> AdataModel:
        for adata in adataList:
            if adata.adata_name == name:
                return adata
        return 0


    def write_adata(self):
        try:
            name = st.session_state.ti_new_adata_name

            st.session_state.adata_state.add_adata(AdataModel(
                work_id=st.session_state.current_workspace.id,
                adata_name=name,
                filename=os.path.join(os.getenv('WORKDIR'), "adata", f"{name}.h5ad"),
                notes=st.session_state.adata_state.current.notes,
            ))
            
        except Exception as e:
            st.error(e)
            print("Error: ", e)


    def add_experiment(self):
        with st.sidebar:
            try:
                with st.form(key="new_workspace_form"):
                    st.subheader("Create New Adata")
                    st.text_input(label="Name", key="ti_new_adata_name")
                    st.text_input(label="Notes", key="ti_new_adata_notes")
                    st.form_submit_button(label="Save", on_click=self.write_adata)
            except Exception as e:
                print("Error: ", e)
                st.error(e)

    def show(self):
        with st.sidebar:
            def set_adata():
                st.session_state.adata_state.switch_adata(st.session_state.sb_adata_selection)
            def save_file():
                selected_adata = st.session_state.adata_state.get_adata(adata_name=st.session_state.sb_adata_selection)
                if not selected_adata:
                    st.toast("Couldn't find selected adata to save")
                else:
                    sc.write(filename=f"{os.getenv('WORKDIR')}/downloads/{selected_adata.adata_name}.h5ad", adata=st.session_state.current_adata)
                    st.toast("Downloaded file", icon='✅')
            options=[item.adata_name for item in st.session_state.adata_state.load_adata(st.session_state.current_workspace.id)]
            st.selectbox(label="Current Experiment:", options=options, key="sb_adata_selection", on_change=set_adata)
            st.button(label="Download adata file", on_click=save_file, use_container_width=True, key="btn_save_adata")
            st.button(label="Add experiment", on_click=self.add_experiment, use_container_width=True, key="btn_add_adata")
            self.show_notes()
