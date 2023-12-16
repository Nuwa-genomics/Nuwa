import streamlit as st
import pandas as pd
from models.WorkspaceModel import WorkspaceModel
from models.AdataModel import AdataModel
from utils.AdataState import AdataState
from datetime import datetime
from database.database import *
from database.schemas import schemas
from sqlalchemy.orm import Session
from hashlib import sha256
import time
import os
import scanpy as sc

#create databases if not already present
Base.metadata.create_all(engine)

st.set_page_config(page_title='Nuwa', page_icon='🧬')

os.chdir('/app')

common_style = """
            <style>
            footer {visibility: hidden;}
            .st-emotion-cache-1cypcdb {background: linear-gradient(180deg, rgb(5, 39, 103) 0%, #3a0647 70%); box-shadow: 1px 0 10px -2px #000;}
            .st-emotion-cache-86cver {rgba(250, 250, 250, 0.6)}
            </style>
            """
st.markdown(common_style, unsafe_allow_html=True)

with open('css/workspace.css') as f:
    st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)

class Dashboard:
    def __init__(self):
        self.conn: Session = SessionLocal()
        self.get_workspaces()
        self.draw_page()
    
    def get_workspaces(self):
        st.session_state['workspaces'] = self.conn.query(schemas.Workspaces).all()


    def write_workspace_to_db(self):
        try:
            #create workspace dir
            dir = sha256(st.session_state.ti_new_workspace_name.encode('utf-8')).hexdigest()[:16]
            dir = f"{st.session_state.ti_new_workspace_name}_{dir}"
            path = f"/streamlit-volume/{dir}/"
            if not os.path.exists(path):
                os.mkdir(path)
            os.environ['WORKDIR'] = path # set dir
            new_workspace = schemas.Workspaces(
                workspace_name = st.session_state.ti_new_workspace_name,
                description = st.session_state.ti_new_workspace_desc,
                data_dir = path
            )
            self.conn.add(new_workspace)
            self.conn.commit()
            self.conn.refresh(new_workspace)
            st.toast("Successfully created workspace", icon='✅')
        except Exception as e:
            st.error(e)
            st.toast("Failed to create workspace", icon="❌")
            

    def new_workspace(self):
        with st.sidebar:
            with st.form(key="new_workspace_form"):
                st.subheader("Create New Workspace")
                st.text_input(label="Name", key="ti_new_workspace_name")
                st.text_input(label="Description", key="ti_new_workspace_desc")
                st.form_submit_button(label="Save", on_click=self.write_workspace_to_db)


    def draw_page(self):
        st.title("Workspaces")
    
        col1, col2, col3, col4 = st.columns(4, gap="large")
        columns = [col1, col2, col3, col4]

        if len(st.session_state.workspaces) > 0:
            for i, workspace in enumerate(st.session_state.workspaces): 
                with columns[i % 4]:
                    button = st.button(label=workspace.workspace_name, key=f"btn_workspace_{workspace.id}")

                    if button:
                        for i in st.session_state.items():
                            if i[1] and (i[0].__contains__("btn_workspace")):
                                work_id = i[0].split(sep='_')[-1]

                                for workspace in st.session_state.workspaces:
                                    if workspace.id == int(work_id):
                                        st.session_state["current_workspace"] = workspace
                                        #load adata into adata state if exists
                                        adata: schemas.Adata = self.conn.query(schemas.Adata).filter(schemas.Adata.work_id == workspace.id).order_by(schemas.Adata.created).first()
                                        if adata:
                                            anndata = sc.read_h5ad(filename=adata.filename)
                                            active_adata: AdataModel = AdataModel(work_id=adata.work_id, adata=anndata, filename=adata.filename, created=adata.created, adata_name=adata.adata_name, notes = adata.notes, id=adata.id)
                                            st.session_state["adata_state"] = AdataState(active=active_adata, insert_into_db=False)
                                        os.environ['WORKDIR'] = workspace.data_dir #set wd
                                        with st.sidebar:
                                            st.subheader(f"{workspace.workspace_name}")
                                            st.markdown(f"""<p style='font-size: 16px; color: rgba(255, 255, 255, 1)'>{workspace.description}<p>""", unsafe_allow_html=True)
                                            st.markdown(f"""<p style='font-size: 16px; color: rgba(255, 255, 255, 0.4)'>{workspace.created.ctime()}<p>""", unsafe_allow_html=True)
                                            def delete_workspace():
                                                try:
                                                    self.conn.query(schemas.Workspaces).filter(schemas.Workspaces.id == st.session_state.current_workspace.id).delete()
                                                    self.conn.commit()
                                                    st.toast("Successfully deleted workspace", icon="✅")
                                                except Exception as e:
                                                    st.toast("Failed to delete workspace", icon="❌")
                                            st.button(label="🗑️ Delete workspace", on_click=delete_workspace, key="btn_delete_workspace")
                                            st.divider()
        else:
            st.markdown("""<p style='font-size: 20px;'>You don't have any workspaces yet 🧪</p>""", unsafe_allow_html=True)


        with st.sidebar:
            st.button(label="New workspace", on_click=self.new_workspace, key="btn_new_workspace", use_container_width=True)

            st.markdown(f"""<div style='position: fixed; margin-left: 5px; bottom: 5px;'>
                        <div style='font-size: 16px; color: rgba(255, 255, 255, 0.4)'>Nuwa v{os.getenv('NUWA_VERSION')}</div>
                        </div>""", unsafe_allow_html=True)

dashboard = Dashboard()



