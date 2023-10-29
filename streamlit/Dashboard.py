import streamlit as st
import pandas as pd
from models.WorkspaceModel import WorkspaceModel
from datetime import datetime
from SQL.Workspace import *
from database.database import *
from database.schemas import schemas
from sqlalchemy.orm import Session
from hashlib import sha256

#create databases if not already present
Base.metadata.create_all(engine)

st.set_page_config(page_title='Nuwa', page_icon='🧬')

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
            path = f"/streamlit-volume/{dir}"
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
        except Exception as e:
            st.error(e)
        finally:
            st.toast("Successfully created workspace", icon='✅')

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
                                work_id = i[0][-1]

                                for workspace in st.session_state.workspaces:
                                    if workspace.id == int(work_id):
                                        st.session_state["current_workspace"] = workspace
                                        with st.sidebar:
                                            st.subheader(f"{workspace.workspace_name}")
                                            st.write(f"{workspace.description}")
                                            st.write(f"{workspace.created.ctime()}")
                                            st.divider()
                                            
                                
                            

        with st.sidebar:
            st.button(label="New workspace", on_click=self.new_workspace, key="btn_new_workspace", use_container_width=True)

dashboard = Dashboard()



