import hashlib
import streamlit as st
import pandas as pd
from models.WorkspaceModel import WorkspaceModel
from models.AdataModel import AdataModel
from state.AdataState import AdataState
from utils.file_utils import *
from datetime import datetime
from database.database import *
from database.schemas import schemas
from sqlalchemy.orm import Session
from pathlib import Path
import shutil
import zipfile
import hashlib as hash
from hashlib import sha256
import time
import os
import scanpy as sc
from database.db_export import *
import pickle

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
    """
    Dashboard page for all workspaces. Create, export, import and delete workspaces here.

    Notes
    -----
    .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/dashboard_page.png
    """
    def __init__(self):
        self.conn: Session = SessionLocal()
        if not os.path.exists('/streamlit-volume/exported_workspaces'):
            os.mkdir('/streamlit-volume/exported_workspaces') 
        self.set_workspaces()
        self.draw_page()
    
    def set_workspaces(self):
        workspaces = self.conn.query(schemas.Workspaces).all()
        st.session_state['workspaces'] = [WorkspaceModel(id=workspace.id, workspace_name=workspace.workspace_name, data_dir=workspace.data_dir, created=workspace.created, description=workspace.description) for workspace in workspaces]


    def _write_workspace_to_db(self):
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


    def _save_imported_files(self, file_name):
        new_dir = f'/streamlit-volume/{file_name[:-4]}'
        if not os.path.isdir(new_dir):
            os.mkdir(new_dir)
        with zipfile.ZipFile(os.path.join('/tmp', file_name), mode="r") as zip_f:
            zip_f.extractall(new_dir)
        with zipfile.ZipFile(os.path.join(new_dir, file_name), mode="r") as zip_f:
            zip_f.extractall(new_dir)
        os.remove(os.path.join(new_dir, file_name)) # remove zip file
        os.remove(os.path.join('/tmp', file_name)) # remove in tmp dir

        # read database file
        if import_db_from_json(json_filepath=f"{new_dir}/db.json") == 0: # success code for importing to db
            st.toast("Successfully imported workspace", icon="✅")

    def new_workspace(self):
        """
        Create a new workspace for analysis.

        Parameters
        ----------
        name: str
            Name of workspace.

        description: str
            A short description of the workspace, its purpose, type of analysis etc.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/new_workspace.png
        """
        with st.sidebar:
            with st.form(key="new_workspace_form"):
                st.subheader("Create New Workspace")
                st.text_input(label="Name", key="ti_new_workspace_name")
                st.text_input(label="Description", key="ti_new_workspace_desc")
                col1, _, _ = st.columns(3)
                save_btn = col1.form_submit_button(label="Save", use_container_width=True)
                if save_btn:
                    self._write_workspace_to_db()
                    self.set_workspaces()


    def import_workspace(self):
        """
        Import a workspace from a zip file exported from another project. Optionally use a checksum to verify the files have not been modified.

        Parameters
        ----------
        workspace: List[UploadedFile]
            Zip file containing workspace files.

        checksum: str
            Optional sha256sum of zip file.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/import_workspace.png
        """
        with st.sidebar:
            with st.form(key="form_import_workspace"):
                st.subheader("Import workspace")
                f_uploader = st.file_uploader(label="Workspace zip file", type=['zip', 'sha256'], accept_multiple_files=True)
                checksum = st.text_input(label="Checksum", placeholder="optional", max_chars=64, help="A checksum helps to see if the original files have been corrupted or modified.")

                col1, _, _ = st.columns(3)
                save_btn = col1.form_submit_button(label="Save", use_container_width=True)
                if save_btn:
                    if len(f_uploader) > 1:
                        st.toast("Too many files", icon="❌")
                    else:
                        for f in f_uploader:
                            ext = f.name.split('.')[-1]
                            if ext == "zip":
                                # calculate checksum
                                BLOCKSIZE = 65536
                                sha256 = hash.sha256()
                                file_buffer = f.read(BLOCKSIZE)
                                while len(file_buffer) > 0: 
                                    sha256.update(file_buffer)
                                    file_buffer = f.read(BLOCKSIZE)
                                digest = sha256.hexdigest()
                                st.session_state["sha256_to_verify"] = digest
                    
                                zip_filename = f.name
      
                                with zipfile.ZipFile(os.path.join('/tmp', f.name), mode="w") as zip_f:
                                    zip_f.writestr(f.name, data=f.getvalue())
                                
                                
                        if checksum: 
                            if checksum == st.session_state["sha256_to_verify"]: # checksums match
                                st.success("Checksums match")
                                self._save_imported_files(zip_filename)
                            else:
                                st.toast("Checksums don't match", icon="❌")
                            
                            #delete from session state
                            del st.session_state["sha256_true"]
                        else:
                            self._save_imported_files(zip_filename)


    def export_workspace(self):
        """
        Export a workspace which can be shared and re-imported into another project. Includes a checksum.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/export_workspace.png
        """
        with st.sidebar:
                                            
            with st.spinner(text="Zipping archive"):

                output_file = f"/streamlit-volume/exported_workspaces/{st.session_state.current_workspace.workspace_name}"

                if not os.path.exists(output_file):
                    os.mkdir(output_file)

                export_db_as_json(workspace_id=st.session_state.current_workspace.id, 
                    output_file=os.path.join(st.session_state.current_workspace.data_dir, 'db.json'))

                workspace_size = get_dir_size(st.session_state.current_workspace.data_dir) / 1000000
                                                        
                if workspace_size > 1000: #1GB
                    st.toast(f"Workspace contents are large ({round(workspace_size / 1000, ndigits=2)} GB) may take longer", icon="⚠️")
                                                        
                shutil.make_archive(base_name=f"{output_file}/{st.session_state.current_workspace.workspace_name}", 
                    format='zip', root_dir=st.session_state.current_workspace.data_dir)                            
                                                    
                                                    
            BLOCKSIZE = 65536

            with st.spinner(text="Computing hash"):

                sha256 = hash.sha256()
                with open(f'{output_file}/{st.session_state.current_workspace.workspace_name}.zip', 'rb') as f:
                    file_buffer = f.read(BLOCKSIZE)
                    while len(file_buffer) > 0:
                        sha256.update(file_buffer)
                        file_buffer = f.read(BLOCKSIZE)

                    digest = sha256.hexdigest()
                    with open(f'{output_file}/checksum.sha256', 'wb') as sha256_f:
                        sha256_f.write(digest.encode('utf-8'))

                    st.markdown("""<p style='font-size:18px; text-align: center; font-weight: bold;'>Exported files</p>""", unsafe_allow_html=True)
                    st.download_button(label="📁 Download", data=f, file_name=st.session_state.current_workspace.workspace_name, mime="application/zip", use_container_width=True)
                    st.info(f"File saved to:\n{output_file}")

                st.markdown("### 🔑 Sha256 sum:", help="The checksum for verifying if files have been modified. The recipient should paste this string when importing a workspace. Checksums are stored in /streamlit-volume/exported_workspaces.")
                st.code(digest)
                # write hash to txt file

            st.toast("Completed workspace export!", icon="✅")
            st.divider()


    def delete_workspace(self):
        """
        Delete a workspace from the dashboard. This will not delete the files associated with the workspace and can be found in /streamlit-volume/<projectname_id>.
        
        Notes
        -----
        .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/delete_workspace.png
        """
        try:
            self.conn.query(schemas.Workspaces).filter(schemas.Workspaces.id == st.session_state.current_workspace.id).delete()
            self.conn.commit()
            st.toast("Successfully deleted workspace", icon="✅")
        except Exception as e:
            st.toast("Failed to delete workspace", icon="❌")
                                                


    def draw_page(self):
        st.title("Workspaces")
    
        subheader_container = st.container()
        col1, col2, col3, col4 = st.columns(4, gap="large")
        columns = [col1, col2, col3, col4]

        if len(st.session_state.workspaces) > 0:
            subheader_container.subheader("Select a workspace")
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
                                        with open(os.path.join(os.getenv('TMP_DIR'), 'session_state.pkl'), 'wb') as pkl_file: 
                                            pickle.dump(st.session_state, pkl_file)
                                        with st.sidebar:
                                            st.subheader(f"{workspace.workspace_name}")
                                            st.markdown(f"""<p style='font-size: 16px; color: rgba(255, 255, 255, 1)'>{workspace.description}<p>""", unsafe_allow_html=True)
                                            st.markdown(f"""<p style='font-size: 16px; color: rgba(255, 255, 255, 0.4)'>{workspace.created.ctime()}<p>""", unsafe_allow_html=True)

                                            sidebar_col1, sidebar_col2 = st.columns(2, gap="small")
                                            sidebar_col1.button(label="📁 Export", on_click=self.export_workspace, use_container_width=True)
                                            sidebar_col2.button(label="🗑️ Delete", on_click=self.delete_workspace, key="btn_delete_workspace", use_container_width=True)
                                            st.divider()
        else:
            st.markdown("""<p style='font-size: 20px;'>You don't have any workspaces yet 🧪</p>""", unsafe_allow_html=True)


        with st.sidebar:
            
            self.new_workspace()

            self.import_workspace()

            st.markdown(f"""<div style='position: fixed; margin-left: 5px; bottom: 5px;'>
                        <div style='font-size: 16px; color: rgba(255, 255, 255, 0.4)'>Nuwa v{os.getenv('NUWA_VERSION')}</div>
                        </div>""", unsafe_allow_html=True)

dashboard = Dashboard()



