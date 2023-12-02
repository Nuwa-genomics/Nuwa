from scanpy import AnnData
from sqlalchemy import update
from database.schemas import schemas
from models.AdataModel import AdataModel
import scanpy as sc
import streamlit as st
from database.database import SessionLocal
from sqlalchemy.orm import Session
from utils.ScriptState import ScriptState
import os

class AdataState:
    def __init__(self, active: AdataModel):
        self.conn: Session = SessionLocal()
        self.conn
        self.workspace_id = active.work_id
        #set initial value for current adata. This won't include timestamp or id so needs to be reinitialised.
        self.current = active
        #insert active into db
        insert_warning = self.insert_record(active)
        if insert_warning != 0:
            st.warning(insert_warning)
        #check record is inserted
        db_adatas = self.conn.query(schemas.Adata).filter(schemas.Adata.work_id == active.work_id).filter(schemas.Adata.adata_name == active.adata_name)
        assert db_adatas.count() != 0
        #reinitialise current to include additional fields
        current: schemas.Adata = db_adatas.first()
        self.current = AdataModel(work_id=current.work_id, adata_name=current.adata_name, created=current.created, notes=current.notes, id=current.id, filename=current.filename)
        #add original adata to object
        self.current.adata = active.adata
        #set current to newly created adata
        self.current_index = self.get_index_of_current()
        #set script state
        st.session_state["script_state"] = ScriptState(adata_id=db_adatas.first().id)

    def switch_adata(self, adata_name):
        new_current = self.load_adata(workspace_id=self.workspace_id, adata_name=adata_name)
        if new_current != -1:
            self.current = new_current
            self.current_index = self.get_index_of_current()
            st.session_state["script_state"].switch_adata(new_current.id) #swap adata in script state
        else:
            return -1
        
    def get_adata_options(self):
        return [item.adata_name for item in self.load_adata(self.workspace_id)]
        
    def get_index_of_current(self):
        adata_list = self.load_adata(workspace_id=self.workspace_id)
        for i, adata in enumerate(adata_list):
            if adata.adata_name == self.current.adata_name:
                return i

    def load_adata(self, workspace_id, adata_name = None):
        #if no adata is provided fetch all
        if not adata_name:
            adatas = self.conn.query(schemas.Adata).filter(schemas.Adata.work_id == workspace_id).order_by(schemas.Adata.created).all()
            adata_list = [AdataModel(work_id=adata.work_id, id=adata.id, filename=adata.filename, notes=adata.notes, created=adata.created, adata_name=adata.adata_name) for adata in adatas]
            return adata_list
        
        adata = self.conn.query(schemas.Adata).filter(schemas.Adata.work_id == workspace_id).filter(schemas.Adata.adata_name == adata_name).first()
       
        if adata:
            return AdataModel(work_id=adata.work_id, id=adata.id, filename=adata.filename, notes=adata.notes, created=adata.created, adata_name=adata.adata_name, adata=sc.read_h5ad(adata.filename))
        else:
            st.toast("Couldn't load adata", icon="❌")
            return -1
        

    def insert_record(self, adata: AdataModel):
        try:
            if not adata.adata:
                if self.current == None:
                    raise Exception
                adata.adata = self.current.adata.copy() #if adata doesn't come from sender, add it here

            new_adata = schemas.Adata(
                work_id=adata.work_id,
                adata_name=adata.adata_name,
                filename=adata.filename,
                notes=adata.notes
            )

            sc.write(filename=adata.filename, adata=adata.adata)
            
            if os.path.exists(adata.filename):
                if self.conn.query(schemas.Adata).filter(schemas.Adata.adata_name == adata.adata_name).filter(schemas.Adata.work_id == adata.work_id).count() == 0:
                    self.conn.add(new_adata)
                    self.conn.commit()
                    st.toast("Created new adata", icon="✅")
                    return 0
                else:
                    return "Dataset already exists in workspace, using original."
            else:
                print("Error: file doesn't exist")
                raise Exception
        except Exception as e:
            print("Error: ", e)
            st.toast("Couldn't add experiment", icon="❌")

    def update_record(self, adata: AdataModel = None):
        try:
            if adata == None:
                adata = self.current

            update_query = self.conn.query(schemas.Adata).filter(schemas.Adata.adata_name == adata.adata_name).filter(schemas.Adata.work_id == adata.work_id)
            update_query.update({'notes': adata.notes, 'adata_name': adata.adata_name}) #some fields should not need to be updated so not included
            self.conn.commit()
        except Exception as e:
            st.error(e)

    def delete_record(self, adata_name: str = None):
        try:
            if adata_name == None:
                adata = self.current.adata_name

            num_of_records = self.conn.query(schemas.Adata).count()
            if num_of_records <= 1:
                raise Exception

            self.conn.query(schemas.Adata).filter(schemas.Adata.adata_name == adata_name).filter(schemas.Adata.work_id == self.workspace_id).delete()
            self.conn.commit()

            st.toast("Deleted Experiment", icon="✅")
        except Exception as e:
            print("Error: ", e)
            st.toast("Couldn't delete Experiment", icon="❌")

    

    