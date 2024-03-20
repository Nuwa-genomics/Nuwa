from typing import List
from scanpy import AnnData
from sqlalchemy import update
from database.schemas import schemas
from models.AdataModel import AdataModel
import scanpy as sc
import streamlit as st
from database.database import SessionLocal
from sqlalchemy.orm import Session
from state.ScriptState import ScriptState
import os
from models.ErrorMessage import ErrorMessage, WarningMessage

class AdataState:
    def __init__(self, active: AdataModel, insert_into_db=True):
        self.conn: Session = SessionLocal()
        self.workspace_id = active.work_id
        #set initial value for current adata. This won't include timestamp or id so needs to be reinitialised.
        self.current = active
        #insert active into db
        if insert_into_db:
            self.insert_record(active)
            #check record is inserted
            db_adatas = self.conn.query(schemas.Adata) \
            .filter(schemas.Adata.work_id == active.work_id) \
            .filter(schemas.Adata.adata_name == active.adata_name)
            if db_adatas.count() == 0:
                st.toast(ErrorMessage.ADATA_NOT_FOUND.value, icon="❌")
            #reinitialise current to include additional fields
            current: schemas.Adata = db_adatas.first()
            self.current = AdataModel(work_id=current.work_id, adata_name=current.adata_name, created=current.created, notes=current.notes, id=current.id, filename=current.filename)

        # add original adata to object
        self.current.adata = active.adata
        #set current to newly created adata
        self.current_index = self.get_index_of_current()
        #set script state
        st.session_state["script_state"] = ScriptState(adata_id=self.current.id)


    def switch_adata(self, adata_name):
        try:
            new_current = self.load_adata(workspace_id=self.workspace_id, adata_name=adata_name)

            if not new_current:
                return
    
            self.current = new_current
            self.current_index = self.get_index_of_current()
            st.session_state["script_state"].switch_adata(new_current.id) #swap adata in script state

        except Exception as e:
            st.toast(e, icon="❌")
        

    def get_adata_options(self) -> List[str]:
        try:
            adata_list = self.load_adata(self.workspace_id)

            if not adata_list:
                return
            
            return [item.adata_name for item in adata_list]
        
        except Exception as e:
            st.toast(e, icon="❌")
    
        
    def get_index_of_current(self) -> int:
        try:
            adata_list = self.load_adata(workspace_id=self.workspace_id)

            if not adata_list:
                return
            
            for i, adata in enumerate(adata_list):
                if adata.adata_name == self.current.adata_name:
                    return i
        
        except Exception as e:
            st.toast(e, icon="❌")

    def load_adata(self, workspace_id, adata_name = None) -> AdataModel | List[AdataModel]:
        """Returns an adata model of one is specified, else return all adata models."""
        try:
            #if no adata is provided fetch all
            if not adata_name:
                adatas = self.conn.query(schemas.Adata).filter(schemas.Adata.work_id == workspace_id).order_by(schemas.Adata.created).all()
                adata_list = [AdataModel(
                    work_id=adata.work_id, 
                    id=adata.id, 
                    filename=adata.filename, 
                    notes=adata.notes, 
                    created=adata.created, 
                    adata_name=adata.adata_name
                ) for adata in adatas]

                return adata_list
            
            adata = self.conn.query(schemas.Adata) \
            .filter(schemas.Adata.work_id == workspace_id) \
            .filter(schemas.Adata.adata_name == adata_name) \
            .first()
        
            if not adata:
                raise Exception(ErrorMessage.ADATA_NOT_FOUND.value)
            
            return AdataModel(
                work_id=adata.work_id, 
                id=adata.id, 
                filename=adata.filename, 
                notes=adata.notes, 
                created=adata.created, 
                adata_name=adata.adata_name, 
                adata=sc.read_h5ad(adata.filename)
            )
        
        except Exception as e:
            st.toast(e, icon="❌")
    

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
                if self.conn.query(schemas.Adata) \
                    .filter(schemas.Adata.adata_name == adata.adata_name) \
                    .filter(schemas.Adata.work_id == adata.work_id) \
                    .count() == 0:
                    self.conn.add(new_adata)
                    self.conn.commit()
                else:
                    st.toast(WarningMessage.DATASET_ALREADY_EXISTS.value, icon="⚠️")
            else:
                print(f"Error: {ErrorMessage.FILE_DOES_NOT_EXIST}")
                st.toast(ErrorMessage.FILE_DOES_NOT_EXIST, icon="❌")
        except Exception as e:
            print("Error: ", e)
            st.toast(ErrorMessage.CANNOT_ADD_EXPERIMENT, icon="❌")

    def update_record(self, adata: AdataModel = None):
        try:
            if adata == None:
                adata = self.current

            update_query = self.conn.query(schemas.Adata) \
            .filter(schemas.Adata.adata_name == adata.adata_name) \
            .filter(schemas.Adata.work_id == adata.work_id)
            update_query.update({'notes': adata.notes, 'adata_name': adata.adata_name}) #some fields should not need to be updated so not included
            self.conn.commit()
        except Exception as e:
            print("Error: ", e)
            st.toast(ErrorMessage.CANNOT_UPDATE_EXPERIMENT.value, icon="❌")

    def delete_record(self, adata_name: str = None):
        try:
            if adata_name == None:
                adata_name = self.current.adata_name

            # Must be at least 1 experiment in workspace
            num_of_records = self.conn.query(schemas.Adata) \
            .filter(schemas.Adata.work_id == self.workspace_id) \
            .count()

            if num_of_records <= 1:
                st.toast(ErrorMessage.AT_LEAST_ONE_EXPERIMENT_REQUIRED.value, icon="❌")

            self.conn.query(schemas.Adata) \
            .filter(schemas.Adata.adata_name == adata_name) \
            .filter(schemas.Adata.work_id == self.workspace_id) \
            .delete()
            self.conn.commit()

            st.toast("Deleted Experiment", icon="✅")
        except Exception as e:
            print("Error: ", e)
            st.toast(ErrorMessage.CANNOT_DELETE_EXPERIMENT.value, icon="❌")

    

    