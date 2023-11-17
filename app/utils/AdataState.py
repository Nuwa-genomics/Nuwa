from scanpy import AnnData
from sqlalchemy import update
from database.schemas import schemas
from models.AdataModel import AdataModel
import scanpy as sc
import streamlit as st
from database.database import SessionLocal
from sqlalchemy.orm import Session
import os

class AdataState:
    def __init__(self, workspace_id, active="adata_raw"):
        self.conn: Session = SessionLocal()
        self.workspace_id = workspace_id
        raw = self.conn.query(schemas.Adata).filter(schemas.Adata.adata_name == "adata_raw").filter(schemas.Adata.work_id == workspace_id).first()
        raw_filename = os.path.join(os.getenv('WORKDIR'), "adata", "adata_raw.h5ad")
        if not raw:
            #create new record
            new_adata = schemas.Adata(
                work_id=workspace_id,
                adata_name="adata_raw",
                filename=raw_filename,
                notes=""
            )
            self.conn.add(new_adata)
            self.conn.commit()
            self.conn.refresh(new_adata)
        
        raw = self.conn.query(schemas.Adata).filter(schemas.Adata.adata_name == "adata_raw").filter(schemas.Adata.work_id == workspace_id).first()
        self.raw: AdataModel = AdataModel(work_id=raw.work_id, id=raw.id, filename=raw.filename, notes=raw.notes, 
            adata=sc.read_h5ad(raw_filename), created=raw.created, adata_name=raw.adata_name)
        self.current: AdataModel = self.raw

    def __iter__(self):
        for adata in self.adata_list:
            yield adata

    def switch_adata(self, adata_name):
        new_current = self.load_adata(workspace_id=self.workspace_id, adata_name=adata_name)
        if new_current != -1:
            self.current = new_current
        else:
            return -1

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

    def add_adata(self, adata: AdataModel):
        try:
            if not adata.adata:
                if self.current == None:
                    raise Exception
                adata.adata = self.current.adata.copy() #adata doesn't come from sender, so add it here
            
            self.insert_record(adata)
        except Exception as e:
            st.error(e)
            print("Error: ", e)
        

    def insert_record(self, adata: AdataModel = None):
        try:
            if adata == None:
                adata = self.current

            new_adata = schemas.Adata(
                work_id=adata.work_id,
                adata_name=adata.adata_name,
                filename=adata.filename,
                notes=adata.notes
            )

            sc.write(filename=adata.filename, adata=adata.adata)
            
            
            if os.path.exists(adata.filename):
                self.conn.add(new_adata)
                self.conn.commit()
                self.conn.refresh(new_adata)
                st.toast("Created new adata", icon="✅")
            else:
                print("Error: file doesn't exist")
                raise Exception
        except Exception as e:
            print("Error: ", e)
            st.toast("Couldn't add experiment", icon="❌")

    def update_record(self, adata: AdataModel = None):
        if adata == None:
            adata = self.current

        update_query = self.conn.query(schemas.Adata).filter(schemas.Adata.adata_name == adata.adata_name).filter(schemas.Adata.work_id == adata.work_id)
        update_query.update({'notes': adata.notes, 'created': adata.created, 'adata_name': adata.adata_name}) #some fields should not need to be updated so not included
        self.conn.commit()

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

    

    