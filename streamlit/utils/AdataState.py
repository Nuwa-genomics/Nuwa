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
        self.adata_list: list(AdataModel) = []
        self.load_adata(workspace_id)
        self.current: AnnData = self.get_adata(active)

    def __iter__(self):
        for adata in self.adata_list:
            yield adata

    def switch_adata(self, adata_name):
        self.current = self.get_adata(adata_name)

    def load_adata(self, workspace_id):
        adatas = self.conn.query(schemas.Adata).filter(schemas.Adata.work_id == workspace_id).all()
        for adata in adatas:
            new_adata = AdataModel(
                work_id=adata.work_id,
                id=adata.id,
                adata_name=adata.adata_name,
                filename=adata.filename,
                created=adata.created,
                notes=adata.notes,
                adata=sc.read_h5ad(adata.filename)
            )
            self.adata_list.append(new_adata)
        return self.adata_list
        

    def add_adata(self, adata: AdataModel):
        adata.adata = self.current.copy() #adata doesn't come from sender, so add it here
        self.adata_list.append(adata)
        #write adata to db
        self.write_to_db(adata)
        #write adata to file
        sc.write(filename=adata.filename, adata=adata.adata)

    def write_to_db(self, adata: AdataModel = None):
        if adata == None:
            adata = self.current

        new_adata = schemas.Adata(
            work_id=adata.work_id,
            adata_name=adata.adata_name,
            filename=adata.filename,
            notes=adata.notes
        )
        self.conn.add(new_adata)
        self.conn.commit()
        self.conn.refresh(new_adata)

    def update_record(self, adata: AdataModel = None):
        if adata == None:
            adata = self.current

        update_query = self.conn.query(schemas.Adata).filter(schemas.Adata.adata_name == adata.adata_name)
        update_query.update({'notes': adata.notes, 'filename': adata.filename, 'work_id': adata.work_id, 'id': adata.id, 'created': adata.created, 'adata_name': adata.adata_name})
        self.conn.commit()


    def get_adata(self, adata_name=None) -> AdataModel:
        if not adata_name:
            return self.adata_list
        
        for adata in self.adata_list:
            if adata.adata_name == adata_name:
                return adata

    

    