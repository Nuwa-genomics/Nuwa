from scanpy import AnnData
from sqlalchemy import update
from database.schemas import schemas
from models.ScriptModel import ScriptModel
import scanpy as sc
import streamlit as st
from database.database import SessionLocal
from sqlalchemy.orm import Session
import os

class ScriptState:
    def __init__(self, adata_id):
        self.conn: Session = SessionLocal()
        self.adata_id = adata_id
        self.add_imports()
        
        
    def load_script(self):
        scripts = self.conn.query(schemas.Scripts).filter(schemas.Scripts.adata_id == self.adata_id).order_by(schemas.Scripts.created).all()
        if scripts:
            return [script.script for script in scripts]
    
    def add_script(self, script):
        new_script = schemas.Scripts(
            adata_id=self.adata_id,
            script=script
        )
        self.conn.add(new_script)
        self.conn.commit()
        self.conn.refresh(new_script)
        
    def add_imports(self):
        #add scanpy import here if not present
        if self.conn.query(schemas.Scripts).filter(schemas.Scripts.script.contains("import scanpy as sc")).filter(schemas.Scripts.adata_id == self.adata_id).count() == 0:
            self.add_script("import scanpy as sc\nimport numpy as np\nimport matplotlib.pyplot as plt\nimport squidpy as sq\n")
        
    def switch_adata(self, adata_id):
        self.adata_id = adata_id
        self.add_imports()
        