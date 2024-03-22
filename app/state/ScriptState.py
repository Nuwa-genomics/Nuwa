from typing import List
from scanpy import AnnData
from sqlalchemy import update
from database.schemas import schemas
from models.ScriptModel import ScriptModel
from enums.Language import Language
import scanpy as sc
import streamlit as st
from database.database import SessionLocal
from sqlalchemy.orm import Session


class ScriptState:
    def __init__(self, adata_id):
        self.conn: Session = SessionLocal()
        self.adata_id = adata_id
        self.add_imports()
        
        
    def load_script(self, language = "python") -> List[str]:
        scripts = self.conn.query(schemas.Scripts) \
        .filter(schemas.Scripts.adata_id == self.adata_id) \
        .filter(schemas.Scripts.language == language) \
        .order_by(schemas.Scripts.created).all()
        if scripts:
            return [script.script for script in scripts]
    
    def add_script(self, script, language: Language | str):
        
        if isinstance(language, Language):
            language = language.value # convert to string

        new_script = schemas.Scripts(
            adata_id=self.adata_id,
            script=script,
            language=language
        )
        self.conn.add(new_script)
        self.conn.commit()
        self.conn.refresh(new_script)
        
    def add_imports(self):
        # Add imports to python scripts
        # add imports here if scanpy is not present
        if self.conn.query(schemas.Scripts) \
        .filter(schemas.Scripts.script.contains("import scanpy as sc")) \
        .filter(schemas.Scripts.adata_id == self.adata_id) \
        .filter(schemas.Scripts.language == "python").count() == 0:
            self.add_script(script="import scanpy as sc\nimport numpy as np\nimport matplotlib.pyplot as plt\nimport squidpy as sq", language="python")

        # Add imports to R scripts
        # add R imports here if Seurat is not present
        if self.conn.query(schemas.Scripts) \
        .filter(schemas.Scripts.script.contains("library(Seurat)")) \
        .filter(schemas.Scripts.adata_id == self.adata_id) \
        .filter(schemas.Scripts.language == "R").count() == 0:
            self.add_script(script="library(Seurat)", language="R")
            
        
    def switch_adata(self, adata_id):
        self.adata_id = adata_id
        self.add_imports()
        