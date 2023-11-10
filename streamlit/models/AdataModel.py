from typing import Optional
from anndata import AnnData
from datetime import date
import pydantic
from pydantic import BaseModel as PydanticBaseModel

class BaseModel(PydanticBaseModel):
    class Config:
        arbitrary_types_allowed = True

class AdataModel(BaseModel):
    work_id: int
    id: Optional[int] #Uses autoincrement when not included
    adata_name: str
    filename: str
    adata: Optional[AnnData]
    notes: Optional[str]
    created: Optional[date]
