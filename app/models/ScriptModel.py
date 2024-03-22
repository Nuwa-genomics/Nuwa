from typing import Optional
from pydantic import BaseModel, ConfigDict, ValidationError
from datetime import date

class ScriptModel(BaseModel):
    adata_id: int
    id: int
    script: str
    language: str
    created: Optional[date]

