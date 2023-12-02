from typing import Optional
from pydantic import BaseModel, ConfigDict, ValidationError
from datetime import date

class ScriptModel(BaseModel):
    adata_id: int
    id: int
    script: str
    created: Optional[date]