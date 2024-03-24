from typing import Optional
from pydantic import BaseModel, ConfigDict, ValidationError
from datetime import date

class SessionModel(BaseModel):
    id: Optional[int] #Uses autoincrement when not included
    adata_id: int
    session_id: str
    filename: str
    description: Optional[str]
    created: Optional[date]