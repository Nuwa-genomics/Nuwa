from typing import Optional
from pydantic import BaseModel, ConfigDict, ValidationError
from datetime import datetime

class AdataModel(BaseModel):
    id: int
    name: str
    adata: bytes
    created: Optional[datetime]