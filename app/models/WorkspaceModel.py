from typing import Optional
from pydantic import BaseModel, ConfigDict, ValidationError
from datetime import date

class WorkspaceModel(BaseModel):
    id: int
    workspace_name: str
    data_dir: str
    created: Optional[date]
    description: Optional[str]