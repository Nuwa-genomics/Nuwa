from ..database import Base
from sqlalchemy import Column, Integer, String, Boolean, ForeignKey
from sqlalchemy.orm import relationship
from sqlalchemy.sql.expression import text
from sqlalchemy.sql.sqltypes import TIMESTAMP

class Workspaces(Base):
    __tablename__ = "workspaces"

    id = Column(Integer, primary_key=True, nullable=False, autoincrement=True)
    workspace_name = Column(String, nullable=False)
    data_dir = Column(String, nullable=False)
    created = Column(TIMESTAMP(timezone=True), nullable=False, server_default=text('now()'))
    description = Column(String, nullable=True)

class Adata(Base):
    __tablename__ = "adata"

    work_id = Column(Integer, ForeignKey("workspaces.id", ondelete="CASCADE"), nullable=False, primary_key=True)
    id = Column(Integer, primary_key=True, nullable=False, autoincrement=True)
    adata_name = Column(String, nullable=False)
    filename = Column(String, nullable=False) 
    notes = Column(String, nullable=True)
    created = Column(TIMESTAMP(timezone=True), nullable=False, server_default=text('now()'))