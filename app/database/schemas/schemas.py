from ..database import Base
from sqlalchemy import Column, Integer, String, Boolean, ForeignKey
from sqlalchemy.orm import relationship
from sqlalchemy.sql.expression import text
from sqlalchemy.sql.sqltypes import TIMESTAMP

class Workspaces(Base):
    __tablename__ = "workspaces"

    id = Column(Integer, primary_key=True, nullable=False, autoincrement=True)
    workspace_name = Column(String, nullable=False, unique=True)
    data_dir = Column(String, nullable=False, unique=True)
    created = Column(TIMESTAMP(timezone=True), nullable=False, server_default=text('now()'))
    description = Column(String, nullable=True)

class Adata(Base):
    __tablename__ = "adata"

    work_id = Column(Integer, ForeignKey("workspaces.id", ondelete="CASCADE"), nullable=False)
    id = Column(Integer, primary_key=True, nullable=False, autoincrement=True)
    adata_name = Column(String, nullable=False)
    filename = Column(String, nullable=False, unique=True) 
    notes = Column(String, nullable=True)
    created = Column(TIMESTAMP(timezone=True), nullable=False, server_default=text('now()'))
    
class Scripts(Base):
    __tablename__ = "scripts"
    
    adata_id = Column(Integer, ForeignKey("adata.id", ondelete="CASCADE"), nullable=False)
    id = Column(Integer, primary_key=True, nullable=False, autoincrement=True)
    script = Column(String, nullable=False)
    language = Column(String, nullable=False)
    created = Column(TIMESTAMP(timezone=True), nullable=False, server_default=text('now()'))