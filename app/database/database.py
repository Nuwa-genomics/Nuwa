from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
import os
from dotenv import load_dotenv

load_dotenv()

SQLALCHEMY_DATABASE_URL = f"postgresql://\
{os.getenv('POSTGRES_USER')}:\
{os.getenv('POSTGRES_PASSWORD')}@\
{os.getenv('POSTGRES_HOST')}:{os.getenv('POSTGRES_PORT')}/\
{os.getenv('POSTGRES_DB')}"

# Dependency
def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()

engine = create_engine(SQLALCHEMY_DATABASE_URL)

SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

Base = declarative_base()