from tests.test_preprocess import Test_Preprocess
from tests.test_dashboard import Test_Dashboard
from tests.test_upload import Test_Upload

from database.database import SessionLocal
from sqlalchemy.orm import Session
from utils.AdataState import AdataState
from database.schemas import schemas
from models.WorkspaceModel import WorkspaceModel
import os
from random import randrange
import shutil

conn: Session = SessionLocal()

workspace_name = f"test_workspace_{randrange(1, 1000000)}"
dashboard_test = Test_Dashboard(workspace_name=workspace_name)
dashboard_state = dashboard_test.get_final_session_state()

upload_test = Test_Upload(session_state=dashboard_state)
upload_state = upload_test.get_final_session_state()

#pp_test = Test_Preprocess(session_state=upload_state)
#pp_state = pp_test.get_final_session_state()

#tear down

#remove test records from db
conn.query(schemas.Workspaces).filter(schemas.Workspaces.workspace_name == workspace_name).delete()
conn.commit()

#remove files
workspace_dir = os.getenv('WORKDIR')
shutil.rmtree(workspace_dir, ignore_errors=True)



