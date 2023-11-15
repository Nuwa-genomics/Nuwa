from tests.test_preprocess import Test_Preprocess
from tests.test_dashboard import Test_Dashboard
from tests.test_upload import Test_Upload
from tests.test_create_model import Test_Create_Model
from tests.test_train import Test_Train
from tests.test_cluster_analysis import Test_Cluster_Analysis
from tests.test_tranjectory_inference import Test_Trajectory_Inference
from tests.test_spatial_transcriptomics import Test_Spatial_Transcriptomics

from database.database import SessionLocal
from sqlalchemy.orm import Session
from utils.AdataState import AdataState
from database.schemas import schemas
from models.WorkspaceModel import WorkspaceModel
import os
from random import randrange
import shutil
from datetime import datetime, date

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

try:

    conn: Session = SessionLocal()

    print()
    print(f"{bcolors.BOLD}===============Testing Dashboard===============")
    workspace_name = f"test_workspace_{randrange(1, 1000000)}"
    dashboard_test = Test_Dashboard(workspace_name=workspace_name)
    dashboard_state = dashboard_test.get_final_session_state()
    print()
    print(f"{bcolors.OKGREEN}TEST PASSED{bcolors.ENDC}")

    print()
    print(f"{bcolors.BOLD}===============Testing Upload===============")
    upload_test = Test_Upload(session_state=dashboard_state)
    upload_state = upload_test.get_final_session_state()
    print()
    print(f"{bcolors.OKGREEN}TEST PASSED{bcolors.ENDC}")

    print()
    print(f"{bcolors.BOLD}===============Testing Preprocess===============")
    pp_test = Test_Preprocess(session_state=upload_state)
    pp_state = pp_test.get_final_session_state()
    print()
    print(f"{bcolors.OKGREEN}TEST PASSED{bcolors.ENDC}")

    print()
    print(f"{bcolors.BOLD}===============Testing Create Model===============")
    create_model_test = Test_Create_Model(session_state=pp_state)
    create_model_state = create_model_test.get_final_session_state()
    print()
    print(f"{bcolors.OKGREEN}TEST PASSED{bcolors.ENDC}")

    print()
    print(f"{bcolors.BOLD}===============Testing Train===============")
    train_test = Test_Train(session_state=create_model_state)
    train_state = train_test.get_final_session_state()
    print()
    print(f"{bcolors.OKGREEN}TEST PASSED{bcolors.ENDC}")

    print()
    print(f"{bcolors.BOLD}===============Testing Cluster Analysis===============")
    cluster_analysis_test = Test_Cluster_Analysis(session_state=train_state)
    cluster_analysis_state = cluster_analysis_test.get_final_session_state()
    print()
    print(f"{bcolors.OKGREEN}TEST PASSED{bcolors.ENDC}")

    print()
    print(f"{bcolors.BOLD}===============Testing Trajectory Inference===============")
    trajectory_inference_test = Test_Trajectory_Inference(session_state=pp_state)
    trajectory_inference_state = trajectory_inference_test.get_final_session_state()
    print()
    print(f"{bcolors.OKGREEN}TEST PASSED{bcolors.ENDC}")

    print()
    print(f"{bcolors.BOLD}===============Testing Spatial Transcriptomics===============")
    spatial_transcriptomics_test = Test_Spatial_Transcriptomics(session_state=pp_state)
    spatial_transcriptomics_state = spatial_transcriptomics_test.get_final_session_state()
    print()
    print(f"{bcolors.OKGREEN}TEST PASSED{bcolors.ENDC}")
    
    print()
    print("----------------------------------------------")
    print(f"Tests completed: {date.today().strftime('%B %d, %Y')} at {datetime.now().strftime('%H:%M:%S')}\n")
    print(f"{bcolors.BOLD}ALL TESTS SUCCESSFUL.{bcolors.ENDC}")
    print()
    print(f"Removing test data...{bcolors.ENDC}")

except Exception as e:
    print()
    print(f"{bcolors.FAIL}TEST FAILED{bcolors.ENDC}")
    print(e)

finally:

    #tear down

    #remove test records from db
    conn.query(schemas.Workspaces).filter(schemas.Workspaces.workspace_name == workspace_name).delete()
    conn.commit()

    #remove files
    workspace_dir = os.getenv('WORKDIR')
    shutil.rmtree(workspace_dir, ignore_errors=True)
    
    print()
    print(f"{bcolors.BOLD}Tests Complete!{bcolors.ENDC}")
    print("----------------------------------------------")



