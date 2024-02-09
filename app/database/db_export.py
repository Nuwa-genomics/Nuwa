from database.schemas import schemas
from database.database import *
from database.schemas.schemas import *
from sqlalchemy.orm import Session
from sqlalchemy.exc import IntegrityError
from os import PathLike
import streamlit as st
import json


def export_db_as_json(workspace_id: int, output_file: str | PathLike):

    try:
        conn: Session = SessionLocal()
        
        # create workspace dictionary
        workspace_dict = conn.query(schemas.Workspaces).filter(schemas.Workspaces.id == workspace_id).first().__dict__

        # create adatas dict
        adatas_dict = conn.query(schemas.Adata).filter(schemas.Adata.work_id == workspace_id).all()
        if workspace_dict.keys().__contains__("_sa_instance_state"):
            workspace_dict.pop("_sa_instance_state")

        # iterate adatas
        adatas = []
        for adata in adatas_dict:
            adata_dict = adata.__dict__
            if adata_dict.keys().__contains__("_sa_instance_state"):
                adata_dict.pop("_sa_instance_state")

            # iterate scripts
            scripts = []
            scripts_dict = conn.query(schemas.Scripts).filter(schemas.Scripts.adata_id == adata.id).all()
            for script in scripts_dict:
                script_dict = script.__dict__
                if script_dict.keys().__contains__("_sa_instance_state"):
                    script_dict.pop("_sa_instance_state")

                scripts.append(script_dict)
            
            # add scripts to adata
            adata_dict["scripts"] = scripts

            adatas.append(adata_dict)

        # add adatas to workspace
        workspace_dict["adatas"] = adatas

        with open(output_file, 'w') as f:
            json.dump(workspace_dict, f, default=str, indent=4)

    except Exception as e:
        st.toast(f"Error: {e}", icon="❌")
        



def import_db_from_json(json_filepath: str | PathLike):

    try:
        conn: Session = SessionLocal()

        if not os.path.exists(json_filepath):
            st.toast("Cannot import: no database file found", icon="❌")
            return -1
        
        with open(json_filepath, mode="r") as f:
            json_file = json.loads(f.read())

        # add new workspace
        new_workspace = Workspaces(id=json_file['id'], 
            workspace_name=json_file['workspace_name'], 
            description=json_file['description'], 
            data_dir=json_file['data_dir'],
            created=json_file['created']
        )

        conn.add(new_workspace)
        conn.commit()
        conn.refresh(new_workspace)

        # add adatas
        for adata in json_file['adatas']:
            new_adata = Adata(
                work_id=adata["work_id"],
                id=adata["id"],
                adata_name=adata["adata_name"],
                filename=adata["filename"],
                notes=adata["notes"],
                created=adata["created"]
            )

            conn.add(new_adata)
            conn.commit()
            conn.refresh(new_adata)

            for script in adata["scripts"]:
                new_script = Scripts(
                    adata_id=script["adata_id"],
                    id=script["id"],
                    script=script["script"],
                    created=script["created"]
                )

                conn.add(new_script)
                conn.commit()
                conn.refresh(new_script)

    except IntegrityError as unique_ke:
        st.toast("Cannot import: record already exists", icon="❌")
        return -1


    return 0
