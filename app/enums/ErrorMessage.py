from enum import Enum
import streamlit as st

class ErrorMessage(Enum):
    ADATA_NOT_FOUND = "Couldn't find adata"
    AT_LEAST_ONE_EXPERIMENT_REQUIRED = "Workspace must contain at least one experiment"
    CANNOT_DELETE_EXPERIMENT = "Couldn't delete experiment"
    CANNOT_UPDATE_EXPERIMENT = "Couldn't update experiment"
    CANNOT_ADD_EXPERIMENT = "Couldn't add experiment"
    FILE_DOES_NOT_EXIST = "File does not exist"
    NOTES_FAILED_TO_SAVE = "Notes failed to save"

class WarningMessage(Enum):
    DATASET_ALREADY_EXISTS = "Dataset already exists in workspace, using original."
