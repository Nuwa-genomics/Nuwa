import hashlib
import pickle
import os
import streamlit as st
from database.schemas import schemas
from database.database import SessionLocal
from models.SessionModel import SessionModel

def cache_data_to_session(description: str = None):
    try:
        
        state = {}
        for key in st.session_state:
            #streamlit form and button can't be set using session state, so remove them here
            if not (key.__contains__("btn") or key.__contains__("toggle") or key.__contains__("FormSubmitter")):
                state[key] = st.session_state[key]
    
        state['adata_state'].conn = None
        state['script_state'].conn = None

        # create hash to be used as filename
        encoded = pickle.dumps(state)
        hash = hashlib.md5()
        hash.update(encoded)
        state_hash = hash.hexdigest()
        
        filepath = os.path.join(os.getenv('WORKDIR'), 'tmp', state_hash)
   
        # Write to file
        dbfile = open(filepath, mode='wb')
        pickle.dump(state, dbfile)
        dbfile.close()

        # python doesn't copy the objects so db connection in state is destroyed. Add it back here
        st.session_state["adata_state"].conn = SessionLocal()
        st.session_state["script_state"].conn = SessionLocal()

        # write cache file to db
        conn = SessionLocal()
        new_session = schemas.Session(session_id=state_hash, adata_id=st.session_state.adata_state.current.id, filename=filepath, description=description)
        conn.add(new_session)

        # commit cache file to db
        conn.commit()
        
    except Exception as e:
        st.toast(e, icon="❌")

def load_data_from_cache(state_file):
    try:
        dbfile = open(os.path.join(os.getenv('WORKDIR'), 'tmp', state_file), 'rb')    
        session = pickle.load(dbfile)
        for key in session:
            if not (key.__contains__("FormSubmitter") or key.__contains__("file_uploader") or 
                key.__contains__("btn") or key.__contains__("toggle")):
                st.session_state[key] = session[key] # load in keys to session state
        dbfile.close()

        adata_state = st.session_state.adata_state
        adata_state.conn = SessionLocal()
        script_state = st.session_state.script_state
        script_state.conn = SessionLocal()

        st.session_state["adata_state"] = adata_state
        st.session_state["script_state"] = script_state

    except Exception as e:
        st.toast(e, icon="❌")
