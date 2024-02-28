import hashlib
import pickle
import os
import streamlit as st
from database.database import SessionLocal

def cache_data_to_session():
    try:
        
        state = {}
        for key in st.session_state:
            #streamlit form and button can't be set using session state, so remove them here
            if not (key.__contains__("btn") or key.__contains__("toggle") or key.__contains__("form")):
                state[key] = st.session_state[key]
    
        state['adata_state'].conn = None
        state['script_state'].conn = None

        # create hash to be used as filename
        encoded = pickle.dumps(state)
        hash = hashlib.md5()
        hash.update(encoded)
        state_hash = hash.hexdigest()

        # Write to file
        dbfile = open(os.path.join(os.getenv('WORKDIR'), 'tmp', state_hash), 'wb')
        pickle.dump(state, dbfile)
        #python doesn't copy the objects so db connection in state is destroyed. Add it back here
        st.session_state["adata_state"].conn = SessionLocal()
        st.session_state["script_state"].conn = SessionLocal()
        dbfile.close()
    except Exception as e:
        st.toast(e, icon="❌")

def load_data_from_cache(state_file):
    try:
        dbfile = open(os.path.join(os.getenv('TMP_DIR'), state_file), 'rb')    
        session = pickle.load(dbfile)
        for key in session:
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
