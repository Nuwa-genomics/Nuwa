import streamlit as st
from components.sidebar import *
from models.AdataModel import AdataModel
import os
import pandas as pd

st.set_page_config(layout="wide", page_title='Nuwa', page_icon='🧬')

os.chdir('/app')

with open('css/common.css') as f:
    common_style = f"""
                <style>
                {f.read()}
                </style>
                """
    st.markdown(common_style, unsafe_allow_html=True)

class Integrate:
    def __init__(self):
        st.title("Integrate datasets")
        col1, col2, col3 = st.columns(3)
        with col1:
            self.ingest()
            self.bbknn()
        with col2:
            self.quick_map()
        with col3:
            self.umap()
            
        
    def ingest(self):
        with st.form(key="ingest_form"):
            st.subheader("Integrate with Ingest")
            st.markdown(f"""<div style='color: rgb(50, 168, 82); font-weight: bold''>{st.session_state.adata_ref.adata_name} → {st.session_state.adata_target.adata_name}</div>""", unsafe_allow_html=True)
            obs = st.multiselect(label="Obs", options=st.session_state.adata_ref.adata.obs_keys())
            subcol1, _, _, _ = st.columns(4)
            submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True, disabled=(not st.session_state.sync_genes))
            if submit_btn:
                with st.spinner(text="Integrating datasets"):
                    sc.tl.ingest(adata=st.session_state.adata_target.adata, adata_ref=st.session_state.adata_ref.adata, obs=obs)
            
    def bbknn(self):
        with st.form(key="bbknn_form"):
            st.subheader("BBKNN")
            subcol1, _, _, _ = st.columns(4)
            submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True, disabled=(not st.session_state.sync_genes))
            
    def quick_map(self):
        with st.form(key="quick_map_form"):
            st.subheader("Quick map")
            st.write("Map an attribute from one dataset to another.")
            row1_subcol1, row1_subcol2 = st.columns(2)
            row1_subcol1.selectbox(label="Reference adata", options=[adata.adata_name for adata in st.session_state.adata_state.load_adata(st.session_state.current_workspace.id)])
            row1_subcol2.multiselect(label="Obs", options=['obs1', 'obs2', 'obs3'], key="ms_adata_ref_obs")
            st.markdown("""<div style='text-align: center; color: rgb(235, 143, 52);'>Map to ⬇</div>""", unsafe_allow_html=True)
            row2_subcol1, row2_subcol2 = st.columns(2)
            row2_subcol1.selectbox(label="Target adata", options=[adata.adata_name for adata in st.session_state.adata_state.load_adata(st.session_state.current_workspace.id)])
            row2_subcol2.multiselect(label="Obs", options=['obs1', 'obs2', 'obs3'], key="ms_adata_obs")
            subcol1_btn, _, _, _ = st.columns(4)
            subcol1_btn.form_submit_button(label="Run", use_container_width=True, disabled=(not st.session_state.sync_genes))
    
            
    def umap(self):
        with st.form(key="integrate_umap_form"):
            st.subheader("UMAP")
            colors = st.multiselect(label="Color (obs)", options=st.session_state.adata_state.current.adata.obs_keys())
            container = st.container()
            subcol1, _, _, _ = st.columns(4)
            submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
            if submit_btn:
                with st.spinner(text="Computing umap"):
                    adata = st.session_state.adata_state.current.adata
                    sc.pp.neighbors(adata)
                    sc.tl.umap(adata)
                    #ax = sc.pl.umap(adata, color=colors)
                    #empty.pyplot(ax)
                    for color in colors:
                        df_umap = pd.DataFrame({'umap1': adata.obsm['X_umap'][:,0], 'umap2': adata.obsm['X_umap'][:,1], 'color': adata.obs[color]})
                        container.scatter_chart(data=df_umap, x='umap1', y='umap2', color='color', size=18)
            
    


sidebar = Sidebar()

sidebar.show(integrate=True)

sidebar.show_preview(integrate=True)
        
integrate = Integrate()

sidebar.export_script()

sidebar.delete_experiment_btn()