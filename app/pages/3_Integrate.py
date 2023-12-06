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
        col1, col2, col3 = st.columns(3, gap="medium")
        with col1:
            self.auto_integrate_recipies()
            self.ingest()
            self.bbknn()
        with col2:
            self.quick_map()
            self.concat()
        with col3:
            self.umap()
            
            
    def auto_integrate_recipies(self):
        with st.form(key="auto_integrate_recipies_form"):
            st.subheader("Autointegrate recipies")
            recipie = st.selectbox(label="Recipie", options=['rec'])
            subcol1, _, _, _ = st.columns(4)
            submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True, disabled=(not st.session_state.sync_genes))
            if submit_btn:
                with st.spinner(text="Applying integration recipie"):
                    print("hi")
            
        
    def ingest(self):
        with st.form(key="ingest_form"):
            st.subheader("Integrate with Ingest", help="Integrates embeddings and annotations of an adata with a reference dataset adata_ref through projecting on a PCA (or alternate model) that has been fitted on the reference data. The function uses a knn classifier for mapping labels and the UMAP package [McInnes18] for mapping the embeddings.")
            st.markdown(f"""<div style='color: rgb(50, 168, 82); font-weight: bold''>{st.session_state.adata_ref.adata_name} → {st.session_state.adata_target.adata_name}</div>""", unsafe_allow_html=True)
            obs = st.multiselect(label="Obs", options=st.session_state.adata_ref.adata.obs_keys())
            subcol1, _, _, _ = st.columns(4)
            submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True, disabled=(not st.session_state.sync_genes))
            if submit_btn:
                try:
                    with st.spinner(text="Integrating datasets"):
                        sc.pp.pca(st.session_state.adata_ref.adata)
                        sc.pp.neighbors(st.session_state.adata_ref.adata)
                        sc.tl.umap(st.session_state.adata_ref.adata)
                        sc.tl.ingest(adata=st.session_state.adata_target.adata, adata_ref=st.session_state.adata_ref.adata, obs=obs)
                        st.toast("Integration complete", icon="✅")
                except Exception as e:
                    st.toast("Failed to integrate datasets", icon="❌")
                    st.error(e)
                    print("Error: ", e)
            
    def bbknn(self):
        with st.form(key="bbknn_form"):
            st.subheader("BBKNN", help="Batch balanced kNN alters the kNN procedure to identify each cell’s top neighbours in each batch separately instead of the entire cell pool with no accounting for batch. The nearest neighbours for each batch are then merged to create a final list of neighbours for the cell. Aligns batches in a quick and lightweight manner.")
            st.write(f"Apply to {st.session_state.adata_state.current.adata_name}")
            batch_key = st.selectbox(label="Batch key", options=st.session_state.adata_state.current.adata.obs_keys())
            subcol1, _, _, _ = st.columns(4)
            submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
            if submit_btn:
                with st.spinner(text="Computing PCA"):
                    sc.tl.pca(st.session_state.adata_state.current.adata)
                with st.spinner(text="Applying BBKNN"): 
                    sc.external.pp.bbknn(st.session_state.adata_state.current.adata, batch_key=batch_key)
                    
            
    def quick_map(self):
        with st.form(key="quick_map_form"):
            st.subheader("Quick map")
            st.write("Map an attribute from one dataset to another.")
            
            row1_subcol1, row1_subcol2 = st.columns(2)
            source_dataset_name = row1_subcol1.selectbox(label="Reference adata", options=[adata.adata_name for adata in st.session_state.adata_state.load_adata(st.session_state.current_workspace.id)])
            row1_subcol2.text_input(label="Source attributes", placeholder="e.g. uns.louvain_colors", key="ti_adata_src_atts")
            
            st.markdown("""<p style='text-align: center; color: rgb(235, 143, 52);'>Map to ⬇</p>""", unsafe_allow_html=True)
            
            dest_dataset_name = st.selectbox(label="Target adata", options=[adata.adata_name for adata in st.session_state.adata_state.load_adata(st.session_state.current_workspace.id)])
            #row2_subcol2.text_input(label="Destination attributes", placeholder="e.g. uns.louvain_colors", key="ti_adata_dest_atts")

            
            subcol1_btn, _, _, _ = st.columns(4)
            submit_btn = subcol1_btn.form_submit_button(label="Run", use_container_width=True, disabled=(not st.session_state.sync_genes))
            
            if submit_btn:
                try:
                    source_adata = st.session_state.adata_state.load_adata(st.session_state.current_workspace.id, source_dataset_name).adata
                    dest_adata = st.session_state.adata_state.load_adata(st.session_state.current_workspace.id, dest_dataset_name).adata
                    
                    st.session_state.adata_state.insert_record(AdataModel(
                        work_id=st.session_state.current_workspace.id, adata=dest_adata, 
                        filename=os.path.join(os.getenv('WORKDIR'), "adata", f"{dest_dataset_name}.h5ad"), adata_name=dest_dataset_name)
                    )
                    
                    st.toast("Successfully mapped requested fields into dataset.", icon="✅")
                except Exception as e:
                    st.toast("Failed to map attributes into dataset.", icon="❌")
                
                
            
            
    def concat(self):
        with st.form(key="concat_form"):
            st.subheader("Concatenate datasets")
            subcol1, subcol2 = st.columns(2)
            adata1_name = subcol1.selectbox(label="Dataset 1", options=[adata.adata_name for adata in st.session_state.adata_state.load_adata(st.session_state.current_workspace.id)])
            adata2_name = subcol2.selectbox(label="Dataset 2", options=[adata.adata_name for adata in st.session_state.adata_state.load_adata(st.session_state.current_workspace.id)])
            batch_cat1 = subcol1.text_input(label="Batch category 1")
            batch_cat2 = subcol2.text_input(label="Batch category 2")
            batch_key = st.text_input(label="Batch key", value="batch")
            empty = st.empty()
            subcol1_btn, _, _, _ = st.columns(4)
            submit_btn = subcol1_btn.form_submit_button(label="Run", use_container_width=True)
            if submit_btn:
                with st.spinner(text="Concatenating datasets"):
                    try:
                        if batch_cat1 == "" or batch_cat2 == "":
                            st.toast("Batch categories cannot be null", icon="❌")
                        else:
                            adata1: AnnData = st.session_state.adata_state.load_adata(st.session_state.current_workspace.id, adata1_name).adata
                            adata2: AnnData = st.session_state.adata_state.load_adata(st.session_state.current_workspace.id, adata2_name).adata
                            adata_concat = adata1.concatenate(adata2, batch_key=batch_key, batch_categories=[batch_cat1, batch_cat2])

                            #insert as new adata
                            adata_name = f"concat_{adata1_name}_{adata2_name}"
                            st.session_state.adata_state.insert_record(AdataModel(
                                work_id=st.session_state.current_workspace.id, adata=adata_concat, 
                                filename=os.path.join(os.getenv('WORKDIR'), "adata", f"{adata_name}.h5ad"), adata_name=adata_name)
                            )
                            st.toast("Successfully concatenated dataframes", icon="✅")
                    except Exception as e:
                        st.toast("Couldn't concatenate dataframes", icon="❌")
                        empty.error(e)
                        
                    
                
            
    
            
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