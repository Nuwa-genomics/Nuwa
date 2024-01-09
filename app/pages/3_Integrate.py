from anndata import AnnData
import streamlit as st
from components.sidebar import *
from models.AdataModel import AdataModel
import os
import pandas as pd
import scanorama
import scvi
import torch
from rich import print
from scib_metrics.benchmark import Benchmarker

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

        self.device = "cuda" if torch.cuda.is_available() else "cpu"
        st.title("Integrate datasets")
        col1, col2, col3 = st.columns(3, gap="medium")
        with col1:
            self.auto_integrate_recipies()
            self.ingest()
            self.scanorama_integrate()
        with col2:
            self.quick_map()
            self.concat()
        with col3:
            self.bbknn()
            self.umap()

        st.header("Train deep learning model")
        int_col1, int_col2, int_col3 = st.columns(3)
        with int_col1:
            self.scvi_integrate()
        with int_col2:
            self.scanvi_integrate()
        with int_col3:
            self.scvi_integrate_graphs()
            


    def save_adata(self):
        sc.write(filename=os.path.join(os.getenv('WORKDIR'), 'adata', st.session_state.adata_state.current.adata_name), adata=self.adata)
        st.session_state.adata_state.current.adata = self.adata
            
            
    def auto_integrate_recipies(self):
        with st.form(key="auto_integrate_recipies_form"):
            st.subheader("Autointegrate recipies")
            recipie = st.selectbox(label="Recipie", options=['rec'])
            subcol1, _, _ = st.columns(3)
            submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True, disabled=(not st.session_state.sync_genes))
            if submit_btn:
                with st.spinner(text="Applying integration recipie"):
                    self.save_adata()
            
        
    def ingest(self):
        with st.form(key="ingest_form"):
            st.subheader("Integrate with Ingest", help="Integrates embeddings and annotations of an adata with a reference dataset adata_ref through projecting on a PCA (or alternate model) that has been fitted on the reference data. The function uses a knn classifier for mapping labels and the UMAP package [McInnes18] for mapping the embeddings.")
            st.markdown(f"""<div style='color: rgb(50, 168, 82); font-weight: bold''>{st.session_state.adata_ref.adata_name} → {st.session_state.adata_target.adata_name}</div>""", unsafe_allow_html=True)
            obs = st.multiselect(label="Obs", options=st.session_state.adata_ref.adata.obs_keys())
            subcol1, _, _ = st.columns(3)
            submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True, disabled=(not st.session_state.sync_genes))
            if submit_btn:
                try:
                    with st.spinner(text="Integrating datasets"):
                        sc.pp.pca(st.session_state.adata_ref.adata)
                        sc.pp.neighbors(st.session_state.adata_ref.adata)
                        sc.tl.umap(st.session_state.adata_ref.adata)
                        sc.tl.ingest(adata=st.session_state.adata_target.adata, adata_ref=st.session_state.adata_ref.adata, obs=obs)
                        self.save_adata()
                        st.toast("Integration complete", icon="✅")
                except Exception as e:
                    st.toast("Failed to integrate datasets", icon="❌")
                    st.error(e)
                    print("Error: ", e)
                    
                    
    def scanorama_integrate(self):
        with st.form(key="scanorama_integrate_form"):
            st.subheader("Integrate with Scanorama")
            subcol1, _, _ = st.columns(3)
            submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True, disabled=(not st.session_state.sync_genes))
            if submit_btn:
                with st.spinner(text="Integrating with Scanorama"):
                    scanorama.integrate_scanpy(adatas, dimred=50)
                    self.save_adata()
            
            
    def bbknn(self):
        with st.form(key="bbknn_form"):
            st.subheader("BBKNN", help="Batch balanced kNN alters the kNN procedure to identify each cell’s top neighbours in each batch separately instead of the entire cell pool with no accounting for batch. The nearest neighbours for each batch are then merged to create a final list of neighbours for the cell. Aligns batches in a quick and lightweight manner.")
            st.write(f"Apply to {st.session_state.adata_state.current.adata_name}")
            batch_key = st.selectbox(label="Batch key", options=st.session_state.adata_state.current.adata.obs_keys())
            subcol1, _, _ = st.columns(3)
            submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
            if submit_btn:
                with st.spinner(text="Computing PCA"):
                    sc.tl.pca(st.session_state.adata_state.current.adata)
                with st.spinner(text="Applying BBKNN"): 
                    sc.external.pp.bbknn(st.session_state.adata_state.current.adata, batch_key=batch_key)

                self.save_adata()
                    
            
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

            
            subcol1_btn, _, _ = st.columns(3)
            submit_btn = subcol1_btn.form_submit_button(label="Run", use_container_width=True, disabled=(not st.session_state.sync_genes))
            
            if submit_btn:
                try:
                    source_adata = st.session_state.adata_state.load_adata(st.session_state.current_workspace.id, source_dataset_name).adata
                    dest_adata = st.session_state.adata_state.load_adata(st.session_state.current_workspace.id, dest_dataset_name).adata
                    
                    st.session_state.adata_state.insert_record(AdataModel(
                        work_id=st.session_state.current_workspace.id, adata=dest_adata, 
                        filename=os.path.join(os.getenv('WORKDIR'), "adata", f"{dest_dataset_name}.h5ad"), adata_name=dest_dataset_name)
                    )

                    self.save_adata()
                    
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
            subcol1_btn, _, _ = st.columns(3)
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
                            self.save_adata()
                            st.toast("Successfully concatenated dataframes", icon="✅")
                    except Exception as e:
                        st.toast("Couldn't concatenate dataframes", icon="❌")
                        empty.error(e)
                        
                    
            
    def umap(self):
        
        def compute_umap():
            with st.spinner(text="Computing umap"):
                adata = st.session_state.adata_state.current.adata
                sc.pp.neighbors(adata)
                sc.tl.leiden(adata)
                sc.tl.umap(adata)
                for color in colors:
                    df_umap = pd.DataFrame({'umap1': adata.obsm['X_umap'][:,0], 'umap2': adata.obsm['X_umap'][:,1], 'color': adata.obs[color]})
                    container.scatter_chart(data=df_umap, x='umap1', y='umap2', color='color', size=10)
                    
        with st.form(key="integrate_umap_form"):
            st.subheader("UMAP")
            default = 'leiden'
            options = [st.session_state.adata_state.current.adata.obs_keys(), 'leiden']
            for i, item in enumerate(st.session_state.adata_state.current.adata.obs_keys()):
                if item.lower().replace("_", "").__contains__(default):
                    options = st.session_state.adata_state.current.adata.obs_keys()
            colors = st.multiselect(label="Color (obs)", options=options, default=default)
            container = st.container()
            subcol1, _, _ = st.columns(3)
            submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
            
            compute_umap()
            
            if submit_btn:
                compute_umap()

    
    def scvi_integrate(self):
        with st.form(key="scvi_integrate_form"):
            st.subheader("Integrate with SCVI", help="Train a deep learning model to integrate a dataset with a batch key. Cell type annotations are not required.")

            batch_key = st.selectbox(label="Batch key", options=st.session_state.adata_state.current.adata.obs_keys())

            layers = []
            for layer in st.session_state.adata_state.current.adata.layers:
                layers.append(str(layer))
            layer = st.selectbox(label="Layer", options=layers)

            st.write("Model params")
            input_col1, input_col2, input_col3 = st.columns(3)
            n_layers = input_col1.number_input(label="n_layers", min_value=1, step=1, format="%i", value=2)
            n_latent = input_col2.number_input(label="n_latent", min_value=1, step=1, format="%i", value=30)
            n_hidden = input_col3.number_input(label="n_hidden", min_value=1, step=1, format="%i", value=128)

            subcol1, _, _ = st.columns(3)
            submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)

            if submit_btn:
                with st.spinner(text="Training model"):
                    scvi.settings.seed = 0

                    #define and train model

                    scvi.model.SCVI.setup_anndata(st.session_state.adata_state.current.adata, layer=layer, batch_key=batch_key)
                    model = scvi.model.SCVI(st.session_state.adata_state.current.adata, n_layers=n_layers, n_latent=n_latent, n_hidden=n_hidden, gene_likelihood="nb")
                    model.train(use_gpu = (self.device == "cuda"))  

                    #evaluate latent representation

                    SCVI_LATENT_KEY = "X_scVI"
                    st.session_state.adata_state.current.adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()

                    sc.pp.neighbors(st.session_state.adata_state.current.adata, use_rep=SCVI_LATENT_KEY)
                    sc.tl.leiden(st.session_state.adata_state.current.adata)

                    SCVI_MDE_KEY = "X_scVI_MDE"
                    st.session_state.adata_state.current.adata.obsm[SCVI_MDE_KEY] = scvi.model.utils.mde(st.session_state.adata_state.current.adata.obsm[SCVI_LATENT_KEY])
                

    def scanvi_integrate(self):
        with st.form(key="scanvi_integrate_form"):
            st.subheader("Integrate with scANVI", help="Integrates a dataset with scANVI using a deep learning model. Requires cell type annotations.")
            input_col1, input_col2 = st.columns(2)
            batch_key = input_col1.selectbox(label="Batch key", options=st.session_state.adata_state.current.adata.obs_keys())
            labels_key = input_col2.selectbox(label="Labels key", options=st.session_state.adata_state.current.adata.obs_keys())
            unlabeled_category = input_col1.text_input(label="Unlabeled category", value="unknown", help="Value used for unlabeled cells in labels_key used to setup AnnData with scvi.")

            st.write("Model params")
            input_col1, input_col2, input_col3 = st.columns(3)
            n_layers = input_col1.number_input(label="n_layers", min_value=1, step=1, format="%i", value=2)
            n_latent = input_col2.number_input(label="n_latent", min_value=1, step=1, format="%i", value=30)
            n_hidden = input_col3.number_input(label="n_hidden", min_value=1, step=1, format="%i", value=128)

            subcol1, _, _ = st.columns(3)
            submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
            if submit_btn:

                scvi.model.SCVI.setup_anndata(st.session_state.adata_state.current.adata, layer="counts", batch_key=batch_key)
                model = scvi.model.SCVI(st.session_state.adata_state.current.adata, n_layers=n_layers, n_latent=n_latent, n_hidden=n_hidden, gene_likelihood="nb")

                scanvi_model = scvi.model.SCANVI.from_scvi_model(
                    model,
                    adata=st.session_state.adata_state.current.adata,
                    labels_key=labels_key,
                    unlabeled_category=unlabeled_category,
                )

                scanvi_model.train(max_epochs=20, n_samples_per_label=100)

                SCANVI_LATENT_KEY = "X_scANVI"
                st.session_state.adata_state.current.adata.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation(st.session_state.adata_state.current.adata)

                SCANVI_MDE_KEY = "X_scANVI_MDE"
                st.session_state.adata_state.current.adata.obsm[SCANVI_MDE_KEY] = scvi.model.utils.mde(st.session_state.adata_state.current.adata.obsm[SCANVI_LATENT_KEY])

                
    def scvi_integrate_graphs(self):
        with st.form(key="scvi_integrate_graphs_form"):
            st.subheader("SCVI integration plots")
            is_embeddings = ("X_scVI" in st.session_state.adata_state.current.adata.obsm_keys())
            input_col1, input_col2 = st.columns(2)
            colors = input_col1.multiselect(label="Colors", options=st.session_state.adata_state.current.adata.obs_keys(), disabled=(not is_embeddings))
            options = []
            for obsm in st.session_state.adata_state.current.adata.obsm_keys():
                if obsm == "X_scVI":
                    options.append("scVI")
                if obsm == "X_scANVI":
                    options.append("scANVI")
            embedding = input_col2.selectbox(label="Embedding", options=options, disabled=(not is_embeddings))
            preserve_neighbours = st.toggle(label="Preserve neighbours", value=True, disabled=(not is_embeddings))
            subcol1, _, _ = st.columns(3)
            submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True, disabled=(not is_embeddings))
            if submit_btn:
                with st.spinner(text="Generating plots"):
                    for color in colors:
                        st.markdown(f"""<div style='display: flex; align-items: center; justify-content: center;'><h1 style='text-align: center; font-size: 1.6rem;'>{color}</h1></div>""", unsafe_allow_html=True)

                        if embedding == "scVI":
                            if preserve_neighbours:
                                df = pd.DataFrame({'umap1': st.session_state.adata_state.current.adata.obsm['X_scVI_MDE'][:,0], 'umap2': st.session_state.adata_state.current.adata.obsm['X_scVI_MDE'][:,1], 'color': st.session_state.adata_state.current.adata.obs[color]})
                                st.scatter_chart(df, x='umap1', y='umap2', color='color', size=10)
                            else:
                                df = pd.DataFrame({'umap1': st.session_state.adata_state.current.adata.obsm['X_scVI'][:,0], 'umap2': st.session_state.adata_state.current.adata.obsm['X_scVI'][:,1], 'color': st.session_state.adata_state.current.adata.obs[color]})
                                st.scatter_chart(df, x='umap1', y='umap2', color='color', size=10)
                        elif embedding == "scANVI":
                            if preserve_neighbours:
                                df = pd.DataFrame({'umap1': st.session_state.adata_state.current.adata.obsm['X_scANVI_MDE'][:,0], 'umap2': st.session_state.adata_state.current.adata.obsm['X_scANVI_MDE'][:,1], 'color': st.session_state.adata_state.current.adata.obs[color]})
                                st.scatter_chart(df, x='umap1', y='umap2', color='color', size=10)
                            else:
                                df = pd.DataFrame({'umap1': st.session_state.adata_state.current.adata.obsm['X_scANVI'][:,0], 'umap2': st.session_state.adata_state.current.adata.obsm['X_scANVI'][:,1], 'color': st.session_state.adata_state.current.adata.obs[color]})
                                st.scatter_chart(df, x='umap1', y='umap2', color='color', size=10)

            
            
    


sidebar = Sidebar()

sidebar.show(integrate=True)

sidebar.show_preview(integrate=True)
        
integrate = Integrate()

sidebar.export_script()

sidebar.delete_experiment_btn()

sidebar.show_version()