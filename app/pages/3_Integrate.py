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
from scvi.model._utils import get_max_epochs_heuristic

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
    """
    Integrate multiple datasets into a single dataset with tools to minimise batch effect.

    Notes
    -----
    .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/integrate_page.png
    """
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

        self.scvi_metrics()
        
            


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
        """
        Integrates embeddings and annotations of an adata with a reference dataset adata_ref through projecting on a PCA (or alternate model) that has been fitted on the reference data. The function uses a knn classifier for mapping labels and the UMAP package [McInnes18]_ for mapping the embeddings.

        Parameters
        ----------
        obs: List[str]
            List of obs keys in adata_ref.obs which need to be mapped to adata.obs.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/ingest.png

        Example
        -------
        import scanpy as sc

        sc.pp.neighbors(adata_ref.adata)
        sc.tl.umap(adata_ref.adata)
        sc.tl.ingest(adata=adata_target, adata_ref=adata_ref.adata, obs='umap')
        """
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
        """
        Integrate spatial datasets with scanorama.

        Example
        -------
        !pip install scanorama # install external libraries
        import scanpy as sc
        """
        with st.form(key="scanorama_integrate_form"):
            st.subheader("Integrate spatial data with Scanorama")
            col1, col2 = st.columns(2)
            label = col1.text_input(label="Label", value="library_id")
            index_unique = col2.text_input(label="index_unique", value="-")
            index = 0
            disabled=(not st.session_state.sync_genes)
            if 'uns' in st.session_state.adata_ref and 'uns' in st.session_state.adata_target:
                for i, key in enumerate(st.session_state.adata_ref.uns_keys()):
                    if key.lower() == "spatial":
                        index = i
                spatial_key1 = col1.selectbox(label="Spatial key 1", options=st.session_state.adata_ref.uns_keys(), index=index)

                for i, key in enumerate(st.session_state.adata_target.uns_keys()):
                    if key.lower() == "spatial":
                        index = i
                spatial_key2 = col2.selectbox(label="Spatial key 2", options=st.session_state.adata_target.uns_keys(), index=index)
            else:
                st.info("No 'uns' key found in datasets")
                disabled = True

            subcol1, _, _ = st.columns(3)
            submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True, disabled=disabled)
            if submit_btn:
                with st.spinner(text="Integrating with Scanorama"):
                    adatas = [st.session_state.adata_ref.adata, st.session_state.adata_target.adata]

                    adatas_cor = scanorama.correct_scanpy(adatas, return_dimred=True)

                    adata_spatial = sc.concat(
                        adatas_cor,
                        label=label,
                        uns_merge="unique",
                        keys=[
                            k
                            for d in [
                                adatas_cor[0].uns[spatial_key1],
                                adatas_cor[1].uns[spatial_key2],
                            ]
                            for k, v in d.items()
                        ],
                        index_unique=index_unique,
                    )
                    
                    #self.save_adata() TODO: find a way to save adata_spatial, how do we display to user?
            
            
    def bbknn(self):
        """
        Batch balanced kNN [Polanski19] alters the kNN procedure to identify each cell's top neighbours in each batch separately instead of the entire cell pool with no accounting for batch. The nearest neighbours for each batch are then merged to create a final list of neighbours for the cell. Aligns batches in a quick and lightweight manner.
        
        Parameters
        ----------
        batch_key: str
            Provide a batch key for BBKNN.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/bbknn.png

        Example
        -------
        !pip install bbknn # install external libraries
        import scanpy as sc

        sc.tl.pca(adata)
        sc.external.pp.bbknn(adata, batch_key='BATCH')
        """
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
        """
        Map an attribute of one dataset onto another. 

        Parameters
        ----------
        source_adata: str
            Source dataset for the attributes.

        source_attributes: str
            Attributes to map (e.g. uns.louvain_colors).

        target_adata: str
            Target dataset for mapping the attributes (the original attributes will be replaced). 

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/quick_map.png

        Example
        -------
        target_adata.obs = source_adata.obs
        """
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
        """
        Concatenate dataframes with batches.

        Parameters
        ----------
        adata1: str
            First dataset.

        adata2: str
            Second dataset.

        batch_category1: str
            First batch annotation.

        batch_category2: str
            Second batch annotation.

        batch_key: str
            Add the batch annotation to obs using this key.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/concat.png

        Example
        -------
        adata_concat = adata1.concatenate(adata2, batch_key="BATCH", batch_categories=['1', '2'])
        """
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
        """
        Compute UMAP coordinates as part of the integration process.

        Parameters
        ----------
        color: List[str]
            Color of UMAP clusters. If multiple colors are supplied, it will generate separate plots for each color.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/umap_integrate.png

        Example
        -------
        import scanpy as sc

        sc.pp.neighbors(adata)
        sc.tl.leiden(adata)
        sc.tl.umap(adata)
        sc.pl.umap(adata, color='BATCH')
        """
        
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
        """
        Train a deep learning model using SCVI tools to integrate a dataset with a batch key. Cell type annotations are not required.

        Parameters
        ----------
        batch_key: str
            Provide a batch key.

        layer: str
            Provide an optional layer key.

        n_layers: int
            Number of layers in the neural network.

        n_latent: int
            Dimensionality of the latent space in the neural network.

        n_hidden: int
            Number of nodes per hidden layer.

        batch_size: int
            Training batch size.

        max_epochs: int
            Number of epochs to train the neural network.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/scvi_integrate.png

        Example
        -------
        import scanpy as sc

        batch_key = "BATCH"

        scvi.settings.seed = 42
        scvi.model.SCVI.setup_anndata(adata, batch_key=batch_key)
        model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, n_hidden=128, gene_likelihood="nb")
        model.train(use_gpu = True, batch_size=128, max_epochs=400)  

        #evaluate latent representation
        SCVI_LATENT_KEY = "X_scVI"
        adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()

        sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY)
        sc.tl.leiden(adata)

        SCVI_MDE_KEY = "X_scVI_MDE"
        adata.obsm[SCVI_MDE_KEY] = scvi.model.utils.mde(adata.obsm[SCVI_LATENT_KEY])    
        """
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
            input_col1, input_col2 = st.columns(2)
            batch_size = input_col1.number_input(label="batch_size", min_value=1, step=8, value=128, format="%i")
            max_epochs = input_col2.number_input(label="max_epochs", min_value=1, step=5, value=get_max_epochs_heuristic(st.session_state.adata_state.current.adata.n_obs), format="%i")

            subcol1, _, _ = st.columns(3)
            submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)

            if submit_btn:
                with st.spinner(text="Training model"):
                    scvi.settings.seed = 42

                    #define and train model

                    scvi.model.SCVI.setup_anndata(st.session_state.adata_state.current.adata, layer=layer, batch_key=batch_key)
                    model = scvi.model.SCVI(st.session_state.adata_state.current.adata, n_layers=n_layers, n_latent=n_latent, n_hidden=n_hidden, gene_likelihood="nb")
                    model.train(use_gpu = (self.device == "cuda"), batch_size=batch_size, max_epochs=max_epochs)  

                    #evaluate latent representation

                    SCVI_LATENT_KEY = "X_scVI"
                    st.session_state.adata_state.current.adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()

                    sc.pp.neighbors(st.session_state.adata_state.current.adata, use_rep=SCVI_LATENT_KEY)
                    sc.tl.leiden(st.session_state.adata_state.current.adata)

                    SCVI_MDE_KEY = "X_scVI_MDE"
                    st.session_state.adata_state.current.adata.obsm[SCVI_MDE_KEY] = scvi.model.utils.mde(st.session_state.adata_state.current.adata.obsm[SCVI_LATENT_KEY])
                

    def scanvi_integrate(self):
        """
        Integrates a dataset with scANVI using a deep learning model. Requires cell type annotations.

        Parameters
        ----------
        batch_key: str
            Provide a batch key.

        labels_key: str
            Provide a labels key for known cell type annotations.

        unlabeled_category: str
            Value used for unlabeled cells in labels_key used to setup AnnData with scvi.

        n_samples_per_label: int
            Number of subsamples for each label class to sample per epoch. By default, there is no label subsampling.

        n_layers: int
            Number of layers in the neural network.

        n_latent: int
            Dimensionality of the latent space in the neural network.

        n_hidden: int
            Number of nodes per hidden layer.

        max_epochs: int
            Number of epochs to train the neural network.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/scanvi_integrate.png

        Example
        -------
        !pip install scvi-tools # Install external libraries
        import scvi
        
        scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")
        model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, n_hidden=128, gene_likelihood="nb")

        scanvi_model = scvi.model.SCANVI.from_scvi_model(model, adata=adata, labels_key="cell_type", unlabeled_category="unknown")

        scanvi_model.train(use_gpu = True, max_epochs=400, n_samples_per_label=100)

        SCANVI_LATENT_KEY = "X_scANVI"
        adata.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation(adata)

        SCANVI_MDE_KEY = "X_scANVI_MDE"
        adata.obsm[SCANVI_MDE_KEY] = scvi.model.utils.mde(adata.obsm[SCANVI_LATENT_KEY])

        """
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
            input_col1, input_col2 = st.columns(2)
            max_epochs = input_col1.number_input(label="max epochs", value=get_max_epochs_heuristic(st.session_state.adata_state.current.adata.n_obs), format="%i", step=1)
            n_samples_per_label = input_col2.number_input(label="n_samples_per_label", min_value=1, value=100, format="%i", step=1)

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

                scanvi_model.train(use_gpu = (self.device == "cuda"), max_epochs=max_epochs, n_samples_per_label=n_samples_per_label)

                SCANVI_LATENT_KEY = "X_scANVI"
                st.session_state.adata_state.current.adata.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation(st.session_state.adata_state.current.adata)

                SCANVI_MDE_KEY = "X_scANVI_MDE"
                st.session_state.adata_state.current.adata.obsm[SCANVI_MDE_KEY] = scvi.model.utils.mde(st.session_state.adata_state.current.adata.obsm[SCANVI_LATENT_KEY])

                
    def scvi_integrate_graphs(self):
        """
        Plot integrated dataset with scVI or scANVI embeddings.

        Parameters
        ----------
        colors: str
            Key for colour.

        embedding: str
            Embedding to choose from.

        preserve_neighbours: bool
            Uses minimum distortion embedding (MDE) to minimally distort relationships between pairs of items, while possibly satisfying some constraints.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/scvi_integrate_plot.png

        Example
        -------
        import scvi
        import matplotlib.pyplot as plt

        # .. Do integration with scvi or scanvi

        # in this example we use scvi embeddings

        scvi_df = pd.DataFrame({'umap1': adata.obsm['X_scVI'][:,0], 'umap2': adata.obsm['X_scVI'][:,1], 'color': adata.obs[color]})
        scvi_df_mde = pd.DataFrame({'umap1': adata.obsm['X_scVI_MDE'][:,0], 'umap2': adata.obsm['X_scVI_MDE'][:,1], 'color': adata.obs[color]})
        
        fig, (ax1, ax2) = plt.subplots(1, 2)
        ax1.scatter(x=scvi_df['umap1'], y=scvi_df['umap2'], color=scvi_df['color'])
        ax1.title('scvi embedding')
        ax1.xlabel('umap1')
        ax1.ylabel('umap2')
        ax2.scatter(x=scvi_df_mde['umap1'], y=scvi_df_mde['umap2'], color=scvi_df_mde['color'])
        ax2.title('scvi embedding with MDE')
        ax2.xlabel('umap1')
        ax2.ylabel('umap2')

        """
        with st.form(key="scvi_integrate_graphs_form"):
            st.subheader("SCVI integration plots")
            is_embeddings = ("X_scVI" in st.session_state.adata_state.current.adata.obsm_keys())
            if not is_embeddings:
                st.info("No embeddings found for scVI.")
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

    def scvi_metrics(self):
        """
        use the scib-metrics package to assess the quality of the integration.

        Parameters
        ----------
        batch_key: str
            Use a batch key.

        label_key: str
            Use a label key.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/scvi_integrate_metrics.png

        Example
        -------
        import scanpy as sc
        from scib_metrics.benchmark import Benchmarker

        SCANVI_LATENT_KEY = "X_scANVI"
        SCVI_LATENT_KEY = "X_scVI"

        latent_keys = []
        if SCVI_LATENT_KEY in adata.obsm_keys():
            latent_keys.append(SCVI_LATENT_KEY) 
            if SCANVI_LATENT_KEY in adata.obsm_keys():
                latent_keys.append(SCANVI_LATENT_KEY)
            sc.tl.pca(adata)
            latent_keys.append("X_pca")

            bm = Benchmarker(
                adata,
                batch_key="_scvi_batch",
                label_key="_scvi_batch",
                embedding_obsm_keys=latent_keys,
                n_jobs=-1,
            )
                    
            bm.benchmark()
            df = bm.get_results(min_max_scale=False)
            print(df)
        """
        with st.form(key="scvi_metrics_form"):
            st.subheader("scVI integration metrics")
            col1, col2, col3, _, _, _ = st.columns(6)
            is_embeddings = ("X_scVI" in st.session_state.adata_state.current.adata.obsm_keys())
            if not is_embeddings:
                info_col1, _, _ = st.columns(3)
                info_col1.info("No embeddings found for scVI.")
            batch_key = col1.selectbox(label="batch_key", options=st.session_state.adata_state.current.adata.obs_keys(), disabled=(not is_embeddings))
            label_key = col2.selectbox(label="Label key", options=st.session_state.adata_state.current.adata.obs_keys(), disabled=(not is_embeddings))
            subcol1, _, _, _, _, _, _, _, _ = st.columns(9)
            submit_btn = subcol1.form_submit_button(label="Compute", use_container_width=True, disabled=(not is_embeddings))
            if submit_btn:
                with st.spinner(text="Computing metrics"):

                    SCANVI_LATENT_KEY = "X_scANVI"
                    SCVI_LATENT_KEY = "X_scVI"

                    latent_keys = []
                    if SCVI_LATENT_KEY in st.session_state.adata_state.current.adata.obsm_keys():
                        latent_keys.append(SCVI_LATENT_KEY) 
                    if SCANVI_LATENT_KEY in st.session_state.adata_state.current.adata.obsm_keys():
                        latent_keys.append(SCANVI_LATENT_KEY)
                    sc.tl.pca(st.session_state.adata_state.current.adata)
                    latent_keys.append("X_pca")

                    bm = Benchmarker(
                        st.session_state.adata_state.current.adata,
                        batch_key=batch_key,
                        label_key=label_key,
                        embedding_obsm_keys=latent_keys,
                        n_jobs=-1,
                    )
                    
                    bm.benchmark()
                    df = bm.get_results(min_max_scale=False)
                    st.dataframe(df)
            


sidebar = Sidebar()

sidebar.show(integrate=True)

sidebar.show_preview(integrate=True)
        
integrate = Integrate()

sidebar.export_script()

sidebar.delete_experiment_btn()

sidebar.show_version()