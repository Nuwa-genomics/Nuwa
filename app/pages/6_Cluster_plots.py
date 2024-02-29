import numpy as np
import matplotlib.pyplot as plt
from torch.utils.data import DataLoader

from ml.citeseq.dataset import TabularDataset
from ml.citeseq.train import get_encodings
from ml.citeseq.model import CiteAutoencoder
from ml.solo_scvi.solo_model import solo_model
from ml.DeepST.deepst.main import *
import plotly.express as px
from ml.solo_scvi.utils import *

import pandas as pd

import scanpy as sc
import streamlit as st
import umap.umap_ as umap
from utils.plotting import get_color_embeddings_from_key

from components.sidebar import *
from models.AdataModel import AdataModel
from database.schemas import schemas
from state.AdataState import AdataState
from database.database import SessionLocal
import os
from state.StateManager import StateManager

st.set_page_config(layout="wide", page_title='Nuwa', page_icon='🧬')

os.chdir('/app')

with open('css/cluster.css') as f:
    st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)

with open('css/common.css') as f:
    st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)

common_style = """
    <style>
    footer {visibility: hidden;}
    .st-emotion-cache-1cypcdb {background: linear-gradient(180deg, rgb(5, 39, 103) 0%, #3a0647 70%); box-shadow: 1px 0 10px -2px #000;}
    </style>
"""
st.markdown(common_style, unsafe_allow_html=True)


class Cluster_plots:
    """
    Run clustering using a variety of dimensionality reduction algorithms and deep learning models.

    Notes
    -----
    .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/cluster_plots_page.png
    """
    def __init__(self, adata):

        st.title("Cluster plots")

        self.col1, self.col2, self.col3 = st.columns(3)

        self.adata = adata

        self.conn: SessionLocal = SessionLocal()

        self.PLOT_HEIGHT = 800
        self.MARKER_SIZE = 32

        if "cluster_plots" not in st.session_state:
            st.session_state["cluster_plots"] = dict(pca=None, tsne=None, autoencoder=None, umap=None, variance_ratio=None)


    def save_adata(self, name):
        try:
            new_adata = schemas.Adata(
                work_id=int(st.session_state.current_workspace.id),
                adata_name=f"{name}",
                filename=os.path.join(os.getenv('WORKDIR'), "adata", f"{name}.h5ad"),
                notes=st.session_state.adata_state.current.notes
            )
            self.conn.add(new_adata)
            self.conn.commit()
            self.conn.refresh(new_adata)
        except Exception as e:
            print(e)

    def _autoencoder_cluster_plot(self, colors):

        trained_model = st.session_state["trained_model"]
        if isinstance(trained_model, CiteAutoencoder):
            device = st.session_state["device"]
            
            if st.session_state["cluster_plots"]["autoencoder"] == None:
                sc.pp.neighbors(self.adata, n_neighbors=10, n_pcs=40)
                sc.tl.leiden(self.adata)

            rna_df_original = self.adata.to_df()
            rna_df = rna_df_original.reset_index(drop=True)

            test_ds = TabularDataset(rna_df.to_numpy(dtype=np.float32))
            test_dl = DataLoader(test_ds, batch_size=64, shuffle=False)

            encodings = get_encodings(trained_model, test_dl, device)
            encodings = encodings.cpu().numpy()
            metadata_df = self.adata.obs
            cell_ids = rna_df_original.index.values
            plot_df = metadata_df.loc[metadata_df.index.isin(cell_ids)]
            embedding2d = umap.UMAP(random_state=0).fit_transform(encodings)
            embedding3d = umap.UMAP(random_state=0, n_components=3).fit_transform(encodings)

            st.session_state.adata_state.current.adata.obsm["X_citeseq_2d"] = embedding2d
            st.session_state.adata_state.current.adata.obsm["X_citeseq_3d"] = embedding3d
                                        
            plot_df["UMAP1"] = embedding2d[:, 0]
            plot_df["UMAP2"] = embedding2d[:, 1]
                                    
            st.session_state["cluster_plots"]["autoencoder"] = dict(df=plot_df, x="UMAP1", y="UMAP2", color_keys=np.array([colors]).flatten(), size=self.MARKER_SIZE)

    def autoencoder_cluster_plot(self):
        """
        Scatter plot showing embeddings for trained model.

        Parameters
        ----------
        color: str
            Colour based on an obs value (available for citeseq autoencoder).

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/autoencoder_citeseq.png
        .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/autoencoder_solo.png

        Example
        -------
        # For Citeseq see https://nuwa-genomics.github.io/Nuwa/reference/Create_CiteSeq_model/init_model.html#python-example

        # For Solo see https://nuwa-genomics.github.io/Nuwa/reference/Create_Solo_model/init_model.html#python-example
        """
        with self.col1:
            try:
                if "trained_model" not in st.session_state:
                    st.subheader("Autoencoder")
                    st.divider()
                    st.warning("Model hasn't been trained to provide data")
                else:
                    trained_model = st.session_state["trained_model"]
                    

                    if isinstance(trained_model, CiteAutoencoder):
                        with st.form(key="citeseq_form"):
                            st.subheader("Autoencoder Clusters")
                            
                            options = np.append(self.adata.obs_keys(), self.adata.var_names)
                            colors = st.multiselect(label="Colour", options=options, default='leiden')

                            if "losses" in st.session_state:
                                losses = st.session_state["losses"]
                                training_loss = losses["train"]
                                valid_loss = losses["valid"]
                                losses_df = pd.DataFrame({'Train': training_loss, 'Valid': valid_loss})
                                st.markdown("""<p style='font-size: 16px; text-align: center; font-weight: bold;'>Training/validation losses</p>""", unsafe_allow_html=True)
                                st.line_chart(losses_df, use_container_width=True, height=290)
                                    
                            subcol1, _, _ = st.columns(3)
                            submit_btn = subcol1.form_submit_button(label="Update", use_container_width=True)
                                        
                            if submit_btn:
                                df = st.session_state["cluster_plots"]["autoencoder"]["df"]
                                st.session_state["cluster_plots"]["autoencoder"] = dict(df=df, x="UMAP1", y="UMAP2", color_keys=np.array([colors]).flatten(), size=self.MARKER_SIZE)
                                st.toast("Upadated colours", icon="✅")

                    elif(isinstance(trained_model, solo_model)):
                        with st.form(key="solo_doublet_form"):
                            st.subheader("Doublet prediction")
                            
                            with st.spinner("Loading Solo doublet predictions"):
                                
                                df_solo = pd.DataFrame({'UMAP1': self.adata.obsm['X_umap'][:,0], 'UMAP2': self.adata.obsm['X_umap'][:,1]})
                            
                                st.session_state["cluster_plots"]["autoencoder"] = dict(df=df_solo, x="UMAP1", y="UMAP2", color_keys=['prediction'], size=self.MARKER_SIZE)
                                
                                solo_df, vae_df = get_solo_model_history(solo=trained_model.solo, vae=trained_model.vae)
                                solo, vae = st.tabs(['Solo', 'Vae'])
                                with solo:
                                    st.line_chart(solo_df, use_container_width=True, height=360, color=['#52f27d', '#f25272'])
                                with vae:
                                    st.line_chart(vae_df, use_container_width=True, height=360, color=['#52f27d', '#f25272'])

                                subcol1, _, _ = st.columns(3)
                                submit_btn = subcol1.form_submit_button(label="Filter doublets", use_container_width=True)
                                
                                if submit_btn:
                                    self.adata = self.adata[self.adata.obs.prediction == 'singlet']
                                    self.save_adata(name="adata_solo")
                                    st.toast("Removed doublet predictions", icon='✅')


                    else:
                        st.error("Unknown model")

            except Exception as e:
                st.error(e)


    def _do_pca(self, zero_center = True):
        sc.tl.pca(self.adata, zero_center=zero_center, svd_solver='arpack', random_state=42)
        use_raw = True
        if self.adata.raw:
            for var in self.adata.var_names:
                if not (self.adata.raw.var_names.__contains__(var)):
                    use_raw = False
                    st.info("Gene ID format has been changed, setting 'use_raw' to false.")
                    # TODO: decide how to handle use_raw for streamlit plots
                    break


    def pca_graph(self):
        """
        Compute PCA coordinates with optional colour.

        Parameters
        ----------
        gene: List[str]
            Genes to supply as color argument in PCA.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/pca_cluster.png

        Example
        -------
        import scanpy as sc

        genes = ['brca1', 'brca2']
        sc.tl.pca(adata, svd_solver='arpack')
        sc.pl.pca(adata, color=genes)
        """
        with self.col3:
            try:

                with st.form(key="pca_cluster_form"):
                    st.subheader("PCA")
                    pca_colors = st.multiselect(label="Color", options=np.append(self.adata.obs_keys(), self.adata.var_names), default=(self.adata.var_names[0]), key="ms_pca_gene", max_selections=24)

                    zero_center = st.toggle(label="zero center", value=True)
                    df_pca = pd.DataFrame({'PCA1': self.adata.obsm['X_pca'][:,0], 'PCA2': self.adata.obsm['X_pca'][:,1]})  
                    st.session_state["cluster_plots"]["pca"] = dict(df=df_pca, color_keys=np.array([pca_colors]).flatten(), x="PCA1", y="PCA2", size=self.MARKER_SIZE)

                    subcol1, _, _ = st.columns(3)
                    submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
                    if submit_btn:
                        with st.spinner("Computing PCA coordinates"):
                            self._do_pca(zero_center=zero_center)
                            df_pca = pd.DataFrame({'PCA1': self.adata.obsm['X_pca'][:,0], 'PCA2': self.adata.obsm['X_pca'][:,1]})  
                            st.session_state["cluster_plots"]["pca"] = dict(df=df_pca, color_keys=np.array([pca_colors]).flatten(), x="PCA1", y="PCA2", size=self.MARKER_SIZE)
                            st.toast("Recomputed PCA", icon="✅")


            except Exception as e:
                st.error(e)


    def _do_var_ratio(self, n_pcs = 30, log = False):

        sc.pl.pca_variance_ratio(adata, log=log, n_pcs=n_pcs)

        pcs = []
        pc_val = []

        for i, pc in enumerate(adata.uns['pca']['variance_ratio']):
            pcs.append(i+1)
            pc_val.append(pc)

        pcs = np.array(pcs)
        pc_val = np.array(pc_val)
        df = pd.DataFrame({'PC': pcs, 'value': pc_val})
        df["PC"] = df["PC"].astype("category")

        fig = px.scatter(df, x="PC", y="value", color="PC", height=self.PLOT_HEIGHT)
        fig.update_layout(showlegend=False)
                        
        st.session_state["cluster_plots"]["variance_ratio"] = fig


    def variance_ratio_graph(self):
        """
        Compare variance between principle components.

        Parameters
        ----------
        n_pcs: int
            Number of principle components to use.

        log: bool
            Use a log scale.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/variance_ratio.png

        Example
        -------
        import scanpy as sc
        
        sc.tl.pca(adata, svd_solver='arpack')
        sc.pl.pca_variance_ratio(adata, log=True)
        """
        with self.col2:
            try:
 
                with st.form(key="variance_ratio_form"):
                    st.subheader("PCA variance ratio")
                    subcol1, subcol2 = st.columns(2)
                    n_pcs = subcol1.number_input(label="n_pcs", min_value=1, max_value=50, value=30)
                    log = subcol1.toggle(label="Log", value=True)
                    subcol1, _, _ = st.columns(3)
                    submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)

                    if submit_btn:
                        with st.spinner("Computing variance ratio"):
                            self._do_var_ratio(n_pcs=n_pcs, log=log)
            
            except Exception as e:
                print("Error ", e)
                st.error(e)

    def _do_tSNE(self, perplexity):
        sc.pp.neighbors(self.adata, n_neighbors=10, n_pcs=40)
        sc.tl.leiden(self.adata) 
        sc.tl.tsne(self.adata, perplexity=perplexity)  


    def tsne_graph(self):
        """
        
        Parameters
        ----------
        color: str
            Colour based on an obs value.
        
        perplexity: int
            The perplexity is related to the number of nearest neighbors that is used in other manifold learning algorithms. 
            Larger datasets usually require a larger perplexity. Consider selecting a value between 5 and 50. 
            The choice is not extremely critical since t-SNE is quite insensitive to this parameter.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/tsne.png

        Example
        -------
        import scanpy as sc

        # compute neighbours to use for colours
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
        sc.tl.leiden(adata) 
        sc.tl.tsne(adata, perplexity=30)  
        sc.pl.tsne(adata, color="leiden")
        """
        with self.col2:
            try:
  
                with st.form(key="tsne_form"):

                    st.subheader("tSNE")
                    
                    options = np.append(self.adata.obs_keys(), self.adata.var_names)
                    tsne_colors = st.multiselect(label="Colour", options=options, default='leiden')
                    perplexity = st.number_input(label="Perplexity", min_value=1, value=30)

                    df_tsne = pd.DataFrame({'tSNE1': self.adata.obsm['X_tsne'][:,0], 'tSNE2': self.adata.obsm['X_tsne'][:,1]})  
                    st.session_state["cluster_plots"]["tsne"] = dict(df=df_tsne, color_keys=np.array([tsne_colors]).flatten(), x="tSNE1", y="tSNE2", size=self.MARKER_SIZE)

                    subcol1, _, _ = st.columns(3)
                    submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)

                    if submit_btn:
                        with st.spinner("Computing tSNE coordinates"):
                            # Perplexity may have changed in form so recompute tSNE before plotting
                            self._do_tsne(perplexity=perplexity)
                            df_tsne = pd.DataFrame({'tSNE1': self.adata.obsm['X_tsne'][:,0], 'tSNE2': self.adata.obsm['X_tsne'][:,1]})  
                            st.session_state["cluster_plots"]["tsne"] = dict(df=df_tsne, color_keys=np.array([tsne_colors]).flatten(), x="tSNE1", y="tSNE2", size=self.MARKER_SIZE)
                        
            except Exception as e:
                st.error(e)


    def _do_umap(self, resolution = 1, n_neighbors = 10, n_pcs = 40):
        #TODO: the use_raw may need to be set to false if gene symbols change
        sc.pp.neighbors(self.adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=42)
        sc.tl.leiden(self.adata, random_state=42, resolution=resolution) 
        sc.tl.paga(self.adata)
        sc.pl.paga(self.adata, plot=False)
        sc.tl.umap(self.adata, init_pos='paga')


    def neighbourhood_graph(self):
        """
        Compute neighbourhood graph with color mapping for selected genes.

        Parameters
        ----------
        genes: List[str]
            Genes to supply as color argument in PCA.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/neighbourhood_graph.png

        Example
        -------
        import scanpy as sc 

        genes = ['brca1', 'brca2']
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
        sc.tl.leiden(adata) 
        sc.tl.paga(adata)
        sc.pl.paga(adata, use_raw=use_raw, plot=False)
        sc.tl.umap(adata, init_pos='paga')  
        sc.pl.umap(adata)     
        """

        with self.col3:
            try:

                with st.form(key="nhood_graph_form"):
                    st.subheader("Neighbourhood graph")
                    
                    options = np.append(self.adata.obs_keys(), self.adata.var_names)
                    umap_colors = st.multiselect(label='Gene', options=options, default="leiden") 
                    col1, col2, col3 = st.columns(3)

                    resolution = col1.number_input(label="Resolution", value=1.0, step=0.1)
                    n_neighbors = col2.number_input(label="n_neighbors", value=10, step=1, min_value=1, format="%i")
                    n_pcs = col3.number_input(label="n_pcs", value=40, step=1, min_value=1, format="%i")

                    df_umap = pd.DataFrame({'UMAP1': self.adata.obsm['X_umap'][:,0], 'UMAP2': self.adata.obsm['X_umap'][:,1]})  
                    st.session_state["cluster_plots"]["umap"] = dict(df=df_umap, color_keys=np.array([umap_colors]).flatten(), x="UMAP1", y="UMAP2", size=self.MARKER_SIZE)
                        
                    subcol1, _, _ = st.columns(3)
                        
                    submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
                        
                    if submit_btn:
                        with st.spinner("Computing UMAP coordinates"):
                            self._do_umap(resolution=resolution, n_pcs=n_pcs, n_neighbors=n_neighbors)
                            df_umap = pd.DataFrame({'UMAP1': self.adata.obsm['X_umap'][:,0], 'UMAP2': self.adata.obsm['X_umap'][:,1]})  
                            st.session_state["cluster_plots"]["umap"] = dict(df=df_umap, color_keys=np.array([umap_colors]).flatten(), x="UMAP1", y="UMAP2", size=self.MARKER_SIZE)

            except Exception as e:
                st.error(e)


    def plots(self):
        st.subheader("Plots")
        autoencoder_plot, tSNE_plot, pca_plot, nhood_plot, variance_ratio_plot = st.tabs(['Autoencoder', 'tSNE', 'PCA', 'UMAP', 'Variance ratio'])
        with autoencoder_plot:
            if st.session_state["cluster_plots"]["autoencoder"] == None:
                st.info("You must train a model first.")
            else:
                params = st.session_state["cluster_plots"]["autoencoder"]
                self._plot_charts(params)

        with tSNE_plot:
            if st.session_state["cluster_plots"]["tsne"] == None:
                st.info("You must run tSNE.")
            else:
                params = st.session_state["cluster_plots"]["tsne"]
                self._plot_charts(params)

        with pca_plot:
            if st.session_state["cluster_plots"]["pca"] == None:
                st.info("You must run PCA first.")
            else:
                params = st.session_state["cluster_plots"]["pca"]
                self._plot_charts(params)

        with variance_ratio_plot:
            if st.session_state["cluster_plots"]["variance_ratio"] == None:
                st.info("You must run Variance ratio first.")
            else:
                fig = st.session_state["cluster_plots"]["variance_ratio"] 
                st.plotly_chart(fig, use_container_width=True)

        with nhood_plot:
            if st.session_state["cluster_plots"]["umap"] == None:
                st.info("You must run UMAP first.")
            else:
                params = st.session_state["cluster_plots"]["umap"]
                self._plot_charts(params)


    def _plot_charts(self, params):
        color_keys = params["color_keys"]
        for color in color_keys:
            color_embed = get_color_embeddings_from_key(key=color, adata=self.adata)
            st.markdown(f"""<p style='text-align: center; size: 16px; font-weight: bold;'>{color}</p>""", unsafe_allow_html=True)
            df = params["df"]
            df["color"] = color_embed
            st.scatter_chart(data=df, x=params['x'], y=params['y'], color='color', size=params["size"], use_container_width=True, height=self.PLOT_HEIGHT)


try:
    adata = st.session_state.adata_state.current.adata

    sidebar = Sidebar()
    sidebar.show()

    analysis = Cluster_plots(adata)

    if st.session_state["cluster_plots"]["pca"] == None:
        my_bar = st.progress(3, text="Preparing embeddings from neural network")
        analysis._autoencoder_cluster_plot(colors='leiden')
        my_bar.progress(50, "Computing tSNE")
        analysis._do_tSNE(perplexity=30)
        my_bar.progress(70, "Computing neighbourhood graph")
        analysis._do_umap()
        my_bar.progress(85, "Computing PCA")
        analysis._do_pca()
        my_bar.progress(95, "Computing variance ratio")
        analysis._do_var_ratio()
        my_bar.progress(100, "Complete")
        my_bar.empty()

    analysis.autoencoder_cluster_plot()
    analysis.neighbourhood_graph()
    analysis.tsne_graph()
    analysis.pca_graph()
    analysis.variance_ratio_graph()
    st.divider()
    analysis.plots()

    sidebar.show_preview()
    sidebar.export_script()
    sidebar.delete_experiment_btn()
    sidebar.show_version()


except Exception as e:
    if(st.session_state == {}):
        StateManager().load_session()
        st.rerun()
    else:
        st.toast(e, icon="❌")