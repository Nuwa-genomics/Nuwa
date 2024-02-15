import numpy as np
import matplotlib.pyplot as plt
from torch.utils.data import DataLoader

from ml.citeseq.dataset import TabularDataset
from ml.citeseq.train import get_encodings
from ml.citeseq.model import CiteAutoencoder
from ml.solo_scvi.solo_model import solo_model
from ml.DeepST.deepst.main import *

import pandas as pd

import scanpy as sc
import streamlit as st
import umap.umap_ as umap

from components.sidebar import *
from models.AdataModel import AdataModel
from database.schemas import schemas
from state.AdataState import AdataState
from database.database import SessionLocal
import os

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


class Cluster_analysis:
    """
    Run clustering using a variety of dimensionality reduction algorithms and deep learning models.

    Notes
    -----
    .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/cluster_analysis_page.png
    """
    def __init__(self, adata):

        st.title("Cluster analysis")

        self.col1, self.col2 = st.columns(2, gap="medium")

        self.adata = adata

        self.conn: SessionLocal = SessionLocal()


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

    def autoencoder_cluster_plot(self):
        """
        Scatter plot showing embeddings for trained model.

        Parameters
        ----------
        color: str
            Colour based on an obs value (available for citeseq autoencoder).

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/autoencoder_citeseq.png
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/autoencoder_solo.png

        Example
        -------
        # For Citeseq see https://ch1ru.github.io/Nuwa/reference/Create_CiteSeq_model/init_model.html#python-example

        # For Solo see https://ch1ru.github.io/Nuwa/reference/Create_Solo_model/init_model.html#python-example
        """
        with self.col1:
            try:
                if "trained_model" not in st.session_state:
                    st.error("Model hasn't been trained to provide data")
                else:
                    trained_model = st.session_state["trained_model"]
                    device = st.session_state["device"]

                    if isinstance(trained_model, CiteAutoencoder):
                        with st.form(key="citeseq_form"):
                            st.subheader("Autoencoder Clusters")
                            
                            with st.spinner(text="Preparing embeddings"): 
                            
                                rna_df_original = self.adata.to_df()
                                rna_df = rna_df_original.reset_index(drop=True)

                                test_ds = TabularDataset(rna_df.to_numpy(dtype=np.float32))
                                test_dl = DataLoader(test_ds, batch_size=64, shuffle=False)
                                
                                
                                sc.pp.neighbors(self.adata, n_neighbors=10, n_pcs=40)
                                sc.tl.leiden(self.adata)

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
                            
                                colors_var = reversed(self.adata.obs.columns)
                                
                                citeseq_container = st.empty()

                                color = st.selectbox(label="Colour", options=(colors_var), key="sb_auto_colors")
                                citeseq_container.scatter_chart(plot_df, x="UMAP1", y="UMAP2", color='leiden', size=18)
                                
                                subcol1, _, _, _ = st.columns(4)
                                submit_btn = subcol1.form_submit_button(label="Update colours")
                                
                                if submit_btn:
                                    with st.spinner(text="Replotting chart"):
                                        citeseq_container.scatter_chart(plot_df, x="UMAP1", y="UMAP2", color=color, size=18)

                    elif(isinstance(trained_model, solo_model)):
                        with st.form(key="solo_doublet_form"):
                            st.subheader("Doublet prediction")
                            
                            with st.spinner("Loading Solo doublet predictions"):
                                
                                df_solo = pd.DataFrame({'umap1': self.adata.obsm['X_umap'][:,0], 'umap2': self.adata.obsm['X_umap'][:,1], 'color': self.adata.obs['prediction']})
                            
                                st.scatter_chart(data=df_solo, x='umap1', y='umap2', color='color', size=18)
                                
                                subcol1, _, _, _ = st.columns(4)
                                submit_btn = subcol1.form_submit_button(label="Filter doublets", use_container_width=True)
                                
                                if submit_btn:
                                    self.adata = self.adata[self.adata.obs.prediction == 'singlet']
                                    self.save_adata(name="adata_solo")
                                    st.toast("Removed doublet predictions", icon='✅')


                    else:
                        st.error("Unknown model")

            except Exception as e:
                st.error(e)


    def pca_graph(self):
        """
        Compute PCA coordinates with optional colour.

        Parameters
        ----------
        gene: List[str]
            Genes to supply as color argument in PCA.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/pca_cluster.png

        Example
        -------
        import scanpy as sc

        genes = ['brca1', 'brca2']
        sc.tl.pca(adata, svd_solver='arpack')
        sc.pl.pca(adata, color=genes)
        """
        with self.col1:
            try:
                with st.form(key="pca_cluster_form"):
                    st.subheader("PCA")
                    genes = st.multiselect(label="Gene", options=self.adata.var_names, default=(self.adata.var_names[0]), key="ms_pca_gene", max_selections=24)
                    pca_empty = st.empty()
                    info_container = st.empty()
                    
                    with st.spinner(text="Computing PCA"):
                        sc.tl.pca(self.adata, svd_solver='arpack')
                        colors = st.session_state.ms_pca_gene or self.adata.var_names[0]
                        use_raw = True
                        if self.adata.raw:
                            for var in self.adata.var_names:
                                if not (self.adata.raw.var_names.__contains__(var)):
                                    use_raw = False
                                    info_container.info("Gene ID format has been changed, setting 'use_raw' to false.")
                                    break
                        pca_empty.empty()
                        for color in colors:
                            df = pd.DataFrame({'pca1': self.adata.obsm['X_pca'][:,0], 'pca2': self.adata.obsm['X_pca'][:,1], f'{color}': self.adata.to_df()[color].values})
                            pca_empty.markdown(f"""<p style='text-align: center; size: 16px; font-weight: bold;'>{color}</p>""", unsafe_allow_html=True)
                            pca_empty.scatter_chart(df, x='pca1', y='pca2', size=10, color=f'{color}')
                            

                        
                    subcol1, _, _, _ = st.columns(4)
                    submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
                    if submit_btn:
                        with st.spinner(text="Computing PCA"):
                            sc.tl.pca(self.adata, svd_solver='arpack')
                            colors = genes or [self.adata.var_names[0]]
                            use_raw = True
                            if self.adata.raw:
                                for var in self.adata.var_names:
                                    if not (self.adata.raw.var_names.__contains__(var)):
                                        use_raw = False
                                        info_container.info("Gene ID format has been changed, setting 'use_raw' to false.")
                                        break
                            df = pd.DataFrame({'pca1': self.adata.obsm['X_pca'][:,0], 'pca2': self.adata.obsm['X_pca'][:,1], f'{color}': self.adata.to_df()[color].values})
                            pca_empty.empty()
                            for color in colors:
                                pca_empty.markdown(f"""<p style='text-align: center; size: 16px; font-weight: bold;'>{color}</p>""", unsafe_allow_html=True)
                                pca_empty.scatter_chart(df, x='pca1', y='pca2', size=10, color=f'{color}')


            except Exception as e:
                st.error(e)

    def variance_ratio(self):
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
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/variance_ratio.png

        Example
        -------
        import scanpy as sc
        
        sc.tl.pca(adata, svd_solver='arpack')
        sc.pl.pca_variance_ratio(adata, log=True)
        """
        with self.col1:
            try:
                def do_var_ratio():
                    with st.spinner(text="Computing variance ratio"):
                            
                            empty = st.empty()
                            sc.pl.pca_variance_ratio(adata, log=log)

                            pcs = []
                            pc_val = []

                            for i, pc in enumerate(adata.uns['pca']['variance_ratio']):
                                pcs.append(i+1)
                                pc_val.append(pc)

                            pcs = np.array(pcs)
                            pc_val = np.array(pc_val)
                            df = pd.DataFrame({'PC': pcs, 'value': pc_val})
                            df["PC"] = df["PC"].astype("category")

                            empty.empty()
                            empty.scatter_chart(df[:n_pcs], x='PC', y='value', color="PC")

                with st.form(key="variance_ratio_form"):
                    st.subheader("PCA variance ratio")
                    subcol1, subcol2 = st.columns(2)
                    n_pcs = subcol1.number_input(label="n_pcs", min_value=1, max_value=50, value=30)
                    log = subcol1.checkbox(label="Log", value=True)
                    subcol1, _, _, _ = st.columns(4)
                    submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
                    do_var_ratio()
                    if submit_btn:
                        do_var_ratio()
            
            except Exception as e:
                print("Error ", e)
                st.error(e)



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
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/tsne.png

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
                    perplexity = st.slider(label="Perplexity", min_value=1, max_value=100, value=30)
                    colors_var = reversed(self.adata.obs.columns)
                    
                    tsne_color = st.selectbox(label="Colour", options=(colors_var), key="sb_tsne_colors")
                    tsne_container = st.empty()
                    
                    with st.spinner(text="Computing tSNE"):
                        sc.pp.neighbors(self.adata, n_neighbors=10, n_pcs=40)
                        sc.tl.leiden(self.adata) 
                        sc.tl.tsne(self.adata, perplexity=30)   
                        df_tsne = pd.DataFrame({'tsne1': self.adata.obsm['X_tsne'][:,0], 'tsne2': self.adata.obsm['X_tsne'][:,1], 'color': self.adata.obs[f'leiden']})  
                        tsne_container.scatter_chart(data=df_tsne, x='tsne1', y='tsne2', color='color', size=18)  
                    
                    subcol1, _, _, _ = st.columns(4)
                    submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)

                    if submit_btn:
                        sc.tl.tsne(self.adata, perplexity=perplexity)   
                        df_tsne = pd.DataFrame({'tsne1': self.adata.obsm['X_tsne'][:,0], 'tsne2': self.adata.obsm['X_tsne'][:,1], 'color': self.adata.obs[f'{tsne_color}']})  
                        tsne_container.scatter_chart(data=df_tsne, x='tsne1', y='tsne2', color='color', size=18)  
                        
            except Exception as e:
                st.error(e)



    def neighbourhood_graph(self):
        """
        Compute neighbourhood graph with color mapping for selected genes.

        Parameters
        ----------
        genes: List[str]
            Genes to supply as color argument in PCA.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/neighbourhood_graph.png

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
        plt.style.use('dark_background')
        with self.col2:
            try:
                with st.form(key="nhood_graph_form"):
                    st.subheader("Neighbourhood graph")
                    
                    umap_options = np.append(self.adata.var_names, 'leiden')
                    genes = st.multiselect(label='Gene', options=(umap_options), key="ms_umap_select_gene", default=["leiden", self.adata.var_names[0]], max_selections=24) 
                        
                    nhood_container = st.container()
                    info_container = st.empty() 
                    
                    with st.spinner(text="Computing neighbourhood graph"):
                                
                        colors = genes or 'leiden'
                        use_raw = True
                        if self.adata.raw:
                            for var in self.adata.var_names:
                                if not (self.adata.raw.var_names.__contains__(var)):
                                    use_raw = False
                                    info_container.info("Gene ID format has been changed, setting 'use_raw' to false.")
                                    break

                        sc.pp.neighbors(self.adata, n_neighbors=10, n_pcs=40)
                        sc.tl.leiden(self.adata) 
                        sc.tl.paga(self.adata)
                        sc.pl.paga(self.adata, use_raw=use_raw, plot=False)
                        sc.tl.umap(self.adata, init_pos='paga')

                        for color in colors:
                            if color in self.adata.obs_keys():
                                df = pd.DataFrame({'umap1': self.adata.obsm['X_umap'][:,0], 'umap2': self.adata.obsm['X_umap'][:,1], f'{color}': self.adata.obs[color].values})
                            else:
                                df = pd.DataFrame({'umap1': self.adata.obsm['X_umap'][:,0], 'umap2': self.adata.obsm['X_umap'][:,1], f'{color}': self.adata.to_df()[color].values})
                            nhood_container.markdown(f"""<p style='text-align: center; size: 16px; font-weight: bold;'>{color}</p>""", unsafe_allow_html=True)
                            nhood_container.scatter_chart(df, x='umap1', y='umap2', size=10, color=f'{color}')
                        
                    subcol1, _, _, _ = st.columns(4)
                        
                    submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
                        
                    if submit_btn:
                        with st.spinner(text="Computing neighbourhood graph"):

                            colors = genes or 'leiden'
                            use_raw = True
                            if self.adata.raw:
                                for var in self.adata.var_names:
                                    if not (self.adata.raw.var_names.__contains__(var)):
                                        use_raw = False
                                        info_container.info("Gene ID format has been changed, setting 'use_raw' to false.")
                                        break
                                
                            sc.pp.neighbors(self.adata, n_neighbors=10, n_pcs=40)
                            sc.tl.leiden(self.adata) 
                            sc.tl.paga(self.adata)
                            sc.pl.paga(self.adata, use_raw=use_raw, plot=False)
                            sc.tl.umap(self.adata, init_pos='paga')

                            for color in colors:
                                if color in self.adata.obs_keys():
                                    df = pd.DataFrame({'umap1': self.adata.obsm['X_umap'][:,0], 'umap2': self.adata.obsm['X_umap'][:,1], f'{color}': self.adata.obs[color].values})
                                else:
                                    df = pd.DataFrame({'umap1': self.adata.obsm['X_umap'][:,0], 'umap2': self.adata.obsm['X_umap'][:,1], f'{color}': self.adata.to_df()[color].values})
                                nhood_container.markdown(f"""<p style='text-align: center; size: 16px; font-weight: bold;'>{color}</p>""", unsafe_allow_html=True)
                                nhood_container.scatter_chart(df, x='umap', y='umap2', size=10, color=f'{color}')

            except Exception as e:
                st.error(e)



try:
    adata = st.session_state.adata_state.current.adata

    sidebar = Sidebar()
    sidebar.show()

    analysis = Cluster_analysis(adata)

    analysis.autoencoder_cluster_plot()
    analysis.tsne_graph()
    analysis.neighbourhood_graph()
    analysis.pca_graph()
    analysis.variance_ratio()

    sidebar.show_preview()
    sidebar.export_script()
    sidebar.delete_experiment_btn()
    sidebar.show_version()


except KeyError as ke:
    print("KeyError: ", ke)
    st.error("Couldn't find adata object in session, have you uploaded one?")



