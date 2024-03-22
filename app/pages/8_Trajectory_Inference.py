import pickle
import scanpy as sc
import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time

from models.AdataModel import AdataModel
from components.sidebar import *
import os
from utils.paga import *
from state.StateManager import StateManager

st.set_page_config(layout="wide", page_title='Nuwa', page_icon='🧬')

os.chdir('/app')

with open('css/common.css') as f:
    common_style = f"""
                <style>
                {f.read()}
                </style>
                """
    st.markdown(common_style, unsafe_allow_html=True)

class Trajectory_inference:
    """
    Trajectory inference gives key insights into cell differentiation through ordering the stages of development into a continuous sequence of clusters.

    Notes
    -----
    .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/trajectory_inference_page.png
    """
    def __init__(self, adata: AnnData):
        try:
            self.adata = adata
            self.adata.raw = adata
            #convert to float64
            self.adata.X = self.adata.X.astype('float64')
            #compute pca
            sc.tl.pca(self.adata, svd_solver='arpack')
            #compute neighbours
            sc.pp.neighbors(self.adata, n_neighbors=4, n_pcs=20)
            sc.tl.draw_graph(self.adata)
            sc.tl.louvain(adata)

        except Exception as e:
            st.error(e)
        

    def draw_graph(self):
        """
        Plot graph with force-directed graph drawing [Islam11], [Jacomy14], [Chippada18].

        Parameters
        ----------
        color: str
            Obs value to color samples.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/traj_infer_graph.png
        
        Example
        -------
        import scanpy as sc

        #compute pca
        sc.tl.pca(adata, svd_solver='arpack')
        #compute neighbours
        sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)
        sc.tl.draw_graph(adata)
        #compute louvain clusters
        sc.tl.louvain(adata, resolution=1.0)
        sc.pl.draw_graph(adata, color='louvain', legend_loc='on data')
        """
        with self.col1:
            try:
                st.subheader("Force-directed graph")
                with st.spinner(text="Plotting clusters"):
                    subcol1, _, _, _ = st.columns(4, gap='large')
                    key_index = self.adata.obs_keys().index('louvain') #created in constructor
                    if key_index == -1: #in case louvain was not created
                        key_index = 0
                    subcol1.selectbox(label='color', options=self.adata.obs.columns, key='sb_plt_colors', index=key_index)
                    st.slider(label="Point size", min_value=1, max_value=100, key="sl_traj_graph_size", value=50)
                    #denoise
                    #sc.tl.diffmap(self.adata)
                    #sc.pp.neighbors(self.adata, n_neighbors=10, use_rep='X_diffmap')
                    #sc.tl.draw_graph(self.adata)
                    self.ax_df = pd.DataFrame({'fa1': self.adata.obsm['X_draw_graph_fa'][:,0], 'fa2': self.adata.obsm['X_draw_graph_fa'][:,1], 'color': self.adata.obs[st.session_state.sb_plt_colors]})
                    st.scatter_chart(self.ax_df, x='fa1', y='fa2', color='color', height=600, size=st.session_state.sl_traj_graph_size)

            except Exception as e:
                st.error(e)



    def paga_clustering(self):
        """
        Computes a PAGA graph using louvain/leiden embeddings. View expression levels for individual genes across nodes of the graph. You may also recompute a scatter chart based on PAGA initialization. 

        Parameters
        ----------
        algorithm: str
            Clustering algorithm to using as PAGA groups.

        resolution: float
            Resolution to use in clustering algorithm. A higher resolution produces more clusters.

        threshold: float
            Do not draw edges for weights below this threshold. Set to 0 if you want all edges. 
            Discarding low-connectivity edges helps in getting a much clearer picture of the graph.

        root: str | int
            This is the index of the root node or a list of root node indices. If this is a non-empty vector then the supplied node IDs are used as the roots of the trees (or a single tree if the graph is connected). 
            If this is None or an empty list, the root vertices are automatically calculated based on topological sorting.

        n_neighbours: int
            The size of local neighborhood (in terms of number of neighboring data points) used for manifold approximation. 
            Larger values result in more global views of the manifold, while smaller values result in more local data being preserved. In general values should be in the range 2 to 100.

        n_pcs: int
            Number of PCs.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/paga1.png
        .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/paga2.png
        .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/paga3.png

        Example
        -------
        import scanpy as sc

        adata = sc.datasets.paul15()

        # louvain cluster
        sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)
        sc.tl.louvain(adata, resolution=resolution)

        # run PAGA             
        sc.tl.paga(adata, groups=algorithm)
        sc.pl.paga(adata, color=['louvain', '0610007L01Rik'])

        # recompute paga
        sc.tl.draw_graph(adata, init_pos='paga')
        sc.pl.draw_graph(adata, color=['louvain', '0610007L01Rik'], legend_loc='on data')
        """
        with self.col2:
            try:
                with st.form(key="louvain_cluster_form"):
                    st.subheader("PAGA")

                    #algorithm
                    col1, col2 = st.columns(2)
                    algorithm = col1.radio(label="Clustering algorithm", options=['leiden', 'louvain'])
                    resolution = col2.number_input(label="Resolution", min_value=0.01, value=1.00, step=0.10)
                    #paga
                    col1, col2 = st.columns(2)
                    threshold = col1.number_input(label="Paga threshold", min_value=0.01, value=0.03, step=0.01)
                    root = col2.text_input(label="Root")
                    #neighbours
                    col1, col2 = st.columns(2)
                    n_neighbours = col1.number_input(label="n_neighbours", value=4, step=1, format="%i")
                    n_pcs = col2.number_input(label="n_pcs", value=20, step=1, format="%i")
                    #genes
                    options = st.multiselect(label="Genes", options=self.adata.var_names, default=[self.adata.var_names[0]])
                    subcol1, _, _, _ = st.columns(4)
                    submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
                    
                    if submit_btn:
                        with st.spinner(text="Computing paga clusters"):
                            #louvain cluster
                            sc.pp.neighbors(self.adata, n_neighbors=n_neighbours, n_pcs=n_pcs)

                            if algorithm == "louvain":
                                sc.tl.louvain(self.adata, resolution=resolution)
                                options = np.append(['louvain'], options)
                            elif algorithm == "leiden":
                                sc.tl.leiden(self.adata, resolution=resolution)
                                options = np.append(['leiden'], options)
                            
                            sc.tl.paga(self.adata, groups=algorithm)

                            #paga
                            st.session_state.trajectory_plots["paga"] = compute_paga(self.adata, color=options, height=500)
                            if not root:
                                root = 0
                            sc.pl.paga(adata, threshold=threshold, root=root)

                            #recompute paga
                            sc.tl.draw_graph(self.adata, init_pos='paga')
                            sc.pl.draw_graph(self.adata, color=options, legend_loc='on data')
                            st.session_state.trajectory_plots["paga_cluster_colors"] = options
                            st.session_state.trajectory_plots["paga_cluster_embedding"] = self.adata.obsm["X_draw_graph_fa"]

                            st.toast("Added PAGA embeddings", icon="✅")

                            
                            #write to script state
                            # st.session_state["script_state"].add_script("sc.tl.louvain(adata, resolution=1.0)")
                            # st.session_state["script_state"].add_script("sc.tl.paga(adata, groups='louvain')")
                            # st.session_state["script_state"].add_script(f"sc.pl.paga(adata, colors={options}, cmap='viridis', node_size_scale=3)")
                            # st.session_state["script_state"].add_script("plt.show()")
                        

            except Exception as e:
                st.error(e)


    def draw_plots(self):
        st.subheader("Plots")
        paga, paga_embed, dpt = st.tabs(['Paga', 'Paga Embedding', 'Diffusion Pseudotime'])
        with paga:
            if st.session_state.trajectory_plots["paga"] == None:
                st.info("You must run paga first.")
            else:
                figs = st.session_state.trajectory_plots["paga"]
                col1, col2 = st.columns(2)
                columns = [col1, col2]
                for i, fig in enumerate(figs):
                    with columns[i % 2]:
                        st.plotly_chart(fig, use_container_width=True)
        
        with paga_embed:
            if not isinstance(st.session_state.trajectory_plots["paga_cluster_colors"], np.ndarray):
                st.info("You must run paga first.") 
            else:
                colors = st.session_state.trajectory_plots["paga_cluster_colors"]
                col1, col2 = st.columns(2)
                columns = [col1, col2]
                for i, color in enumerate(colors):
                    if color in self.adata.obs:
                        df = pd.DataFrame({'fa1': st.session_state.trajectory_plots["paga_cluster_embedding"][:, 0], 'fa2': st.session_state.trajectory_plots["paga_cluster_embedding"][:, 1], f'{color}': self.adata.obs['louvain']})
                        with columns[i % 2]:
                            st.markdown(f"""<p style='font-size: 16px; font-weight: bold;'>{color}</p>""", unsafe_allow_html=True)
                            st.scatter_chart(df, x='fa1', y='fa2', color=f'{color}', size=10, height=500, use_container_width=True)
                    else:
                        df = pd.DataFrame({'fa1': st.session_state.trajectory_plots["paga_cluster_embedding"][:,0], 'fa2': st.session_state.trajectory_plots["paga_cluster_embedding"][:,1], f'{color}': self.adata.to_df()[color].values})
                        with columns[i % 2]:
                            st.markdown(f"""<p style='font-size: 16px; font-weight: bold;'>{color}</p>""", unsafe_allow_html=True)
                            st.scatter_chart(df, x='fa1', y='fa2', color=f'{color}', size=10, height=500, use_container_width=True)

        with dpt:
            if st.session_state.trajectory_plots["dpt_cluster_color"] == None:
                st.info("You must run dpt first.") 
            else:
                color = st.session_state.trajectory_plots["dpt_cluster_color"]
                col1, col2 = st.columns(2)
                
                with col1:
                    df = pd.DataFrame({'fa1': st.session_state.trajectory_plots["paga_cluster_embedding"][:, 0], 'fa2': st.session_state.trajectory_plots["paga_cluster_embedding"][:, 1], f'{color}': self.adata.obs[color]})
                    st.markdown(f"""<p style='font-size: 16px; font-weight: bold;'>{color}</p>""", unsafe_allow_html=True)
                    st.scatter_chart(df, x='fa1', y='fa2', color=f'{color}', size=10, height=500, use_container_width=True)

                with col2:
                    df = pd.DataFrame({'fa1': st.session_state.trajectory_plots["paga_cluster_embedding"][:, 0], 'fa2': st.session_state.trajectory_plots["paga_cluster_embedding"][:, 1], 'dpt_pseudotime': self.adata.obs['dpt_pseudotime']})
                    st.markdown(f"""<p style='font-size: 16px; font-weight: bold;'>Diffusion pseudotime</p>""", unsafe_allow_html=True)
                    st.scatter_chart(df, x='fa1', y='fa2', color='dpt_pseudotime', size=10, height=500, use_container_width=True)
                


    def draw_page(self):
        st.title("Trajectory Inference")
        #graphs
        if 'trajectory_plots' not in st.session_state:
            st.session_state["trajectory_plots"] = dict(paga=None, paga_cluster_colors=None, dpt_cluster_color=None, paga_cluster_embedding=None)
        #columns
        self.col1, self.col2 = st.columns(2, gap="medium")
        self.draw_graph()
        self.paga_clustering()
        self.diffusion_pseudotime()
        #self.show_path()
        self.draw_plots()


    def diffusion_pseudotime(self):
        """
        Infer progression of cells through geodesic distance along the graph [Haghverdi16]_ [Wolf19]_.

        Parameters
        ----------
        algorithm: str
            Clustering algorithm to use as embedding.

        root: int | str
            This is the index of the root node or a list of root node indices. use this node as the root cell for diffusion pseudotime. 

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/dpt1.png
        .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/dpt2.png

        Example
        -------
        import scanpy as sc

        sc.tl.louvain(adata)
        adata.uns['iroot'] = np.flatnonzero(adata.obs['louvain'] == '16')[0]
        sc.tl.dpt(adata)
        sc.pl.draw_graph(adata, color=['louvain', 'dpt_pseudotime'], legend_loc='on data', cmap='viridis')
        """
        with self.col2:
            try:
                with st.form(key="dpt_form"):
                    st.subheader("Diffusion Pseudotime")
                    algorithm = st.radio(label="Clustering algorithm", options=['leiden', 'louvain'])
                    root = st.selectbox(label='Root cell', options=(self.adata.obs['louvain'].unique()))
                    subcol1, _, _, _ = st.columns(4)
                    dpt_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
                    if dpt_btn:
                        with st.spinner(text="Computing dpt"):
                            if algorithm not in self.adata.obs:
                                if algorithm == "leiden":
                                    sc.tl.leiden(self.adata)
                                if algorithm == "louvain":
                                    sc.tl.louvain(self.adata)
                            
                            st.session_state["trajectory_plots"]["dpt_cluster_color"] = algorithm

                            self.adata.uns['iroot'] = np.flatnonzero(self.adata.obs[algorithm]  == root)[0]
                            
                            sc.tl.dpt(self.adata)
                            #sc.pp.log1p(self.adata)
                            #sc.pp.scale(self.adata)

                            sc.pl.draw_graph(self.adata, color=['louvain', 'dpt_pseudotime'], legend_loc='on data', cmap='viridis')

                            st.toast("Finished dpt analysis", icon="✅")
                            
                            #add to script state
                            
                            # st.session_state["script_state"].add_script("adata.uns['iroot'] = np.flatnonzero(adata.obs['louvain']  == st.session_state.sb_root_cell)[0]")
                            # st.session_state["script_state"].add_script("sc.tl.dpt(adata)")
                            # st.session_state["script_state"].add_script("sc.pp.log1p(adata)")
                            # st.session_state["script_state"].add_script("sc.pp.scale(adata)")
                            # st.session_state["script_state"].add_script("sc.pl.draw_graph(adata, color=['louvain', 'dpt_pseudotime'], legend_loc='on data', cmap='viridis')")
                            # st.session_state["script_state"].add_script("plt.show()")

                            #     sc.tl.dpt(adata_raw, n_branchings=1, n_dcs=10)
                            #     adata_raw.obs.dpt_groups = adata_raw.obs.dpt_groups[:60]
                            #     ax2 = sc.pl.dpt_timeseries(adata_raw)
                            #     st.pyplot(ax2)

            except Exception as e:
                st.error(e)

    def show_path(self):
        with self.col2:
            try:
                st.subheader("Paths")


                paths = [('erythrocytes', [16, 12, 7, 13, 18, 6, 5, 10]),
                        ('neutrophils', [16, 0, 4, 2, 14, 19]),
                        ('monocytes', [16, 0, 4, 11, 1, 9, 24])]


                self.adata.obs['distance'] = self.adata.obs['dpt_pseudotime']
                self.adata.obs['clusters'] = self.adata.obs['louvain']  # just a cosmetic change
                self.adata.uns['clusters_colors'] = self.adata.uns['louvain_colors']

                gene_names = ['Gata2', 'Gata1', 'Klf1', 'Epor', 'Hba-a2',  # erythroid
                'Elane', 'Cebpe', 'Gfi1',                    # neutrophil
                'Irf8', 'Csf1r', 'Ctsg'] 

        
                _, axs = plt.subplots(ncols=3, figsize=(6, 2.5), gridspec_kw={'wspace': 0.05, 'left': 0.12})
                plt.subplots_adjust(left=0.05, right=0.98, top=0.82, bottom=0.2)


                for ipath, (descr, path) in enumerate(paths):
                    ax = sc.pl.paga_path(
                        self.adata, path, gene_names,
                        show_node_names=False,
                        ax=axs[ipath],
                        ytick_fontsize=12,
                        left_margin=0.15,
                        n_avg=50,
                        annotations=['distance'],
                        show_yticks=True if ipath==0 else False,
                        show_colorbar=False,
                        color_map='Greys',
                        groups_key='clusters',
                        color_maps_annotations={'distance': 'viridis'},
                        title='{} path'.format(descr),
                        return_data=True,
                        show=False)
                    #data.to_csv('./write/paga_path_{}.csv'.format(descr))
                
                st.pyplot(axs)
                plt.savefig('./figures/paga_path_paul15.pdf')

            except Exception as e:
                st.error(e)
                

try:
    adata = st.session_state.adata_state.current.adata
    sidebar = Sidebar()

    sidebar.show()

    tji = Trajectory_inference(adata)

    tji.draw_page()



except Exception as e:
    if(st.session_state == {}):
        StateManager().load_session()
        st.rerun()
    else:
        st.toast(e, icon="❌")





