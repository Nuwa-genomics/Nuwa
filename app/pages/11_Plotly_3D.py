import streamlit as st
import os
import pandas as pd
import plotly.express as px
import scanpy as sc

st.set_page_config(layout="wide", page_title='Nuwa', page_icon='🧬')

os.chdir('/app')

with open('css/common.css') as f:
    common_style = f"""
                <style>
                {f.read()}
                </style>
                """
    st.markdown(common_style, unsafe_allow_html=True)

class Plotly_3D:
    """
    An interactive 3D chart to view cluster embeddings.

    Notes
    -----
    .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/plotly_page.png
    """
    def __init__(self, adata, df):
        self.adata = adata
        self.df = df

        self.plot_chart()

    def plot_chart(self):
        """
        Render plotly chart from given embeddings.

        Parameters
        ----------
        dataset: str
            The dataset for analysis.
        
        cluster_map: str
            The algorithm or obsm name containing coordinates. Currently does not support tSNE for 3D plotting.

        color: str
            The colour for clusters.

        point_size: float
            Adjustable point size for scatter chart.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/plotly1.png
        .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/plotly2.png

        Example
        -------
        import scanpy as sc
        import plotly.express as px

        # compute neighbourhood cluster embeddings with leiden
        sc.pp.neighbors(adata)
        sc.tl.leiden(adata)
        sc.tl.umap(adata, n_components=3)

        # create dataframe
        color = 'cell_type'
        df = pd.DataFrame({'umap1': adata.obsm['X_UMAP'][:, 0], 'umap2': adata.obsm['X_UMAP'][:, 1], 'UMAP3': adata.obsm['X_UMAP'][:, 2], 'color': adata.obs[color]})

        # plot chart
        fig = px.scatter_3d(df, x=df.columns[0], y=df.columns[1], z=df.columns[2], color=df.columns[3], width=2000, height=900)
        fig.update_traces(marker_size=1.0)
        fig.update_layout(title=dict(text=f"{df.columns[0][:-1]} clusters with {df.columns[3]}", font=dict(size=50), automargin=True, yref='paper'))
        fig.update_layout(legend= {'itemsizing': 'constant'})
        fig.show()
        """
        with st.spinner(text="Plotting chart"):
            fig = px.scatter_3d(self.df, x=self.df.columns[0], y=self.df.columns[1], z=self.df.columns[2], color=self.df.columns[3], width=2000, height=900)
            fig.update_traces(marker_size=st.session_state.plotly_point_size)
            fig.update_layout(title=dict(text=f"{self.df.columns[0][:-1]} clusters with {self.df.columns[3]}", font=dict(size=50), automargin=True, yref='paper'))
            fig.update_layout(legend= {'itemsizing': 'constant'})
            st.plotly_chart(fig, use_container_width=True)


try:
    
    adata = st.session_state.adata_state.current.adata

    with st.spinner(text="Computing embeddings"):
        sc.pp.neighbors(st.session_state.adata_state.current.adata)
        sc.tl.leiden(st.session_state.adata_state.current.adata)
        sc.tl.umap(st.session_state.adata_state.current.adata, n_components=3)


    def set_adata():
        if st.session_state.adata_state.switch_adata(st.session_state.sb_adata_selection) == -1:
            st.error("Couldn't switch adata")

    def change_embeddings():
        with st.spinner("Mapping embeddings"):
            algo = st.session_state.plotly_algorithm
            color = st.session_state.plotly_embedding_color
            st.session_state['plotly_df'] = pd.DataFrame({f'{algo}1': st.session_state.adata_state.current.adata.obsm[f'X_{algo}'][:, 0], 
                f'{algo}2': st.session_state.adata_state.current.adata.obsm[f'X_{algo}'][:, 1], 
                f'{algo}3': st.session_state.adata_state.current.adata.obsm[f'X_{algo}'][:, 2], 
                f'{color}': st.session_state.adata_state.current.adata.obs[color]})



    with st.sidebar:
        st.selectbox(label="Dataset:", options=st.session_state.adata_state.get_adata_options(), key="plotly_dataset", index=st.session_state.adata_state.get_index_of_current())
        cluster_options = ['umap', 'pca'] #tsne currently does not support 3rd dimension
        if 'X_citeseq_3d' in adata.obsm:
            cluster_options.append('citeseq_3d')
        st.selectbox(label="Cluster map", options=cluster_options, key="plotly_algorithm", index=0)
        cluster_color_index = 0
        for i, name in enumerate(adata.obs_keys()):
            if name == "leiden":
                cluster_color_index = i
                break
        st.selectbox(label="Embedding colors", options=adata.obs_keys(), key="plotly_embedding_color", index=cluster_color_index)
        st.slider(label="Point size", min_value=0.5, max_value=5.0, step=0.1, value=1.0, key="plotly_point_size")
        st.button(label="Plot chart", use_container_width=True, on_click=change_embeddings)
            

    if 'plotly_df' in st.session_state:
        plotly_3d = Plotly_3D(adata, st.session_state.plotly_df)
    else:
        st.info("Run compute plot with selected embeddings to generate a 3D plot.")

    

except Exception as e:
    print("Error: ", e)
    st.error(e)