import streamlit as st
import scanpy as sc
from anndata import AnnData
from components.sidebar import *
import os
import numpy as np
import pandas as pd
import scanpy as sc
import gseapy
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import altair as alt
import plotly.graph_objects as go

st.set_page_config(layout="wide", page_title='Nuwa', page_icon='🧬')

with open('css/common.css') as f:
    common_style = f"""
                <style>
                {f.read()}
                </style>
                """
    st.markdown(common_style, unsafe_allow_html=True)


class Differential_gene_expression:
    """
    Differential gene expression looks at how genes are expressed compared to the rest of the dataset. This includes useful matrix plots, dot plots and violin vlots to visualise variable expression. You can also choose which clusters and statistical tests to run the DE analysis.

    Notes
    -----
    .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/de_page.png
    """
    def __init__(self, adata):
        self.adata = adata
        self.draw_page()

    def save_adata(self):
        sc.write(filename=os.path.join(os.getenv('WORKDIR'), 'adata', st.session_state.adata_state.current.adata_name), adata=self.adata)
        st.session_state.adata_state.current.adata = self.adata

    def draw_page(self):
        col1, col2, col3 = st.columns(3)
        with col1:
            self.upload_marker_genes()
            self.preview_marker_genes()
        with col2:
            self.add_embeddings()
            self.stat_tests()
            

        self.visualize()
        self.rank_genes_groups()
        self.show_top_ranked_genes()
        with col3:
            self.do_umap()
            
        

    def stat_tests(self):
        """
        Compute ranking of differential genes between groups using statistical tests such as t-test.

        Parameters
        ----------
        method: str
            Statistical test to find differential clusters and results in adata for plotting. Options include t-test, t-test_overestim_var, wilcoxon and logreg.

        group_by: str
            Obs value to group clusters by.
            
        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/stats_test.png
        
        Example
        -------
        import scanpy as sc

        sc.tl.rank_genes_groups(adata, groupby=group_by, method = "wilcoxon")
        sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save=True)
        """
        try:
            with st.form(key="stat_tests_form"):
                st.subheader("Run statistical tests")
                method = st.selectbox(label="Method", options=(['logreg', 't-test', 't-test_overestim_var', 'wilcoxon']), key='rb_method')
                group_by = st.selectbox(label="Group by", options=st.session_state.adata_state.current.adata.obs_keys())
                marker_genes_container = st.empty()
                subcol1, _, _ = st.columns(3)
                submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
                if submit_btn:
                    marker_genes_container.empty()
                    with st.spinner(text="Computing tests"):
                        sc.tl.rank_genes_groups(self.adata, groupby=group_by, method = str.lower(method), key_added=method)
                        #sc.pl.rank_genes_groups(st.session_state.adata_state.current.adata, n_genes=25, sharey=False, save=True)
                        self.save_adata()

        except Exception as e:
            st.error(e)

    def upload_marker_genes(self):
        """
        Upload a csv file of marker genes to be used for cell type identification.
        """
        try:
            st.subheader("Upload marker genes file", help="")
            marker_genes_upload = st.file_uploader(label="Upload csv file", type=['csv', 'tsv'])
            if marker_genes_upload:
                st.session_state['marker_genes_df'] = pd.read_csv(marker_genes_upload)
            else:
                st.code("Example format:\nmarkers,cell_type\nIL7R CCR7, Naive CD4+ T", language="csv")

        except Exception as e:
            print("Error: ", e)
            st.error(e)

    def preview_marker_genes(self):
        try:
            st.subheader("Preview")
            if 'marker_genes_df' not in st.session_state:
                st.info("Upload a file")
                #st.markdown("""<div style='display: flex; flex-direction: column;'><p style='text-align: center; margin-top: 100px; color: #b0b0b0; whitespace: pre-line; justify-content: center; align-content: center; align-items: center;'>Upload a file</p></div>""", unsafe_allow_html=True)
            else:
                st.dataframe(st.session_state.marker_genes_df, width=500, height=220)
        except Exception as e:
            print("Error: ", e)

    def do_umap(self):
        """
        Perform a UMAP to visualise clusters.

        Parameters
        ----------
        color: str
            Obs key to use as a color.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/umap_de.png

        Example
        -------
        import scanpy as sc

        sc.pp.neighbors(adata, random_state=42)
        sc.tl.umap(adata, random_state=42)
        sc.pl.umap(adata, color="leiden")
        """
        try:
            with st.form(key="do_umap"):
                st.subheader("UMAP")
                def run_umap(adata):
                    with st.spinner(text="Running UMAP"):
                        sc.pp.neighbors(adata, random_state=42)
                        sc.tl.umap(adata, random_state=42)
                        df_umap = pd.DataFrame({'umap1': adata.obsm['X_umap'][:,0], 'umap2': adata.obsm['X_umap'][:,1], 'color': adata.obs[f'{st.session_state.sb_umap_color_dge}']})  
                        umap_empty.empty()
                        umap_empty.scatter_chart(data=df_umap, x='umap1', y='umap2', color='color', size=18)
                    
                
                index = 0      
                for i, item in enumerate(self.adata.obs_keys()):
                    if item.lower().replace("_", "").__contains__("leiden") or item.lower().replace("_", "").__contains__("louvain") or item.lower().replace("_", "").__contains__("celltype"): #give precedence to batch if present since it is relevant to preprocessing
                        index = i           
                umap_color = st.selectbox(label="Color", options=self.adata.obs_keys(), key="sb_umap_color_dge", index=index)
                subcol1, _, _ = st.columns(3)
                umap_dge_btn = subcol1.form_submit_button("Apply", use_container_width=True)
                umap_empty = st.empty()
                
                run_umap(self.adata)
                
                if umap_dge_btn:
                    run_umap(self.adata)

        except Exception as e:
            print("error: ", e)
            st.error(e)

    def add_embeddings(self):
        """
        Add leiden or louvain embeddings to dataset if not already present or for a specific resolution.

        Parameters
        ----------
        algorithm: str, default = "Leiden"
            Clustering algorithm to use.
        
        resolution: float, default = 0.6
            Resolution to use in clustering algorithm (a higher resultion finds more clusters).

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/cluster_embeddings.png

        Example
        -------
        import scanpy as sc

        resolution = 0.6
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
        # add leiden embeddings
        sc.tl.leiden(adata, resolution=resolution, key_added=f"leiden_{resolution}") 
        # add louvain embeddings
        sc.tl.louvain(adata, resolution=resolution, key_added=f"louvain_{resolution}")
        """
        try:
            with st.form(key="add_embeddings_form"):
                st.subheader("Cluster embeddings", help="Embed clusters in dataset for DGE.")
                algorithm = st.selectbox(label="Algorithm", options=['Leiden', 'Louvain'])
                resolution = st.number_input(label="Resolution", min_value=0.1, value=0.6, format="%.1f")
                subcol1, _, _ = st.columns(3)
                submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
                if submit_btn:
                    with st.spinner(text=f"Computing {algorithm} clusters"):
                        sc.pp.neighbors(self.adata, n_neighbors=10, n_pcs=40)
                        if algorithm == "Leiden":
                            sc.tl.leiden(self.adata, resolution=resolution, key_added=f"leiden_{resolution}") 
                        if algorithm == "Louvain":
                            sc.tl.louvain(self.adata, resolution=resolution, key_added=f"louvain_{resolution}")
                    
                    self.save_adata()

        except Exception as e:
            st.error(e)



    def visualize(self):
        """
        Visualize differential expression between individual genes or clusters.

        

        """
        try:
 
            with st.form(key="visualize_form"):

                st.subheader("Visualize plots")

                input_col1, input_col2, input_col3, _, _ = st.columns(5, gap="large")
                method = input_col1.radio(label="Method", options=['t-test', 't-test_overestim_var', 'wilcoxon', 'logreg'])
                n_genes = input_col2.text_input(label="n_genes", value=5, help="Number of genes to display in each cluster.")
                group_by = input_col2.selectbox(label="Group by", options=st.session_state.adata_state.current.adata.obs_keys())
                subcol1, _, _, _, _, _, _, _, _ = st.columns(9)
                viz_submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
                
                if viz_submit_btn:
                    plt.style.use('dark_background')
                    for obs in st.session_state.adata_state.current.adata.obs_keys():
                        if obs.__contains__("leiden") or obs.__contains__("louvain"):
                            st.session_state.adata_state.current.adata.obs[obs] = st.session_state.adata_state.current.adata.obs[obs].astype('category')
                    n_genes = int(n_genes)
                    heatmap, dotplot, stacked_violins, matrix_plot = st.tabs(['Heatmap', 'Dotplot', 'Stacked violins', 'Matrix plot'])
                    with st.spinner(text="Generating plots"):
                        with heatmap:
                            heatmap_ax = sc.pl.rank_genes_groups_heatmap(st.session_state.adata_state.current.adata, n_genes=n_genes, key=method, groupby=group_by, show_gene_labels=True)
                            st.pyplot(heatmap_ax)
                        with dotplot:
                            dotplot_ax = sc.pl.rank_genes_groups_dotplot(st.session_state.adata_state.current.adata, n_genes=n_genes, key=method, groupby=group_by)
                            st.pyplot(dotplot_ax)
                        with stacked_violins:
                            stacked_violins_ax = sc.pl.rank_genes_groups_stacked_violin(st.session_state.adata_state.current.adata, n_genes=n_genes, key=method, groupby=group_by)
                            st.pyplot(stacked_violins_ax)
                        with matrix_plot:
                            matrix_plot_ax = sc.pl.rank_genes_groups_matrixplot(st.session_state.adata_state.current.adata, n_genes=n_genes, key=method, groupby=group_by)
                            st.pyplot(matrix_plot_ax)
   

        except Exception as e:
            st.error(e)


    def rank_genes_groups(self):
        """
        Rank gene groups by a cluster key according to their differential expression accross clusters or genes.

        Parameters
        ----------
        method: str
            Statistical test to find differential clusters and results in adata for plotting. Options include t-test, t-test_overestim_var, wilcoxon and logreg.

        group_by: str
            Cluster name to group by.

        number_of_genes: int
            Number of genes to display.

        reference: str, default = "rest"
            Group to compare clusters to.

        compare_group_1: str
            First comparison group (violin plot).

        compare_group_2: str
            Second comparison group (violin plot). 

        genes: List[str]
            List of genes to compare across clusters.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/elbow_plot.png
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/violin_plot_de.png
        .. image:: https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/violin_plot_de_specific.png

        Example
        -------

        """
        try:
            elbow_plots, violin_plots, violin_plots_specific_genes = st.tabs(['Elbow plot', 'Violin plot', 'Violin plot (specific genes)'])

            with elbow_plots:
                with st.form(key="rank_genes_groups_elbow"):
                    st.subheader("Rank genes groups")
                    subcol1, subcol2, subcol3, _, _ = st.columns(5, gap="large")
                    method = subcol1.radio(label="Method", options=['t-test', 't-test_overestim_var', 'wilcoxon', 'logreg'])
                    group_by = subcol2.selectbox(label="Group by", options=st.session_state.adata_state.current.adata.obs_keys())
                    n_genes = subcol2.number_input(label="Number of genes", min_value=1, value=25, format="%i")
                    reference = subcol3.text_input(label="Reference", value="rest")
                    use_raw = subcol3.toggle(label="Use raw", value=False)
                    subcol1, _, _, _, _, _, _, _, _ = st.columns(9)
                    submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
                    if submit_btn:
                        with st.spinner(text="Calculating plots"):
                            sc.tl.rank_genes_groups(st.session_state.adata_state.current.adata, groupby=group_by, method=method, use_raw=use_raw, reference=reference)
                            #n_graphs = len(st.session_state.adata_state.current.adata.uns['rank_genes_groups']['names'][0])
                            subcol1_graph, subcol2_graph = st.columns(2, gap="large")
                            columns = [subcol1_graph, subcol2_graph]
                            reference = st.session_state.adata_state.current.adata.uns['rank_genes_groups']['params']['reference']

                            names_df_all = pd.DataFrame(st.session_state.adata_state.current.adata.uns['rank_genes_groups']['names'])
                            scores_df_all = pd.DataFrame(st.session_state.adata_state.current.adata.uns['rank_genes_groups']['scores'])
                            pvals_df_all = pd.DataFrame(st.session_state.adata_state.current.adata.uns['rank_genes_groups']['pvals'])
                            pvals_adj_df_all =pd.DataFrame(st.session_state.adata_state.current.adata.uns['rank_genes_groups']['pvals_adj'])
                            logfoldchanges_df_all = pd.DataFrame(st.session_state.adata_state.current.adata.uns['rank_genes_groups']['logfoldchanges'])

                            for i, column in enumerate(names_df_all):
                                with columns[i % 2]:
                                    st.markdown(f"""<div style='margin-left: 20px; display: flex; align-items: center; justify-content: center;'><h1 style='text-align: center; font-size: 2rem;'>{i} vs {reference}</h1></div>""", unsafe_allow_html=True)
                                    df = pd.DataFrame({'gene': names_df_all[column], 'score': scores_df_all[column], 'p value': pvals_df_all[column], 'p value adj': pvals_adj_df_all[column], 'logfoldchanges': logfoldchanges_df_all[column]})
                                    df["gene"] = df["gene"].astype("category")
                                    altair_chart = alt.Chart(df[:n_genes]).mark_circle(size=60).encode(
                                        x=alt.X('gene', type='nominal', sort=None),
                                        y='score',
                                        color=alt.Color('gene', legend=None),
                                        tooltip=['gene', 'score', 'p value', 'p value adj']
                                    ).interactive()
                                    st.altair_chart(altair_chart, use_container_width=True)

            with violin_plots:
                col_group1, _, _, _ = st.columns(4)
                col_group1.selectbox(label="Group by", options=st.session_state.adata_state.current.adata.obs_keys(), key="sb_violin_cluster_group")
                with st.form(key="rank_genes_groups_violin"):
                    st.subheader("Violin plots")
                    subcol1, subcol2, subcol3, _, _ = st.columns(5, gap="large")
                    num_genes = subcol2.number_input(label="Number of genes", min_value=1, value=6, format="%i")
                    cluster1 = subcol1.selectbox(label="Compare group 1", options=np.sort(st.session_state.adata_state.current.adata.obs[st.session_state.sb_violin_cluster_group].unique()), key="sb_cluster1_violin")
                    cluster2 = subcol1.selectbox(label="Compare group 2", options=np.append('rest', np.sort(st.session_state.adata_state.current.adata.obs[st.session_state.sb_violin_cluster_group].unique())), key="sb_cluster2_violin")
                    subcol1, _, _, _, _, _, _, _, _ = st.columns(9)
                    submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)

                    if submit_btn:
                        with st.spinner(text="Calculating plots"):
                            sc.tl.rank_genes_groups(st.session_state.adata_state.current.adata, n_genes=n_genes, groupby=st.session_state.sb_violin_cluster_group, reference=cluster2)
                            sc.pl.rank_genes_groups_violin(st.session_state.adata_state.current.adata, n_genes=n_genes, show=False)
                            subcol1_graph, subcol2_graph = st.columns(2, gap="large")
                            columns = [subcol1_graph, subcol2_graph]

                            st.markdown(f"""<div style='margin-left: 20px; display: flex; align-items: center; justify-content: center;'><h1 style='text-align: center; font-size: 2rem;'>{cluster1} vs {cluster2}</h1></div>""", unsafe_allow_html=True)

                            df = pd.DataFrame()

                            gene_list_top_ranking = st.session_state.adata_state.current.adata.uns['rank_genes_groups']['names'][:num_genes]

                            for gene in gene_list_top_ranking:
                                gene = gene[0]
                                expr = st.session_state.adata_state.current.adata.to_df()[gene]
                                cluster = st.session_state.adata_state.current.adata.obs[st.session_state.sb_violin_cluster_group]
                                sub_df = pd.DataFrame({'gene': np.repeat(gene, len(expr)), 'expression': expr, 'cluster': cluster})
                                df = pd.concat([df, sub_df])
                                del sub_df
                                del cluster
                                del expr
                                del gene


                            fig = go.Figure()

                            fig.add_trace(go.Violin(x=df['gene'][df['cluster'] == f'{cluster1}' ], 
                                y=df['expression'][ df['cluster'] == f'{cluster1}' ],
                                legendgroup=f'{cluster1}', scalegroup=f'{cluster1}', name=f'{cluster1}',
                                side='negative',
                                bandwidth=0.4,
                                line_color='grey')
                            )

                         
                            if cluster2 != "rest":
                                fig.add_trace(go.Violin(x=df['gene'][df['cluster'] == f'{cluster2}'],
                                    y=df['expression'][df['cluster'] == f'{cluster2}'],
                                    legendgroup=f'{cluster2}', scalegroup=f'{cluster2}', name=f'{cluster2}',
                                    side='positive',
                                    bandwidth=0.4,
                                    line_color='#fc0377')
                                )
                            elif cluster2 == "rest":
                                fig.add_trace(go.Violin(x=df['gene'][df['cluster'] != f'{cluster1}'],
                                    y=df['expression'][df['cluster'] != f'{cluster1}'],
                                    legendgroup=f'{cluster2}', scalegroup=f'{cluster2}', name=f'{cluster2}',
                                    side='positive',
                                    bandwidth=0.4,
                                    line_color='#fc0377')
                                )
                            
                            fig.update_traces(meanline_visible=True)
                            fig.update_layout(violingap=0, violinmode='overlay', xaxis_title="Gene", yaxis_title="Expression", legend_title=st.session_state.sb_violin_cluster_group)

                            st.plotly_chart(fig, use_container_width=True)

            with violin_plots_specific_genes:
                with st.form(key="violin_plots_specific_genes"):
                    st.subheader("Measure expression in genes across clusters")
                    subcol1, subcol2, subcol3, _ = st.columns(4, gap="large")
                    cluster = subcol1.selectbox(label="Group", options=st.session_state.adata_state.current.adata.obs_keys())
                    genes = subcol2.multiselect(label="Genes", options=np.sort(st.session_state.adata_state.current.adata.var_names))
                    subcol1, _, _, _, _, _, _, _, _ = st.columns(9)
                    submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
                    line_colors = ['#d1454c', '#8ee065', '#eda621', '#f071bf', '#9071f0', '#71e3f0', '#2f39ed', '#ed2f7b']

                    if submit_btn:

                        for i, gene in enumerate(genes):
                            fig = go.Figure()

                            groups = np.sort(st.session_state.adata_state.current.adata.obs[cluster].unique())
                            for j, group in enumerate(groups):
                                fig.add_trace(go.Violin(x=st.session_state.adata_state.current.adata.obs[cluster][st.session_state.adata_state.current.adata.obs[cluster] == group],
                                    y=st.session_state.adata_state.current.adata.to_df()[gene],
                                    legendgroup=group, scalegroup=group, name=group,
                                    bandwidth=0.4,
                                    line_color=line_colors[j % 8]
                                    )
                                )
                                
                            fig.update_traces(meanline_visible=True)
                            fig.update_layout(violingap=0, violinmode='overlay', xaxis_title=cluster, yaxis_title="Expression", legend_title=cluster)

                            st.markdown(f"""<div style='margin-left: 20px; display: flex; align-items: center; justify-content: center;'><h1 style='text-align: center; font-size: 2rem;'>{gene}</h1></div>""", unsafe_allow_html=True)

                            st.plotly_chart(fig, use_container_width=True)


        except Exception as e:
            print("Error ", e)
            st.error(e)


    def show_top_ranked_genes(self):
        try:
            with st.form(key="top_ranked_genes_form"):
                st.subheader("Top ranked genes")
                subcol1, _, _, _, _ = st.columns(5, gap="large")
                num_of_genes = subcol1.number_input(label="Number of genes", min_value=1, step=1, format="%i", value=10)
                show_p_value = st.toggle(label="Show P value", value=False)
                subcol1, _, _, _, _, _, _, _, _ = st.columns(9)
                submit_btn = subcol1.form_submit_button(label="Run", use_container_width=True)
                empty = st.empty()
                if submit_btn:
                    result = st.session_state.adata_state.current.adata.uns['rank_genes_groups']
                    groups = result['names'].dtype.names
                    df_p_values = pd.DataFrame(
                        {group + '_' + key[:1]: result[key][group]
                        for group in groups for key in ['names', 'pvals']}).head(num_of_genes)


                    df = pd.DataFrame(st.session_state.adata_state.current.adata.uns['rank_genes_groups']['names']).head(num_of_genes)

                    empty.empty()
                    if show_p_value:
                        empty.dataframe(df_p_values)
                    else:
                        empty.dataframe(df)

        except Exception as e:
            print("Error: ", e)
            st.error(e)


try:
    sidebar = Sidebar()
    sidebar.show()

    st.title("Differential Gene Expression")

    sidebar.show_preview()

    dge = Differential_gene_expression(st.session_state.adata_state.current.adata.copy())

    sidebar.export_script()
    sidebar.delete_experiment_btn()
    sidebar.show_version()
    
except KeyError as ke:
    print("KeyError: ", ke)
    st.error("Couldn't find adata object in session, have you uploaded one?")
except Exception as e:
    print("Error: ", e)
    st.error(e)