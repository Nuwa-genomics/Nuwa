from typing import Optional
import plotly.graph_objects as go
import numpy as np
import scanpy as sc
from anndata import AnnData
from scipy.sparse import issparse
import pandas as pd
import streamlit as st

def highest_expr_genes_box_plot(adata: AnnData, n_top: int, gene_symbols: Optional[str] = None):
    """
    Shows highest expressed n genes in a box plot. This is an alternative to scanpy's highest_expr_genes plotting function with seaborn.

    Parameters
    ----------
    adata: Anndata
        Anndata object to perform analysis on.
    n_top: int
        Number of genes to display.
    gene_symbols: Optional[str], default=None
        Alternative gene symbols to use.

    Returns
    -------
    fig: Figure
        Plotly box plot figure of highest expr genes.
    """
     # compute the percentage of each gene per cell
    norm_dict = sc.pp.normalize_total(adata, target_sum=100, inplace=False)

    # identify the genes with the highest mean
    if issparse(norm_dict['X']):
        mean_percent = norm_dict['X'].mean(axis=0).A1
        top_idx = np.argsort(mean_percent)[::-1][:n_top]
        counts_top_genes = norm_dict['X'][:, top_idx].A
    else:
        mean_percent = norm_dict['X'].mean(axis=0)
        top_idx = np.argsort(mean_percent)[::-1][:n_top]
        counts_top_genes = norm_dict['X'][:, top_idx]
    columns = (
        adata.var_names[top_idx]
        if gene_symbols is None
        else adata.var[gene_symbols][top_idx]
    )
    counts_top_genes = pd.DataFrame(
        counts_top_genes, index=adata.obs_names, columns=columns
    )

    print(counts_top_genes)

    fig = go.Figure()
    # Use x instead of y argument for horizontal plot
    for column in reversed(counts_top_genes.columns):
        fig.add_trace(go.Box(x=counts_top_genes[column], name=column, marker_size=2, line_width=1, boxpoints=False))

    fig.update_layout(width=400, showlegend=False, title=f"Top {n_top} most expressed genes", yaxis={'dtick':1})

    return fig





    