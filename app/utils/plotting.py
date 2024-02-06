from typing import Optional
import plotly.graph_objects as go
import numpy as np
import scanpy as sc
from anndata import AnnData
from scipy.sparse import issparse
import pandas as pd
import streamlit as st
from squidpy._constants._pkg_constants import Key
from anndata import AnnData
import scanpy as sc
import plotly.express as px
import streamlit as st
from squidpy._constants._constants import RipleyStat
from typing import (
    TYPE_CHECKING,
    Literal,
    Any,
    Mapping,
    Sequence,
)


from pathlib import Path
from types import MappingProxyType
import matplotlib.pyplot as plt
import seaborn as sns

from squidpy.gr._utils import (
    _assert_categorical_obs,
    _assert_non_empty_sequence,
    _get_valid_values,
)
from squidpy.pl._color_utils import Palette_t, _get_palette

def get_data(adata: AnnData, cluster_key: str, func_name: str, attr: str = "uns", **kwargs: Any) -> Any:
    key = getattr(Key.uns, func_name)(cluster_key, **kwargs)

    try:
        if attr == "uns":
            return adata.uns[key]
        elif attr == "obsm":
            return adata.obsm[key]
        else:
            raise ValueError(f"attr must be either 'uns' or 'obsm', got {attr}")
    except KeyError:
        raise KeyError(
            f"Unable to get the data from `adata.uns[{key!r}]`. "
            f"Please run `squidpy.gr.{func_name}(..., cluster_key={cluster_key!r})` first."
        ) from None

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


def plot_ripley(adata: AnnData, cluster_key: str, plot_sims=True, mode: Literal["F", "G", "L"] = "F", height=None, width=None):

    res = get_data(adata, cluster_key=cluster_key, func_name="ripley", mode=mode)

    mode = RipleyStat(mode)

    if TYPE_CHECKING:
        assert isinstance(mode, RipleyStat)
    
    if plot_sims:
        fig = px.line(data_frame=res["sims_stat"], x="bins", y="stats", title=f"Ripley's {mode.s}", height=height, width=width)
        return fig 
    else:
        fig = px.line(data_frame=res[f"{mode.s}_stat"], x="bins", y="stats", color=cluster_key, title=f"Ripley's {mode.s}", height=height, width=width)
        return fig

def plot_co_occurrence(
    adata: AnnData,
    cluster_key: str,
    clusters: str,
    height = None,
    width = None,
    **kwargs: Any,
) -> None:
    """
    Plot co-occurrence probability ratio for each cluster.
    """
    _assert_categorical_obs(adata, key=cluster_key)
    occurrence_data = get_data(adata, cluster_key=cluster_key, func_name="co_occurrence")

    out = occurrence_data["occ"]
    interval = occurrence_data["interval"][1:]
    categories = adata.obs[cluster_key].cat.categories

    clusters = categories if clusters is None else clusters
    clusters = _assert_non_empty_sequence(clusters, name="clusters")
    clusters = sorted(_get_valid_values(clusters, categories))

    fig, axs = plt.subplots(
        1,
        len(clusters)
    )
    axs = np.ravel(axs)  # make into iterable
    figs = []

    for g, ax in zip(clusters, axs):
        idx = np.where(categories == g)[0][0]
        df = pd.DataFrame(out[idx, :, :].T, columns=categories).melt(var_name=cluster_key, value_name="probability")
        df["distance"] = np.tile(interval, len(categories))

        fig = px.line(data_frame=df, x="distance", y="probability", color=cluster_key, height=height, width=width)
        figs.append(fig)

    return figs


def plot_centrality_scores(
    adata: AnnData,
    cluster_key: str,
    score: str | Sequence[str] | None = None,
    legend_kwargs: Mapping[str, Any] = MappingProxyType({}),
    palette: Palette_t = None,
    size = 13,
    height = None,
    width = None,
    **kwargs: Any,
) -> None:
    """
    Plot centrality scores.

    The centrality scores are computed by :func:`squidpy.gr.centrality_scores`.

    Parameters
    ----------
    %(adata)s
    %(cluster_key)s
    score
        Whether to plot all scores or only selected ones.
    legend_kwargs
        Keyword arguments for :func:`matplotlib.pyplot.legend`.
    %(cat_plotting)s

    Returns
    -------
    %(plotting_returns)s
    """
    _assert_categorical_obs(adata, key=cluster_key)
    df = get_data(adata, cluster_key=cluster_key, func_name="centrality_scores").copy()

    legend_kwargs = dict(legend_kwargs)
    if "loc" not in legend_kwargs:
        legend_kwargs["loc"] = "center left"
        legend_kwargs.setdefault("bbox_to_anchor", (1, 0.5))

    scores = df.columns.values
    df[cluster_key] = df.index.values

    clusters = adata.obs[cluster_key].cat.categories
    palette = _get_palette(adata, cluster_key=cluster_key, categories=clusters, palette=palette)

    score = scores if score is None else score
    score = _assert_non_empty_sequence(score, name="centrality scores")
    score = sorted(_get_valid_values(score, scores))

    figs = []
    for g in score:
        fig = px.scatter(data_frame=df, x=g, y=cluster_key, color=cluster_key, height=height, width=width)
        fig.update_traces(marker=dict(size=size), selector=dict(mode='markers'))
        fig.update_layout(showlegend=False)
        figs.append(fig)

    return figs



    