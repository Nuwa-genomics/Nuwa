from typing import Iterable, Optional
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


def plot_top_ranked_genes(adata, cluster_name, method, n_rows=5, height = None, width = None):

    sc.tl.rank_genes_groups(adata, groupby=cluster_name, method=method, key_added=method)

    names = pd.DataFrame(adata.uns[method]['names']).head(n_rows)
    pvalues = pd.DataFrame(adata.uns[method]['pvals']).head(n_rows)
    clusters = names.columns

    #reshape for heatmap. Heatmap plots from botttom up so reverse
    names = reversed(names.values.reshape(n_rows,len(clusters)))
    pvalues = reversed(pvalues.values.reshape(n_rows, len(clusters)))

    fig = go.Figure(data=go.Heatmap(
                        z=pvalues,
                        text=names,
                        x=clusters,
                        y=np.arange(0, n_rows - 1),
                        hovertemplate = "%{text}: <br>P score: %{z} </br> Cluster: %{x}",
                        texttemplate="%{text}",
                        textfont={"size":18}))


    fig.layout.width = width
    fig.layout.height = height 
    fig.update_layout(title="Top ranked genes with p values",
                    yaxis={"title": 'Row'},
                    xaxis={"title": cluster_name})

    return fig


def plot_doubletdetection_threshold_heatmap(
    clf,
    show=False,
    save=None,
    log10=True,
    log_p_grid=None,
    voter_grid=None,
    height=None,
    width=None,
    v_step=2,
    p_step=5,
):
    """Produce a plot showing number of cells called doublet across
       various thresholds

    Args:
        clf (BoostClassifier object): Fitted classifier
        show (bool, optional): If True, runs plt.show()
        save (str, optional): If provided, the figure is saved to this
            filepath.
        log10 (bool, optional): Use log 10 if true, natural log if false.
        log_p_grid (ndarray, optional): log p-value thresholds to use.
            Defaults to np.arange(-100, -1). log base decided by log10
        voter_grid (ndarray, optional): Voting thresholds to use. Defaults to
            np.arange(0.3, 1.0, 0.05).
        p_step (int, optional): number of xlabels to skip in plot
        v_step (int, optional): number of ylabels to skip in plot


    Returns:
        matplotlib figure
    """
    # Ignore numpy complaining about np.nan comparisons
    with np.errstate(invalid="ignore"):
        all_log_p_values_ = np.copy(clf.all_log_p_values_)
        if log10:
            all_log_p_values_ /= np.log(10)
        if log_p_grid is None:
            log_p_grid = np.arange(-100, -1)
        if voter_grid is None:
            voter_grid = np.arange(0.3, 1.0, 0.05)
        doubs_per_t = np.zeros((len(voter_grid), len(log_p_grid)))
        for i in range(len(voter_grid)):
            for j in range(len(log_p_grid)):
                voting_average = np.mean(
                    np.ma.masked_invalid(all_log_p_values_) <= log_p_grid[j], axis=0
                )
                labels = np.ma.filled((voting_average >= voter_grid[i]).astype(float), np.nan)
                doubs_per_t[i, j] = np.nansum(labels)

    df = pd.DataFrame(data=doubs_per_t, index=voter_grid, columns=log_p_grid)
    fig = px.imshow(df, height=height, width=width, aspect="auto", labels=dict(x="log10 p-value", y="Voting threshold", color="Predicted doublets"), title="Threshold Diagnostics")
    fig.update_coloraxes(colorbar_title_side="right")
    return fig

def get_color_embeddings_from_key(key: str, adata: AnnData):
    if key in adata.obs_keys():
        return adata.obs[key].values
    else:
        return adata.to_df()[key].values
    