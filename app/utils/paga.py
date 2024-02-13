from types import MappingProxyType
import scanpy as sc
from anndata import AnnData
from typing import Optional, Union, Mapping, Sequence, Any, List
import numpy as np
from pandas.api.types import CategoricalDtype
from matplotlib import pyplot as pl, rcParams, ticker
from matplotlib import patheffects
from matplotlib.axes import Axes
from matplotlib.colors import is_color_like, Colormap
from scipy.sparse import issparse
from sklearn.utils import check_random_state
import warnings
import collections.abc as cabc
from pathlib import Path
import scipy
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
import networkx as nx

from scanpy.plotting import _utils
from scanpy.plotting._utils import matrix, _IGraphLayout, _FontWeight, _FontSize
from scanpy import _utils as _sc_utils, logging as logg
from scanpy._settings import settings


def draw_network_graph(adjacency, pos, colors, title, names, height=None, width=None):

    G = nx.Graph()

    # add nodes to networkx graph
    for id, p in enumerate(pos):
        G.add_node(id, pos=p)
        
    # add edges to networkx graph
    cx = adjacency.tocoo()    
    for v0,v1,w in zip(cx.row, cx.col, cx.data):
        G.add_edge(v0, v1, weight=w)

    edge_x = []
    edge_y = []

    # add edge positions
    for edge in G.edges():
        x0, y0 = G.nodes[edge[0]]['pos']
        x1, y1 = G.nodes[edge[1]]['pos']
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)

    # create weights list
    weights = []
    for i in G.edges.data("weight", default=1):
        (_, _, w) = i
        weights.append(w)

    # reshape into x, y and z positions
    edge_x = np.array(edge_x).reshape(int(len(edge_x)/3), 3)
    edge_y = np.array(edge_y).reshape(int(len(edge_y)/3), 3)

    # create edge traces with weights
    edge_trace = [go.Scatter(
    x=edge_x[i], y=edge_y[i],
    showlegend=False if isinstance(colors[0], str) else True,
    line=dict(width=weights[i]*3, color='#888'),
    hoverinfo='skip',
    mode='lines') for i, _ in enumerate(edge_x)]

    # create x and y node lists
    node_x = []
    node_y = []
    for node in G.nodes():
        x, y = G.nodes[node]['pos']
        node_x.append(x)
        node_y.append(y)

    # create node trace
    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers+text',
        text=names,
        textposition="top center",
        hoverinfo='skip',
        marker=dict(
            showscale=(not isinstance(colors[0], str)),
            # colorscale options
            #'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu' |
            #'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet' |
            #'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
            colorscale='YlGnBu',
            reversescale=True,
            color=[],
            size=12,
            colorbar=dict(
                thickness=15,
                title='Expression level',
                xanchor='left',
                titleside='right'
            ),
            line_width=2))


    node_hover_template = []

    for node, adjacencies in enumerate(G.adjacency()):

        node_hover_template.append(f'<b>{names[node]}</b>')

    # set colours and legend for categorical vs continuous data
    if isinstance(colors[0], str):
        node_trace.marker.color = colors
        node_trace.showlegend = False
    else:
        coloraxes_colors = [c for c in colors]
        node_trace.marker.color = coloraxes_colors

    # add hovertemplate
    node_trace.hovertemplate = node_hover_template

    # create figure of edge traces only
    fig = go.Figure(data=edge_trace,
                layout=go.Layout(
                    title=title,
                    titlefont_size=16,
                    showlegend=False,
                    hovermode='closest',
                    height=height,
                    width=width,
                    margin=dict(b=20,l=5,r=5,t=40),
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                    )
    
    # now add node trace
    fig.add_trace(node_trace)
    
    return fig



def _compute_pos(
    adjacency_solid,
    layout=None,
    random_state=0,
    init_pos=None,
    adj_tree=None,
    root=0,
    layout_kwds: Mapping[str, Any] = MappingProxyType({}),
):
    import random
    import networkx as nx

    random_state = check_random_state(random_state)

    nx_g_solid = nx.Graph(adjacency_solid)
    if layout is None:
        layout = 'fr'
    if layout == 'fa':
        try:
            from fa2 import ForceAtlas2
        except ImportError:
            logg.warning(
                "Package 'fa2' is not installed, falling back to layout 'fr'."
                'To use the faster and better ForceAtlas2 layout, '
                "install package 'fa2' (`pip install fa2`)."
            )
            layout = 'fr'
    if layout == 'fa':
        # np.random.seed(random_state)
        if init_pos is None:
            init_coords = random_state.random_sample((adjacency_solid.shape[0], 2))
        else:
            init_coords = init_pos.copy()
        forceatlas2 = ForceAtlas2(
            # Behavior alternatives
            outboundAttractionDistribution=False,  # Dissuade hubs
            linLogMode=False,  # NOT IMPLEMENTED
            adjustSizes=False,  # Prevent overlap (NOT IMPLEMENTED)
            edgeWeightInfluence=1.0,
            # Performance
            jitterTolerance=1.0,  # Tolerance
            barnesHutOptimize=True,
            barnesHutTheta=1.2,
            multiThreaded=False,  # NOT IMPLEMENTED
            # Tuning
            scalingRatio=2.0,
            strongGravityMode=False,
            gravity=1.0,
            # Log
            verbose=False,
        )
        if 'maxiter' in layout_kwds:
            iterations = layout_kwds['maxiter']
        elif 'iterations' in layout_kwds:
            iterations = layout_kwds['iterations']
        else:
            iterations = 500
        pos_list = forceatlas2.forceatlas2(
            adjacency_solid, pos=init_coords, iterations=iterations
        )
        pos = {n: [p[0], -p[1]] for n, p in enumerate(pos_list)}
    elif layout == 'eq_tree':
        nx_g_tree = nx.Graph(adj_tree)
        pos = _utils.hierarchy_pos(nx_g_tree, root)
        if len(pos) < adjacency_solid.shape[0]:
            raise ValueError(
                'This is a forest and not a single tree. '
                'Try another `layout`, e.g., {\'fr\'}.'
            )
    else:
        # igraph layouts
        random.seed(random_state.bytes(8))
        g = _sc_utils.get_igraph_from_adjacency(adjacency_solid)
        if 'rt' in layout:
            g_tree = _sc_utils.get_igraph_from_adjacency(adj_tree)
            pos_list = g_tree.layout(
                layout, root=root if isinstance(root, list) else [root]
            ).coords
        elif layout == 'circle':
            pos_list = g.layout(layout).coords
        else:
            # I don't know why this is necessary
            # np.random.seed(random_state)
            if init_pos is None:
                init_coords = random_state.random_sample(
                    (adjacency_solid.shape[0], 2)
                ).tolist()
            else:
                init_pos = init_pos.copy()
                # this is a super-weird hack that is necessary as igraph’s
                # layout function seems to do some strange stuff here
                init_pos[:, 1] *= -1
                init_coords = init_pos.tolist()
            try:
                pos_list = g.layout(
                    layout, seed=init_coords, weights='weight', **layout_kwds
                ).coords
            except AttributeError:  # hack for empty graphs...
                pos_list = g.layout(layout, seed=init_coords, **layout_kwds).coords
        pos = {n: [p[0], -p[1]] for n, p in enumerate(pos_list)}
    if len(pos) == 1:
        pos[0] = (0.5, 0.5)
    pos_array = np.array([pos[n] for count, n in enumerate(nx_g_solid)])
    return pos_array

def compute_paga(
    adata: AnnData,
    threshold: Optional[float] = None,
    color: Optional[Union[str, Mapping[Union[str, int], Mapping[Any, float]]]] = None,
    layout: Optional[_IGraphLayout] = None,
    layout_kwds: Mapping[str, Any] = MappingProxyType({}),
    init_pos: Optional[np.ndarray] = None,
    root: Union[int, str, Sequence[int], None] = 0,
    labels: Union[str, Sequence[str], Mapping[str, str], None] = None,
    single_component: bool = False,
    solid_edges: str = 'connectivities',
    dashed_edges: Optional[str] = None,
    transitions: Optional[str] = None,
    fontsize: Optional[int] = None,
    fontweight: str = 'bold',
    fontoutline: Optional[int] = None,
    text_kwds: Mapping[str, Any] = MappingProxyType({}),
    node_size_scale: float = 1.0,
    node_size_power: float = 0.5,
    edge_width_scale: float = 1.0,
    min_edge_width: Optional[float] = None,
    max_edge_width: Optional[float] = None,
    arrowsize: int = 30,
    title: Optional[str] = None,
    left_margin: float = 0.01,
    random_state: Optional[int] = 0,
    pos: Union[np.ndarray, str, Path, None] = None,
    normalize_to_color: bool = False,
    cmap: Union[str, Colormap] = None,
    cax: Optional[Axes] = None,
    colorbar=None,  # TODO: this seems to be unused
    cb_kwds: Mapping[str, Any] = MappingProxyType({}),
    frameon: Optional[bool] = None,
    add_pos: bool = True,
    export_to_gexf: bool = False,
    use_raw: bool = True,
    colors=None,  # backwards compat
    groups=None,  # backwards compat
    height=None,
    width=None
):
    """\
    Plot the PAGA graph through thresholding low-connectivity edges.

    Compute a coarse-grained layout of the data. Reuse this by passing
    `init_pos='paga'` to :func:`~scanpy.tl.umap` or
    :func:`~scanpy.tl.draw_graph` and obtain embeddings with more meaningful
    global topology [Wolf19]_.

    This uses ForceAtlas2 or igraph's layout algorithms for most layouts [Csardi06]_.

    Parameters
    ----------
    adata
        Annotated data matrix.
    threshold
        Do not draw edges for weights below this threshold. Set to 0 if you want
        all edges. Discarding low-connectivity edges helps in getting a much
        clearer picture of the graph.
    color
        Gene name or `obs` annotation defining the node colors.
        Also plots the degree of the abstracted graph when
        passing {`'degree_dashed'`, `'degree_solid'`}.

        Can be also used to visualize pie chart at each node in the following form:
        `{<group name or index>: {<color>: <fraction>, ...}, ...}`. If the fractions
        do not sum to 1, a new category called `'rest'` colored grey will be created.
    labels
        The node labels. If `None`, this defaults to the group labels stored in
        the categorical for which :func:`~scanpy.tl.paga` has been computed.
    pos
        Two-column array-like storing the x and y coordinates for drawing.
        Otherwise, path to a `.gdf` file that has been exported from Gephi or
        a similar graph visualization software.
    layout
        Plotting layout that computes positions.
        `'fa'` stands for “ForceAtlas2”,
        `'fr'` stands for “Fruchterman-Reingold”,
        `'rt'` stands for “Reingold-Tilford”,
        `'eq_tree'` stands for “eqally spaced tree”.
        All but `'fa'` and `'eq_tree'` are igraph layouts.
        All other igraph layouts are also permitted.
        See also parameter `pos` and :func:`~scanpy.tl.draw_graph`.
    layout_kwds
        Keywords for the layout.
    init_pos
        Two-column array storing the x and y coordinates for initializing the
        layout.
    random_state
        For layouts with random initialization like `'fr'`, change this to use
        different intial states for the optimization. If `None`, the initial
        state is not reproducible.
    root
        If choosing a tree layout, this is the index of the root node or a list
        of root node indices. If this is a non-empty vector then the supplied
        node IDs are used as the roots of the trees (or a single tree if the
        graph is connected). If this is `None` or an empty list, the root
        vertices are automatically calculated based on topological sorting.
    transitions
        Key for `.uns['paga']` that specifies the matrix that stores the
        arrows, for instance `'transitions_confidence'`.
    solid_edges
        Key for `.uns['paga']` that specifies the matrix that stores the edges
        to be drawn solid black.
    dashed_edges
        Key for `.uns['paga']` that specifies the matrix that stores the edges
        to be drawn dashed grey. If `None`, no dashed edges are drawn.
    single_component
        Restrict to largest connected component.
    fontsize
        Font size for node labels.
    fontoutline
        Width of the white outline around fonts.
    text_kwds
        Keywords for :meth:`~matplotlib.axes.Axes.text`.
    node_size_scale
        Increase or decrease the size of the nodes.
    node_size_power
        The power with which groups sizes influence the radius of the nodes.
    edge_width_scale
        Edge with scale in units of `rcParams['lines.linewidth']`.
    min_edge_width
        Min width of solid edges.
    max_edge_width
        Max width of solid and dashed edges.
    arrowsize
       For directed graphs, choose the size of the arrow head head's length and
       width. See :py:class: `matplotlib.patches.FancyArrowPatch` for attribute
       `mutation_scale` for more info.
    export_to_gexf
        Export to gexf format to be read by graph visualization programs such as
        Gephi.
    normalize_to_color
        Whether to normalize categorical plots to `color` or the underlying
        grouping.
    cmap
        The color map.
    cax
        A matplotlib axes object for a potential colorbar.
    cb_kwds
        Keyword arguments for :class:`~matplotlib.colorbar.Colorbar`,
        for instance, `ticks`.
    add_pos
        Add the positions to `adata.uns['paga']`.
    title
        Provide a title.
    frameon
        Draw a frame around the PAGA graph.
    plot
        If `False`, do not create the figure, simply compute the layout.
    save
        If `True` or a `str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on \\{`'.pdf'`, `'.png'`, `'.svg'`\\}.
    ax
        A matplotlib axes object.

    Returns
    -------
    If `show==False`, one or more :class:`~matplotlib.axes.Axes` objects.
    Adds `'pos'` to `adata.uns['paga']` if `add_pos` is `True`.

    Examples
    --------

    .. plot::
        :context: close-figs

        import scanpy as sc
        adata = sc.datasets.pbmc3k_processed()
        sc.tl.paga(adata, groups='louvain')
        sc.pl.paga(adata)

    You can increase node and edge sizes by specifying additional arguments.

    .. plot::
        :context: close-figs

        sc.pl.paga(adata, node_size_scale=10, edge_width_scale=2)

    Notes
    -----
    When initializing the positions, note that – for some reason – igraph
    mirrors coordinates along the x axis... that is, you should increase the
    `maxiter` parameter by 1 if the layout is flipped.

    .. currentmodule:: scanpy

    See also
    --------
    tl.paga
    pl.paga_compare
    pl.paga_path
    """

    if groups is not None:  # backwards compat
        labels = groups
        logg.warning('`groups` is deprecated in `pl.paga`: use `labels` instead')
    if colors is None:
        colors = color

    groups_key = adata.uns['paga']['groups']

    def is_flat(x):
        has_one_per_category = isinstance(x, cabc.Collection) and len(x) == len(
            adata.obs[groups_key].cat.categories
        )
        return has_one_per_category or x is None or isinstance(x, str)

    if isinstance(colors, cabc.Mapping) and isinstance(
        colors[next(iter(colors))], cabc.Mapping
    ):
        # handle paga pie, remap string keys to integers
        names_to_ixs = {
            n: i for i, n in enumerate(adata.obs[groups_key].cat.categories)
        }
        colors = {names_to_ixs.get(n, n): v for n, v in colors.items()}
    if is_flat(colors):
        colors = [colors]

    if frameon is None:
        frameon = settings._frameon
    # labels is a list that contains no lists
    if is_flat(labels):
        labels = [labels for _ in range(len(colors))]

    if title is None and len(colors) > 1:
        title = [c for c in colors]
    elif isinstance(title, str):
        title = [title for c in colors]
    elif title is None:
        title = [None for c in colors]


    if colorbar is None:
        var_names = adata.var_names if adata.raw is None else adata.raw.var_names
        colorbars = [
            (
                (c in adata.obs_keys() and adata.obs[c].dtype.name != 'category')
                or (c in var_names)
            )
            for c in colors
        ]
    else:
        colorbars = [False for _ in colors]


    if isinstance(root, str):
        if root not in labels:
            raise ValueError(
                'If `root` is a string, '
                f'it needs to be one of {labels} not {root!r}.'
            )
        root = list(labels).index(root)
    if isinstance(root, cabc.Sequence) and root[0] in labels:
        root = [list(labels).index(r) for r in root]


    # define the adjacency matrices
    adjacency_solid = adata.uns['paga'][solid_edges].copy()
    adjacency_dashed = None
    if threshold is None:
        threshold = 0.01  # default threshold
    if threshold > 0:
        adjacency_solid.data[adjacency_solid.data < threshold] = 0
        adjacency_solid.eliminate_zeros()
    if dashed_edges is not None:
        adjacency_dashed = adata.uns['paga'][dashed_edges].copy()
        if threshold > 0:
            adjacency_dashed.data[adjacency_dashed.data < threshold] = 0
            adjacency_dashed.eliminate_zeros()

    # compute positions
    adj_tree = None
    if layout in {'rt', 'rt_circular', 'eq_tree'}:
        adj_tree = adata.uns['paga']['connectivities_tree']
    pos = _compute_pos(
        adjacency_solid,
        layout=layout,
        random_state=random_state,
        init_pos=init_pos,
        layout_kwds=layout_kwds,
        adj_tree=adj_tree,
        root=root,
    )


    figs = []

    for icolor, c in enumerate(colors):

        colors = _paga_graph_colors(
            adata,
            colors=colors if isinstance(colors, cabc.Mapping) else c,
            solid_edges=solid_edges,
            dashed_edges=dashed_edges,
            transitions=transitions,
            threshold=threshold,
            adjacency_solid=adjacency_solid,
            adjacency_dashed=adjacency_dashed,
            root=root,
            labels=labels[icolor],
            fontsize=fontsize,
            fontweight=fontweight,
            fontoutline=fontoutline,
            text_kwds=text_kwds,
            node_size_scale=node_size_scale,
            node_size_power=node_size_power,
            edge_width_scale=edge_width_scale,
            min_edge_width=min_edge_width,
            max_edge_width=max_edge_width,
            normalize_to_color=normalize_to_color,
            frameon=frameon,
            cmap=cmap,
            colorbar=colorbars[icolor],
            cb_kwds=cb_kwds,
            use_raw=use_raw,
            title=title[icolor],
            export_to_gexf=export_to_gexf,
            single_component=single_component,
            arrowsize=arrowsize,
            pos=pos,
        )

        fig = draw_network_graph(adjacency=adjacency_solid, 
            pos=pos, 
            colors=colors, 
            title=title[icolor], 
            names=adata.obs[groups_key].cat.categories,
            height=height,
            width=width
        )
        figs.append(fig)

    return figs


def _paga_graph_colors(
    adata,
    solid_edges=None,
    dashed_edges=None,
    adjacency_solid=None,
    adjacency_dashed=None,
    transitions=None,
    threshold=None,
    root=0,
    colors=None,
    labels=None,
    fontsize=None,
    fontweight=None,
    fontoutline=None,
    text_kwds: Mapping[str, Any] = MappingProxyType({}),
    node_size_scale=1.0,
    node_size_power=0.5,
    edge_width_scale=1.0,
    normalize_to_color='reference',
    title=None,
    pos=None,
    cmap=None,
    frameon=True,
    min_edge_width=None,
    max_edge_width=None,
    export_to_gexf=False,
    colorbar=None,
    use_raw=True,
    cb_kwds: Mapping[str, Any] = MappingProxyType({}),
    single_component=False,
    arrowsize=30,
):
    import networkx as nx

    node_labels = labels  # rename for clarity
    if (
        node_labels is not None
        and isinstance(node_labels, str)
        and node_labels != adata.uns['paga']['groups']
    ):
        raise ValueError(
            'Provide a list of group labels for the PAGA groups {}, not {}.'.format(
                adata.uns['paga']['groups'], node_labels
            )
        )
    groups_key = adata.uns['paga']['groups']
    if node_labels is None:
        node_labels = adata.obs[groups_key].cat.categories

    if (colors is None or colors == groups_key) and groups_key is not None:
        if groups_key + '_colors' not in adata.uns or len(
            adata.obs[groups_key].cat.categories
        ) != len(adata.uns[groups_key + '_colors']):
            _utils.add_colors_for_categorical_sample_annotation(adata, groups_key)
        colors = adata.uns[groups_key + '_colors']
        for iname, name in enumerate(adata.obs[groups_key].cat.categories):
            if name in settings.categories_to_ignore:
                colors[iname] = 'grey'

    nx_g_solid = nx.Graph(adjacency_solid)
    if dashed_edges is not None:
        nx_g_dashed = nx.Graph(adjacency_dashed)

    # convert pos to array and dict
    if not isinstance(pos, (Path, str)):
        pos_array = pos
    else:
        pos = Path(pos)
        if pos.suffix != '.gdf':
            raise ValueError(
                'Currently only supporting reading positions from .gdf files. '
                'Consider generating them using, for instance, Gephi.'
            )
        s = ''  # read the node definition from the file
        with pos.open() as f:
            f.readline()
            for line in f:
                if line.startswith('edgedef>'):
                    break
                s += line
        from io import StringIO

        df = pd.read_csv(StringIO(s), header=-1)
        pos_array = df[[4, 5]].values

    # convert to dictionary
    pos = {n: [p[0], p[1]] for n, p in enumerate(pos_array)}

    # uniform color
    if isinstance(colors, str) and is_color_like(colors):
        colors = [colors for c in range(len(node_labels))]

    # color degree of the graph
    if isinstance(colors, str) and colors.startswith('degree'):
        # see also tools.paga.paga_degrees
        if colors == 'degree_dashed':
            colors = [d for _, d in nx_g_dashed.degree(weight='weight')]
        elif colors == 'degree_solid':
            colors = [d for _, d in nx_g_solid.degree(weight='weight')]
        else:
            raise ValueError('`degree` either "degree_dashed" or "degree_solid".')
        colors = (np.array(colors) - np.min(colors)) / (np.max(colors) - np.min(colors))

    # plot gene expression
    var_names = adata.var_names if adata.raw is None else adata.raw.var_names
    if isinstance(colors, str) and colors in var_names:
        x_color = []
        cats = adata.obs[groups_key].cat.categories
        for icat, cat in enumerate(cats):
            subset = (cat == adata.obs[groups_key]).values
            if adata.raw is not None and use_raw:
                adata_gene = adata.raw[:, colors]
            else:
                adata_gene = adata[:, colors]
            x_color.append(np.mean(adata_gene.X[subset]))
        colors = x_color

    # plot continuous annotation
    if (
        isinstance(colors, str)
        and colors in adata.obs
        and not isinstance(adata.obs[colors].dtype, CategoricalDtype)
    ):
        x_color = []
        cats = adata.obs[groups_key].cat.categories
        for icat, cat in enumerate(cats):
            subset = (cat == adata.obs[groups_key]).values
            x_color.append(adata.obs.loc[subset, colors].mean())
        colors = x_color

    # plot categorical annotation
    if (
        isinstance(colors, str)
        and colors in adata.obs
        and isinstance(adata.obs[colors].dtype, CategoricalDtype)
    ):
        asso_names, asso_matrix = _sc_utils.compute_association_matrix_of_groups(
            adata,
            prediction=groups_key,
            reference=colors,
            normalization='reference' if normalize_to_color else 'prediction',
        )
        _utils.add_colors_for_categorical_sample_annotation(adata, colors)
        asso_colors = _sc_utils.get_associated_colors_of_groups(
            adata.uns[colors + '_colors'], asso_matrix
        )
        colors = asso_colors

    if len(colors) != len(node_labels):
        raise ValueError(
            f'Expected `colors` to be of length `{len(node_labels)}`, '
            f'found `{len(colors)}`.'
        )

    # count number of connected components
    n_components, labels = scipy.sparse.csgraph.connected_components(adjacency_solid)
    if n_components > 1 and not single_component:
        logg.debug(
            'Graph has more than a single connected component. '
            'To restrict to this component, pass `single_component=True`.'
        )
    if n_components > 1 and single_component:
        component_sizes = np.bincount(labels)
        largest_component = np.where(component_sizes == component_sizes.max())[0][0]
        adjacency_solid = adjacency_solid.tocsr()[labels == largest_component, :]
        adjacency_solid = adjacency_solid.tocsc()[:, labels == largest_component]
        colors = np.array(colors)[labels == largest_component]
        node_labels = np.array(node_labels)[labels == largest_component]
        cats_dropped = (
            adata.obs[groups_key].cat.categories[labels != largest_component].tolist()
        )
        logg.info(
            'Restricting graph to largest connected component by dropping categories\n'
            f'{cats_dropped}'
        )
        nx_g_solid = nx.Graph(adjacency_solid)
        if dashed_edges is not None:
            raise ValueError('`single_component` only if `dashed_edges` is `None`.')

    # edge widths
    base_edge_width = edge_width_scale * 5 * rcParams['lines.linewidth']

    # draw dashed edges
    if dashed_edges is not None:
        widths = [x[-1]['weight'] for x in nx_g_dashed.edges(data=True)]
        widths = base_edge_width * np.array(widths)
        if max_edge_width is not None:
            widths = np.clip(widths, None, max_edge_width)

    # draw solid edges
    if transitions is None:
        widths = [x[-1]['weight'] for x in nx_g_solid.edges(data=True)]
        widths = base_edge_width * np.array(widths)
        if min_edge_width is not None or max_edge_width is not None:
            widths = np.clip(widths, min_edge_width, max_edge_width)
    # draw directed edges
    else:
        adjacency_transitions = adata.uns['paga'][transitions].copy()
        if threshold is None:
            threshold = 0.01
        adjacency_transitions.data[adjacency_transitions.data < threshold] = 0
        adjacency_transitions.eliminate_zeros()
        g_dir = nx.DiGraph(adjacency_transitions.T)
        widths = [x[-1]['weight'] for x in g_dir.edges(data=True)]
        widths = base_edge_width * np.array(widths)
        if min_edge_width is not None or max_edge_width is not None:
            widths = np.clip(widths, min_edge_width, max_edge_width)

    if export_to_gexf:
        if isinstance(colors[0], tuple):
            from matplotlib.colors import rgb2hex

            colors = [rgb2hex(c) for c in colors]
        for count, n in enumerate(nx_g_solid.nodes()):
            nx_g_solid.node[count]['label'] = str(node_labels[count])
            nx_g_solid.node[count]['color'] = str(colors[count])
            nx_g_solid.node[count]['viz'] = dict(
                position=dict(
                    x=1000 * pos[count][0],
                    y=1000 * pos[count][1],
                    z=0,
                )
            )
        filename = settings.writedir / 'paga_graph.gexf'
        logg.warning(f'exporting to {filename}')
        settings.writedir.mkdir(parents=True, exist_ok=True)
        nx.write_gexf(nx_g_solid, settings.writedir / 'paga_graph.gexf')


    # groups sizes
    if groups_key is not None and groups_key + '_sizes' in adata.uns:
        groups_sizes = adata.uns[groups_key + '_sizes']
    else:
        groups_sizes = np.ones(len(node_labels))
    base_scale_scatter = 2000
    base_pie_size = (
        base_scale_scatter / (np.sqrt(adjacency_solid.shape[0]) + 10) * node_size_scale
    )
    median_group_size = np.median(groups_sizes)
    groups_sizes = base_pie_size * np.power(
        groups_sizes / median_group_size, node_size_power
    )

    if fontsize is None:
        fontsize = rcParams['legend.fontsize']
    if fontoutline is not None:
        text_kwds = dict(text_kwds)
        text_kwds['path_effects'] = [
            patheffects.withStroke(linewidth=fontoutline, foreground='w')
        ]
    # usual scatter plot
    if not isinstance(colors[0], cabc.Mapping):
        n_groups = len(pos_array)


    # else pie chart plot
    else:
        for ix, (xx, yy) in enumerate(zip(pos_array[:, 0], pos_array[:, 1])):
            if not isinstance(colors[ix], cabc.Mapping):
                raise ValueError(
                    f'{colors[ix]} is neither a dict of valid '
                    'matplotlib colors nor a valid matplotlib color.'
                )
            color_single = colors[ix].keys()
            fracs = [colors[ix][c] for c in color_single]
            total = sum(fracs)

            if total < 1:
                color_single = list(color_single)
                color_single.append('grey')
                fracs.append(1 - sum(fracs))
            elif not np.isclose(total, 1):
                raise ValueError(
                    f'Expected fractions for node `{ix}` to be '
                    f'close to 1, found `{total}`.'
                )

            cumsum = np.cumsum(fracs)
            cumsum = cumsum / cumsum[-1]
            cumsum = [0] + cumsum.tolist()

            for r1, r2, color in zip(cumsum[:-1], cumsum[1:], color_single):
                angles = np.linspace(2 * np.pi * r1, 2 * np.pi * r2, 20)
                x = [0] + np.cos(angles).tolist()
                y = [0] + np.sin(angles).tolist()

                xy = np.column_stack([x, y])
                s = np.abs(xy).max()

    
    return colors


