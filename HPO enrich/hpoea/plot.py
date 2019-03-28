from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, HoverTool
import networkx as nx

from hpoea.utils.log import get_logger

log = get_logger(__name__)


def dot_plot(enrich_table, figsize=None, show_items=10, x="gene_ratio", **kwargs):
    """ Draw the enrichment dot chart.

    Args
    ----
    enrich_table : pandas.DataFrame
        Enrichment result table, see `hpoea.enrich.GSEA`.
    figsize : tuple, optional
        The size of output figure, (width, height).
    show_items : int, optional
        Number of items to show. default 10.
    x : str, optional
        The x-axis coresponding value.
        ['gene_ratio' | 'background_ratio' | 'odd_ratio' | 'pvalue' | 'padj']
    """
    if figsize:
        width, height = figsize
    else:  # set default figsize
        grid_height = 50
        height = grid_height * show_items
        width = 500

    if x in ['pvalue', 'padj']:
        df = enrich_table.sort_values(by=x, ascending=True)
    else:
        df = enrich_table.sort_values(by=x, ascending=False)
    df = df.iloc[:show_items, :]

    source = ColumnDataSource(data=df)
    p = figure(plot_width=width, plot_height=height,
               y_range=list(df.HPO_term_name))

    kwargs.setdefault("color", "navy")
    kwargs.setdefault("size", 10)
    kwargs.setdefault("alpha", 0.6)
    p.circle(x=x, y='HPO_term_name',
             source=source,
             **kwargs)

    p.add_tools(HoverTool(
        tooltips=[
            ( 'HPO_term_ID', '@HPO_term_ID' ),
            ( 'gene_num', '@gene_num' ),
            ( 'gene_ratio',   '@gene_ratio{%0.2f}' ),
            ( 'background_ratio',   '@background_ratio{%0.2f}' ),
            ( 'odd_ratio',   '@odd_ratio{0.00}' ),
            ( 'pvalue', '@pvalue' ),
            ( 'padj', '@padj' )
        ]
    ))

    p.xaxis.axis_label = " ".join([i.capitalize() for i in x.split("_")])
    p.yaxis.axis_label = "HPO Terms"

    return p


class LineagePlot(object):
    """Plot the tree structure of OBO file
    """
    def __init__(self, obo_file=None):
        if not obo_file:  # download HPO OBO file
            from hpoea.utils.download import get_hpo_obo
            obo_file = get_hpo_obo()
        from hpoea.utils.parse import parse_hpo_obo
        self.obo = parse_hpo_obo(obo_file)

    def plot(self, nodes, ax=None, **kwargs):
        import matplotlib.pyplot as plt

        G = self._get_subgraph(nodes)
        try:
            from networkx.drawing.nx_agraph import graphviz_layout
            pos = graphviz_layout(G)
        except ImportError as e:
            log.warning("pygraphviz is not installed using spring layout.")
            from networkx import spring_layout
            pos = spring_layout(G)

        if not ax:  # create an matplotlib ax
            fig, ax = plt.subplots()
        
        others = list(set(list(self._descendants)+list(self._ancestors)))
        nx.draw_networkx_nodes(G, pos=pos,
            nodelist=others,
            node_color='grey',
            node_size=500,
            alpha=0.6)
        nx.draw_networkx_nodes(G, pos=pos,
            nodelist=list(self._nodes),
            node_color='red',
            node_size=700,
            alpha=0.5)
        nx.draw_networkx_edges(G, pos=pos)
        nx.draw_networkx_labels(G, pos=pos)

        return ax

    def _get_subgraph(self, nodes):
        G = self.obo
        self._nodes = set(nodes)
        descendants = set()
        ancestors = set()
        for n in nodes:
            des_ = nx.descendants(G, n)
            descendants.update(des_)
            ans_ = nx.ancestors(G, n)
            ancestors.update(ans_)
        self._ancestors = ancestors - self._nodes
        self._descendants = descendants - self._nodes
        all_nodes = set()
        all_nodes.update(self._nodes)
        all_nodes.update(descendants)
        all_nodes.update(ancestors)
        subg = G.subgraph(all_nodes)
        return subg