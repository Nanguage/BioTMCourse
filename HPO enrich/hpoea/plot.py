from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, HoverTool


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
            ( 'gene_ratio',   '@gene_ratio{%0.2f}' ),
            ( 'background_ratio',   '@background_ratio{%0.2f}' ),
            ( 'odd_ratio',   '@odd_ratio{0.00}' ),
            ( 'pvalue', '@pvalue' ),
            ( 'padj', '@padj' )
        ]
    ))

    return p