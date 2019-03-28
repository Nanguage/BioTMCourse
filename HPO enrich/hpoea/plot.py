from bokeh.plotting import figure

def bar_plot(enrich_table, figsize, show_items=10):
    """ Draw the enrichment bar chart.

    Args
    ----
    enrich_table : pandas.DataFrame
        Enrichment result table, see `hpoea.enrich.GSEA`.
    figsize : tuple, optional
        The size of output figure, (width, height).
    show_items : int, optional
        Number of items to show. default 10.
    """
    if figsize:
        width, height = figsize
    else:  # set default figsize
        grid_height = 70
        height = grid_height * show_items
        width = 500
    p = figure(plot_width=width, plot_height=height)

    df = enrich_table.iloc[:show_items, :]
    p.circle()
    return p


def dot_plot(enrich_table, figsize=None, show_items=10):
    """ Draw the enrichment dot chart.

    Args
    ----
    enrich_table : pandas.DataFrame
        Enrichment result table, see `hpoea.enrich.GSEA`.
    figsize : tuple, optional
        The size of output figure, (width, height).
    show_items : int, optional
        Number of items to show. default 10.
    """
    if figsize:
        width, height = figsize
    else:  # set default figsize
        grid_height = 70
        height = grid_height * show_items
        width = 500
    p = figure(plot_width=width, plot_height=height)

    df = enrich_table.iloc[:show_items, :]
    p.circle()
    return p