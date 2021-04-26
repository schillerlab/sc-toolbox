from typing import List

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sb


def custom_plot_size(width: int, height: int, dpi: int):
    """
    Create a custom axis object of desired sizes.

    Args:
        width: Desired plot width
        height: Desired plot height
        dpi: Desired plot DPI.

    Returns: Axis of desired sizes
    """
    fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)

    return fig.gca()


def standart_lineplot(
    data,
    order: List,
    xlabel: str,
    ylabel: str,
    typ=None,
    gene=None,
    smooth=None,
    cols=None,
    title=None,
    rotation=None,
    figsize=(15, 5),
    tick_size=None,
    label_size=None,
    order_smooth=3,
    conf_int=None,
    scatter=None,
    save=None,
):
    """

    Args:
        data: Data to plot. Usually AnnData object
        order: Order of x-axis labels from left to right
        xlabel: 
        ylabel:
        typ:
        gene:
        smooth:
        cols:
        title:
        rotation:
        figsize:
        tick_size:
        label_size:
        order_smooth:
        conf_int:
        scatter:
        save:

    Returns:

    """
    if smooth:
        # Possible to set alpha of scatter with scatter_kws={'alpha': 0.1}
        if typ:
            cat = sb.lmplot(
                data=data,
                x=xlabel,
                y=gene,
                ci=conf_int,
                order=order_smooth,
                scatter=scatter,
                hue=typ,
                truncate=True,
                palette=cols,
            )
        else:
            cat = sb.lmplot(data=data, x=xlabel, y=gene, ci=conf_int, order=order_smooth, scatter=scatter, palette=cols)

    else:
        # Removed Parameter order = order, as order should be given numerically anyways.
        if typ:
            cat = sb.catplot(data=data, x=xlabel, y=gene, linestyles="-", kind="point", hue=typ, palette=cols)
        else:
            cat = sb.catplot(data=data, x=xlabel, y=gene, linestyles="-", kind="point", palette=cols)
        if scatter:
            cat2 = sb.stripplot(data=data, x=xlabel, y=gene, palette=cols, hue=typ, size=7)
            if typ:
                cat2.legend_.remove()

    cat.set(xticks=np.unique(data.loc[:, xlabel]))

    cat.set_xticklabels(order)
    cat.fig.set_size_inches(figsize)

    if rotation:
        cat.ax.set_xticklabels(order, rotation="vertical")
    cat.ax.set_title(title, size=label_size)
    cat.ax.set_xlabel(xlabel, size=label_size)
    cat.ax.set_ylabel(ylabel, size=label_size)
    cat.ax.tick_params(labelsize=tick_size)

    if save:
        cat.fig.savefig(f"{save}", bbox_inches="tight")
        print(f"[bold blue]Saving figure to {save}")
    plt.show()
    plt.close()
