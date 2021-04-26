from typing import List, Tuple

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
    hue=None,
    gene=None,
    smooth: bool = None,
    cols=None,
    title=None,
    rotation: int = None,
    figsize: Tuple[int, int] = (15, 5),
    tick_size=None,
    label_size=None,
    order_smooth: int = 3,
    conf_int=None,
    scatter=None,
    save: str = None,
):
    """
    Draws a standard line plot based on Seaborn's lmplot.

    Args:
        data: Data to plot. Usually AnnData object
        order: Order of x-axis labels from left to right
        xlabel: x-axis label
        ylabel: y-axis label
        hue: Subsets of the data which will be drawn on separate facets in the grid. Example: "condition"
        gene:
        smooth: Whether to smoothen (interpolate) the curve
        cols:
        title: Title of the plot
        rotation: Rotation of the x-axis labels
        figsize: Size of the figure as specified in matplotlib.
        tick_size: Size of the ticks as specified in matplotlib.
        label_size: Size of the labels as specified in matplotlib.
        order_smooth: If greater than 1, numpy.polyfit is used to estimate a polynomial regression.
        conf_int: Confidence interval
        scatter:
        save: Path to save the plot to
    """
    if smooth:
        # Possible to set alpha of scatter with scatter_kws={'alpha': 0.1}
        if hue:
            cat = sb.lmplot(
                data=data,
                x=xlabel,
                y=gene,
                ci=conf_int,
                order=order_smooth,
                scatter=scatter,
                hue=hue,
                truncate=True,
                palette=cols,
            )
        else:
            cat = sb.lmplot(data=data, x=xlabel, y=gene, ci=conf_int, order=order_smooth, scatter=scatter, palette=cols)

    else:
        # Removed Parameter order = order, as order should be given numerically anyways.
        if hue:
            cat = sb.catplot(data=data, x=xlabel, y=gene, linestyles="-", kind="point", hue=hue, palette=cols)
        else:
            cat = sb.catplot(data=data, x=xlabel, y=gene, linestyles="-", kind="point", palette=cols)
        if scatter:
            cat2 = sb.stripplot(data=data, x=xlabel, y=gene, palette=cols, hue=hue, size=7)
            if hue:
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


def split_boxplot(
    table,
    order,
    xlabel,
    ylabel,
    column=None,
    hue=None,
    cols=None,
    width=1,
    title=None,
    figsize=(15, 6),
    jitter=None,
    save=None,
):
    """
    Draws a boxsplit split by

    Args:
        table:
        order:
        xlabel:
        ylabel:
        column:
        hue:
        cols:
        width:
        title:
        figsize:
        jitter:
        save:
    """
    fig, ax = plt.subplots()
    fig.set_size_inches(figsize)

    if cols is not None:
        fig = sb.boxplot(data=table, hue=hue, x=xlabel, y=column, order=order, width=width, palette=cols)
    else:
        fig = sb.boxplot(data=table, hue=hue, x=xlabel, y=column, order=order, width=width)

    if jitter is not None:
        fig = sb.swarmplot(data=table, color="black", x=xlabel, y=column, order=order)

    if hue is not None:
        plt.legend(loc="upper right")

    if title:
        fig.set_title(title, size=15)

    fig.set_xlabel(xlabel, size=15)
    fig.set_ylabel(ylabel, size=15)
    fig.tick_params(labelsize=12)

    if save:
        fig.get_figure().savefig("{save}")

    plt.show()
    plt.close()
