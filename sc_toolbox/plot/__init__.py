from enum import Enum
from typing import Dict, List, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sb
from adjustText import adjust_text
from matplotlib import colors
from rich import print


class Colormaps(Enum):
    """Useful Colormaps for e.g. UMAPs."""

    grey_red = colors.LinearSegmentedColormap.from_list("grouping", ["lightgray", "red", "darkred"], N=128)
    grey_green = colors.LinearSegmentedColormap.from_list("grouping", ["lightgray", "limegreen", "forestgreen"], N=128)
    grey_yellow = colors.LinearSegmentedColormap.from_list("grouping", ["lightgray", "yellow", "gold"], N=128)
    grey_violet = colors.LinearSegmentedColormap.from_list(
        "grouping", ["lightgray", "mediumvioletred", "indigo"], N=128
    )
    grey_blue = colors.LinearSegmentedColormap.from_list("grouping", ["lightgray", "cornflowerblue", "darkblue"], N=128)


def custom_plot_size(width: int, height: int, dpi: int):
    """Create a custom axis object of desired sizes.

    Args:
        width: Desired plot width
        height: Desired plot height
        dpi: Desired plot DPI.

    Returns: Axis of desired sizes
    """
    fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)

    return fig.gca()


def standard_lineplot(
    data,
    order: List,
    xlabel: str,
    ylabel: str,
    hue=None,
    gene=None,
    smooth: bool = None,
    palette=None,
    title=None,
    rotation: int = None,
    figsize: Tuple[int, int] = (15, 5),
    tick_size=None,
    label_size=None,
    order_smooth: int = 3,
    confidence_interval=None,
    scatter=None,
    save: str = None,
):
    """Draws a standard line plot based on Seaborn's lmplot.

    Args:
        data: Data frame containing averaged expression values
        order: Order of x-axis labels from left to right
        xlabel: x-axis label
        ylabel: y-axis label
        hue: Subsets of the data which will be drawn on separate facets in the grid. Example: "condition"
        gene: Gene of interest
        smooth: Whether to smoothen (interpolate) the curve
        palette: Color palette. For example a list of colors.
        title: Title of the plot
        rotation: Rotation of the x-axis labels
        figsize: Size of the figure as specified in matplotlib
        tick_size: Size of the ticks as specified in matplotlib
        label_size: Size of the labels as specified in matplotlib
        order_smooth: If greater than 1, numpy.polyfit is used to estimate a polynomial regression
        confidence_interval: Confidence interval
        scatter: Set to true in order to add mean expression per sample in form of scatter point
        save: Path to save the plot to
    """
    if smooth:
        # Possible to set alpha of scatter with scatter_kws={'alpha': 0.1}
        if hue:
            cat = sb.lmplot(
                data=data,
                x=xlabel,
                y=gene,
                ci=confidence_interval,
                order=order_smooth,
                scatter=scatter,
                hue=hue,
                truncate=True,
                palette=palette,
            )
        else:
            cat = sb.lmplot(
                data=data,
                x=xlabel,
                y=gene,
                ci=confidence_interval,
                order=order_smooth,
                scatter=scatter,
                palette=palette,
            )

    else:
        # Removed Parameter order = order, as order should be given numerically anyways.
        if hue:
            cat = sb.catplot(data=data, x=xlabel, y=gene, linestyles="-", kind="point", hue=hue, palette=palette)
        else:
            cat = sb.catplot(data=data, x=xlabel, y=gene, linestyles="-", kind="point", palette=palette)
        if scatter:
            cat2 = sb.stripplot(data=data, x=xlabel, y=gene, palette=palette, hue=hue, size=7)
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
        full_save_name = f"{gene}_{save}"
        cat.fig.savefig(f"{full_save_name}", bbox_inches="tight")
        print(f"[bold blue]Saving figure to {full_save_name}")

    plt.show()
    plt.close()


def average_expression(
    gene_expression,
    genes,
    order: List[str],
    id_label: str = "identifier",
    xlabel: str = "days",
    cluster: str = "all",
    hue=None,
    palette: str = "tab:blue",
    figsize: Tuple[int, int] = (15, 6),
    smooth=None,
    rotation: int = None,
    order_smooth=None,
    conf_int=None,
    scatter=None,
    save: str = None,
):
    """Draw a line plot showing the gene expression over time. Expression values are averaged by individual sample.

    Args:
        gene_expression: Data frame containing gene expression values
        genes: List of genes for which individual line plots will be generated
        order: Order of x-axis labels from left to right
        id_label: Adata column in which sample id information is stored
        xlabel: x-axis label
        cluster: Which clusters to plot. Select 'all" if all clusters should be drawn.
        hue: Which value to color by
        figsize: Size of the figure as specified in matplotlib
        smooth: Set to true for smoothened line plot using polynomial regression
        rotation: set to True to rotate x-axis labels 90 degrees
        order_smooth: If greater than 1, use numpy.polyfit to estimate a polynomial regression
        conf_int: Size of the confidence interval for the regression estimate
        scatter: Set to True to add average expression values per sample ID as dots
        save: Path to save the plot to

    Example smooth:
        .. image:: /_images/average_expression_smooth.png

    Example raw:
        .. image:: /_images/average_expression_raw.png
    """
    for gene in genes:
        meanpid = gene_expression.groupby([id_label, xlabel])[gene].mean().reset_index()

        # cluster_label = ", ".join(cluster)
        cluster_label = ", ".join(cluster) if isinstance(cluster, list) else cluster
        standard_lineplot(
            meanpid,
            order=order,
            xlabel=xlabel,
            ylabel=f"Average expression in cluster {cluster_label}",
            hue=hue,
            gene=gene,
            smooth=smooth,
            palette=palette,
            title=gene,
            rotation=rotation,
            figsize=figsize,
            save=save,
            order_smooth=order_smooth,
            confidence_interval=conf_int,
            scatter=scatter,
        )


def average_expression_per_cluster(
    gene_expression,
    genes,
    order,
    obs=None,
    id_label: str = "identifier",
    xlabel: str = "days",
    cluster: str = "all",
    hue=None,
    figsize: Tuple[int, int] = (15, 6),
    smooth=None,
    rotation=None,
    tick_size: int = 12,
    label_size: int = 15,
    order_smooth=None,
    conf_int=None,
    palette=None,
    scatter=None,
    save: str = None,
):
    """Plots gene expression over time split by cluster identity.

    One line per cluster.

    Args:
        gene_expression: Data frame containing gene expression values
        genes: List of genes for which individual line plots will be generated
        order: Order of x-axis labels from left to right
        obs: Data frame containing meta data information
        xlabel: x-axis label
        cluster: Which clusters to plot. Select 'all" if all clusters should be drawn.
        id_label: Meta data column in which sample id information is stored
        hue: Split expression values by this grouping, one line per category will be drawn
        figsize: Size of the figure as specified in matplotlib
        smooth: Set to True for smoothened line plot using polynomial regression
        rotation: Set to True to rotate x-axis labels 90 degrees
        tick_size: Size of the ticks as specified in matplotlib
        label_size: Size of the labels as specified in matplotlib
        order_smooth: If greater than 1, use numpy.polyfit to estimate a polynomial regression
        conf_int: Size of the confidence interval for the regression estimate
        palette: Color palette that gets passed to Seaborn's lineplot. For example a list of colors.
        scatter: Set to True to add average expression values per sample ID as dots
        save: Path to save the plot to
    """
    for gene in genes:
        meanpid = gene_expression.groupby([id_label, xlabel])[gene].mean().reset_index()

        if hue:
            cell_types = {}
            combis = obs.groupby([id_label, hue]).groups.keys()

            for c in combis:
                cell_types[c[0]] = c[1]
            meanpid[hue] = [cell_types[label] for label in meanpid.identifier]

        # cluster_label = ", ".join(cluster)
        cluster_label = ", ".join(cluster) if isinstance(cluster, list) else cluster
        standard_lineplot(
            meanpid,
            order=order,
            xlabel=xlabel,
            ylabel=f"Average expression in cluster {cluster_label}",
            hue=hue,
            gene=gene,
            smooth=smooth,
            palette=palette,
            title=gene,
            tick_size=tick_size,
            label_size=label_size,
            rotation=rotation,
            figsize=figsize,
            save=save,
            order_smooth=order_smooth,
            confidence_interval=conf_int,
            scatter=scatter,
        )


def average_expression_split_cluster(
    gene_expression,
    genes,
    order,
    id_label="identifier",
    xlabel="days",
    hue="genotype",
    cluster=None,
    figsize=(15, 6),
    smooth=None,
    rotation=None,
    cols=None,
    tick_size=12,
    label_size=15,
    order_smooth=None,
    conf_int=None,
    scatter=None,
    save=None,
):
    """
    Plot average gene expression as line plots for multiple clusters at once.

    Args:
        gene_expression: Data frame containing gene expression values
        genes: List of genes for which individual line plots will be generated
        order: Order of x-axis labels from left to right
        id_label: Meta data column in which sample id information is stored
        xlabel: x-axis label
        hue: Split expression values by this grouping, one line per category, will be drawn
        cluster: Which clusters to plot. Select 'all" if all clusters should be drawn.
        figsize: Size of the figure as specified in matplotlib
        smooth: Set to True for smoothened line plot using polynomial regression
        rotation: x-axis label rotation
        cols: List of colors to use for line plot
        tick_size: Size of the ticks as specified in matplotlib
        label_size: Size of the labels as specified in matplotlib
        order_smooth: If greater than 1, numpy.polyfit is used to estimate a polynomial regression
        conf_int: Size of the confidence interval for the regression estimate
        scatter: Set to True to add average expression values per sample ID as dots
        save: Path to save the plot to

    Example smooth:
        .. image:: /_images/average_expression_per_cluster_smooth.png

    Example raw:
        .. image:: /_images/average_expression_per_cluster_raw.png
    """
    if cluster:
        if isinstance(cluster, list):
            ylab = f"Average expression in {', '.join(cluster)}"
        else:
            ylab = f"Average expression in {cluster}"
    else:
        ylab = "Average expression"

    for gene in genes:
        meanpid = gene_expression.groupby([id_label, hue, xlabel])[gene].mean().reset_index()

        standard_lineplot(
            meanpid,
            order=order,
            xlabel=xlabel,
            ylabel=ylab,
            hue=hue,
            gene=gene,
            smooth=smooth,
            palette=cols,
            title=gene,
            tick_size=tick_size,
            label_size=label_size,
            rotation=rotation,
            figsize=figsize,
            save=save,
            order_smooth=order_smooth,
            confidence_interval=conf_int,
            scatter=scatter,
        )


def average_expression_per_cell(
    gene_expression,
    genes,
    order,
    xlabel: str = "days",
    cluster: str = "all",
    hue=None,
    figsize: Tuple[int, int] = (15, 6),
    smooth=None,
    rotation=None,
    tick_size=12,
    label_size=15,
    order_smooth=None,
    conf_int=None,
    scatter=None,
    cols=None,
    save: str = None,
):
    """
    Plots the average gene expression as a line plot per cell.
    Ideally used when the scatter point should not be sample wise, but cell wise.
    Args:
        gene_expression: Data frame containing gene expression values
        genes: List of genes for which individual line plots will be generated
        order: Order of x-axis labels from left to right
        xlabel: x-axis label
        cluster: Which clusters to plot. Select 'all" if all clusters should be drawn.
        hue: Split expression values by this grouping, one line per category, will be drawn
        figsize: Size of the figure as specified in matplotlib
        smooth: Set to true for smoothened line plot using polynomial regression
        rotation: Set to True to rotate x-axis labels 90 degrees
        tick_size: Size of the ticks as specified in matplotlib
        label_size: Size of the labels as specified in matplotlib
        order_smooth: If greater than 1, use numpy.polyfit to estimate a polynomial regression
        conf_int: Size of the confidence interval for the regression estimate
        scatter: Set to True to add average expression values per sample ID as dots
        cols: List of colors to use for line plot
        save: Path to save the plot to

    """
    for gene in genes:
        cluster_label = ", ".join(cluster) if isinstance(cluster, list) else cluster
        standard_lineplot(
            gene_expression,
            order=order,
            xlabel=xlabel,
            ylabel=f"Average expression in cluster {cluster_label}",
            hue=hue,
            gene=gene,
            smooth=smooth,
            palette=cols,
            title=gene,
            tick_size=tick_size,
            label_size=label_size,
            rotation=rotation,
            figsize=figsize,
            save=save,
            order_smooth=order_smooth,
            confidence_interval=conf_int,
            scatter=scatter,
        )


def gene_expression_dpt_ordered(
    data,
    genes,
    xlabel,
    order=3,
    conf_int=95,
    figsize: Tuple[int, int] = (12, 6),
    condition=None,
    label_size: int = 15,
    cols=None,
    scale=None,
    ylim=None,
    save: str = None,
):
    """
    Plot smoothed expression of all cells ordered by pseudo time.

    Args:
        data: AnnData object
        genes: List of genes for which individual line plots will be generated
        xlabel: x-axis label
        order: Order of x-axis labels from left to right
        conf_int: Size of the confidence interval for the regression estimate
        figsize: Size of the figure as specified in matplotlib
        condition: Split expression values by this grouping, one line per category will be drawn
        label_size: Size of the labels as specified in matplotlib
        cols: List of colors to use for line plot
        scale: Set to True to scale expression value to a range between 0 and 1
        ylim: Upper limit on the y-axis if desired
        save: Path to save the plot to

    Example:
        .. image:: /_images/gene_expression_dpt_ordered.png

    Example with columns:
        .. image:: /_images/gene_expression_dpt_ordered_col.png
    """
    import matplotlib.patches as mpatches

    patches = []
    data = data.copy()
    fig, ax = plt.subplots(figsize=figsize)

    # use rainbow colour palette if no colours are specified
    if cols is None:
        from matplotlib import colors

        bins = len(np.unique(data.loc[:, condition])) if condition else len(genes)
        cmap = plt.cm.rainbow
        cmaplist = [cmap(i) for i in range(cmap.N)]
        cmap = colors.LinearSegmentedColormap.from_list("colours", cmaplist, N=bins)
        cols = [cmap(i) for i in range(bins)]

    # only working for one gene at a time for now
    if condition:
        conditions = np.unique(data.loc[:, condition])
        gene = genes[0]
        data = pd.pivot(data, columns=[condition])
        columns = [
            f"{data.columns.get_level_values(0)[i]}_{data.columns.get_level_values(1)[i]}"
            for i in range(len(data.columns.values))
        ]
        data.columns = columns
        data[xlabel] = data.filter(like=xlabel).sum(axis=1).values

        for i, con in enumerate(conditions):
            col = f"{gene}_{con}"

            if scale:
                data[col] = np.interp(data[col], (data[col].min(), data[col].max()), (0, +1))

            cat = sb.regplot(
                data=data, x=xlabel, y=col, scatter=False, order=order, truncate=True, ax=ax, color=cols[i], ci=conf_int
            )
            patches.append(mpatches.Patch(color=cols[i], label=col))

    else:
        for i, gene in enumerate(genes):
            if scale:
                data[gene] = np.interp(data[gene], (data[gene].min(), data[gene].max()), (0, +1))

            cat = sb.regplot(
                data=data,
                x=xlabel,
                y=gene,
                scatter=False,
                order=order,
                truncate=True,
                ax=ax,
                color=cols[i],
                ci=conf_int,
            )
            patches.append(mpatches.Patch(color=cols[i], label=gene))

    cat.set_ylabel("expression", size=label_size)
    cat.set_xlabel(xlabel, size=label_size)
    cat.tick_params(labelsize=label_size)
    sb.despine()

    plt.legend(handles=patches, loc="center left", bbox_to_anchor=(1.02, 0.5), prop={"size": label_size}, frameon=False)

    if ylim:
        cat.set(ylim=ylim)

    if save:
        plt.savefig(f"{save}", bbox_to_anchor="tight")
        print("[bold blue]Saving figure to {save}")

    plt.show()
    plt.close()


def relative_frequencies_boxplots(
    relative_frequencies: pd.DataFrame,
    cluster,
    cols,
    order,
    xlabel: str = "days",
    hue: str = "batch",
    figsize: Tuple[int, int] = (15, 6),
    width: float = 0.5,
    jitter=None,
    save=None,
) -> None:
    """Plots the relative frequencies as split boxplots.

    Use calc_relative_frequencies to get the required input format.

    Args:
        relative_frequencies: Calculated by calc_relative_frequencies as Pandas DataFrame
        cluster: Cluster to be plotted
        cols: List of colors to use for boxes
        order: Order of x-axis labels from left to right
        xlabel: x-axis label
        hue: Value to color by
        figsize: Size of the figure as specified in matplotlib
        width: Width of the plot as specified in matplotlib
        jitter: Set to True for individual dots per sample
        save: Path to save the plot to

    Example:
        .. image:: /_images/relative_frequencies_boxplots.png
    """
    # Subset according to order
    relative_frequencies = relative_frequencies.loc[relative_frequencies[xlabel].isin(order)]

    split_boxplot(
        relative_frequencies,
        order=order,
        xlabel=xlabel,
        ylabel="relative frequency",
        hue=hue,
        column=cluster,
        cols=cols,
        width=width,
        title=cluster,
        figsize=figsize,
        jitter=jitter,
        save=save,
    )


def split_boxplot(
    table,
    order,
    xlabel: str,
    ylabel: str,
    column=None,
    hue=None,
    cols=None,
    width: float = 1,
    title=None,
    figsize: Tuple[int, int] = (15, 6),
    jitter=None,
    save: str = None,
) -> None:
    """Draws a boxsplit split by hue.

    Args:
        table: Table containing the data to draw the boxplots for
        order: Order of the boxplot labels
        xlabel: x-axis label
        ylabel: y-axis label
        column:
        hue: Value to split relative frequencies by
        cols: List of colors to use for boxes
        width: Width of the desired plot
        title: Title of the plot
        figsize: Size of the figure as specified in matplotlib
        jitter: Set to True for individual dots per sample
        save: Path to save the plot to
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


def marker_dendrogram(
    marker_table: pd.DataFrame,
    threshold: float = 0.7,
    column: str = "cluster",
    log_fc_key: str = "log_FC",
    label_size: int = 10,
    orientation: str = "top",
    figsize: Tuple[int, int] = (10, 4),
    save: str = None,
):
    """Plots a dendogram of used marker genes.

    Args:
        marker_table: A marker table as generated by sct.calc.extended_marker_table
        threshold: Threshold for the log fold change
        column: Column to create pivot by; usually just the clusters
        log_fc_key: Key for the stored log fold changes in the marker table
        label_size: Font size of the labels
        orientation: Orientation of the figure; Currently just 'top' or no orientation
        figsize: Size of the figure as specified in matplotlib
        save: Path to save the plot to

    Example:
        .. image:: /_images/marker_dendrogram.png
    """
    import scipy.cluster.hierarchy as hc

    marker_table = marker_table[marker_table[log_fc_key] > threshold]

    marker_table = marker_table.pivot(index="gene", columns=column, values=log_fc_key)
    marker_table.fillna(value=0, inplace=True)

    corr = 1 - marker_table.corr()
    corr = hc.distance.squareform(corr)  # convert to condensed
    z = hc.linkage(corr, method="complete")
    plt.figure(figsize=figsize)
    rot = 90 if orientation == "top" else 0
    hc.dendrogram(
        z,
        labels=marker_table.columns,
        leaf_rotation=rot,
        color_threshold=0,
        orientation=orientation,
        leaf_font_size=label_size,
        above_threshold_color="black",
    )
    plt.yticks(size=label_size)
    if save is None:
        plt.show()
    else:
        plt.savefig("{save}")
        print(f"[bold blue]Saving figure to {save}")
    plt.close()


def volcano_plot(
    table,
    fdr_thresh: float = None,
    log_fc_thresh: float = 0,
    adj_p_val: str = "adj_p_val",
    log_fc: str = "avg_logFC",
    gene: str = "gene",
    sig_col: str = "tab:orange",
    col: str = "tab:blue",
    figsize: Tuple[int, int] = (8, 6),
    save=None,
):
    """
    Scatter plot of differential gene expression results generated by diffxpy

    Args:
        table: diffxpy generated table of results
        fdr_thresh: -log(FDR) threshold for labeling genes. If set to None, we
            will consider the 99th percentile of -log(FDR) values the threshold.
        log_fc_thresh: absolute(log_fc) threshold for labeling genes.
        adj_p_val: Label of the adjusted p value, these are considered FDRs
        log_fc: Label of the log fold change
        gene: Label of column with gene names
        col: Color of dots
        sig_col: Colour of dots surpassing defined FDR threshold
        figsize: Size of the figure as specified in matplotlib
        save: Path to save the plot to

    Example:
        .. image:: /_images/diffxpy_volcano.png
    """
    table["-log_FDR"] = -np.log(table[adj_p_val])

    # take the 99% quantile by default for highlighting
    if not fdr_thresh:
        fdr_thresh = np.percentile(table.loc[:, "-log_FDR"], 99)
    if not log_fc_thresh:
        log_fc_thresh = 0
    lowqval_highfc_de = table.loc[(table["-log_FDR"] > fdr_thresh) & (abs(table[log_fc]) >= log_fc_thresh)]
    other_de = table.loc[~table.index.isin(lowqval_highfc_de.index)]

    fig, ax = plt.subplots()
    fig.set_size_inches(figsize)

    sb.regplot(x=other_de[log_fc], y=other_de["-log_FDR"], fit_reg=False, scatter_kws={"s": 6})
    sb.regplot(x=lowqval_highfc_de[log_fc], y=lowqval_highfc_de["-log_FDR"], fit_reg=False, scatter_kws={"s": 6})
    ax.set_xlabel("log2 FC", fontsize=20)
    ax.set_ylabel("-log Q-value", fontsize=20)
    ax.tick_params(labelsize=15)
    ax.grid(False)

    # Label names and positions
    x = [i - 0.1 for i in lowqval_highfc_de[log_fc]]
    y = [i + 0.1 for i in lowqval_highfc_de["-log_FDR"]]
    labels = lowqval_highfc_de[gene]
    max_n_labels = 50
    if len(labels) > max_n_labels:
        print(f"[bold yellow]Warning: given your thresholds, more than {max_n_labels} genes would have to be labeled.")
        print(
            "[bold yellow]To prevent overcrowding of your plot, make your thresholds stricter.\n"
            "We will leave out the labels for now."
        )
    else:
        # plot labels, and use adjust_text to make sure that labels don't overlap:
        labels = [
            plt.text(x, y, label, ha="center", va="center") for x, y, label in zip(x, y, labels) if not x == np.inf
        ]
        adjust_text(labels)
    if save:
        fig.savefig(f"{save}")
    else:
        plt.show()
    plt.close()


def cluster_composition_stacked_barplot(
    relative_frequencies: pd.DataFrame,
    xlabel: str = "name",
    figsize: Tuple[int, int] = (6, 10),
    width: float = 0.8,
    order=None,
    error_bar=None,
    label_size: int = 15,
    tick_size: int = 13,
    capsize: int = None,
    margins: Tuple[float, float] = (0.02, 0.04),
    colors=None,
    save: str = None,
):
    """Plot relative frequencies as a stacked barplot.

    Args:
        relative_frequencies: Data frame containing relative Frequencies as calculated by calc_relFreq()
        xlabel: x-axis label
        figsize: Size of the figure as specified in matplotlib
        width: Width of the bars
        order: Order of x-axis labels from left to right
        error_bar: Set to True to add error bars (only possible when grouping the frequencies)
        tick_size: Size of the ticks as specified in matplotlib
        label_size: Size of the labels as specified in matplotlib
        capsize: Size of the horizontal lines of the error bar
        margins: Change margins of the plot if desired
        colors: List of colors to use for the bands
        save: Path to save the plot to

    Example:
        .. image:: /_images/cluster_composition_stacked_barplot.png
    """
    import matplotlib.patches as mpatches

    if not colors:
        raise ValueError("Colors was not passed. Obtain them from e.g. adata.uns['cluster_key_colors']")

    patches = []
    fig, ax = plt.subplots()
    fig.set_size_inches(figsize)
    order = np.unique(relative_frequencies.loc[:, xlabel]) if order is None else order
    ci = 95 if error_bar else None
    ax.margins(margins[0], margins[1])
    cell_types = np.flip([col for col in relative_frequencies.columns if col not in ["identifier", xlabel]])
    # cell_types = np.flip(np.setdiff1d(relFreqs.columns, ["identifier", xlabel]))

    bars = pd.DataFrame(index=order, data=np.zeros(len(order)))
    plot_data = pd.DataFrame(relative_frequencies.loc[:, xlabel])

    for i, typ in enumerate(cell_types):
        sum_up = [
            relative_frequencies.loc[:, typ].values[i] + bars.loc[g].values[0]
            for i, g in enumerate(relative_frequencies.loc[:, xlabel])
        ]
        plot_data[typ] = sum_up
        bars.iloc[:, 0] = (
            bars.iloc[:, 0] + relative_frequencies.loc[:, [typ, xlabel]].groupby(xlabel).mean().loc[order, typ]
        )

    for i, typ in enumerate(reversed(cell_types)):
        fig = sb.barplot(
            data=plot_data, x=xlabel, y=typ, order=order, ci=ci, errcolor="black", color=colors[i], capsize=capsize
        )
        patches.append(mpatches.Patch(color=colors[i], label=typ))

    ax.set_xlabel(xlabel, size=label_size)
    ax.set_ylabel("relative frequency", size=label_size)
    ax.tick_params(labelsize=tick_size)
    ax.set_xticklabels(labels=order, rotation="vertical")

    # Change the bar width
    for bar in fig.patches:
        centre = bar.get_x() + bar.get_width() / 2.0
        bar.set_x(centre - width / 2.0)
        bar.set_width(width)

    plt.legend(handles=patches, loc="center left", bbox_to_anchor=(1.02, 0.5), prop={"size": tick_size}, frameon=False)
    if save:
        plt.savefig(f"{save}")
        print(f"[bold blue]Saving Figure to {save}")
    plt.show()
    plt.close()


def gene_boxplot(
    table,
    palette: List[str],
    xlabel: str = "cell_types",
    hue: str = None,
    figsize: Tuple[int, int] = (10, 5),
    legend=True,
    score="Axin2",
    scatter=None,
    rotate=False,
    width=0.7,
    save=None,
):
    """Plot gene values as split boxplots.

    Args:
        table: Pandas DataFrame
        palette:
        xlabel: x-axis label
        hue:
        figsize: Size of the figure as specified in matplotlib
        legend: Whether to draw a legend or not
        score:
        scatter:
        rotate:
        width: Width of the desired plot
        save: Path to save the plot to

    Example:
        .. image:: /_images/gene_boxplot.png
    """
    sb.set_style("ticks")
    fig, ax = plt.subplots()
    fig.set_size_inches(figsize)

    sf = False if scatter else True
    if hue:
        fig = sb.boxplot(data=table, x=xlabel, y=score, width=width, hue=hue, showfliers=sf, palette=palette)
        if scatter:
            fig = sb.stripplot(data=table, x=xlabel, y=score, palette=["black"], size=4, hue=hue, dodge=True)
    else:
        fig = sb.boxplot(data=table, x=xlabel, y=score, width=width, showfliers=sf, palette=palette)
        if scatter:
            fig = sb.stripplot(data=table, x=xlabel, y=score, palette=["black"], size=4, dodge=True)

    if rotate:
        fig.set_xticklabels(fig.get_xticklabels(), rotation=90)
    else:
        fig.set_xticklabels(fig.get_xticklabels())

    if legend:
        ax.legend(bbox_to_anchor=(1.05, 1.06))
    else:
        ax.legend_.remove()

    plt.setp(ax.artists, edgecolor="black")
    plt.setp(ax.lines, color="black")
    sb.despine()  # to not show ouline box

    if save:
        print(f"Saving to {save}")
        plt.savefig(save, bbox_to_anchor="tight")
    plt.show()


def colors_overview(colors: Dict, ncols: int = 2, figsize: Tuple[int, int] = (8, 5), save: str = None):
    """Draw an overview plot of all used colors.

    Args:
        colors: Dictionary of color name and color
        ncols: How many columns for the plot
        figsize: Size of the figure as specified in matplotlib
        save: Path to save the plot to

    Example:
        .. image:: /_images/colors.png
    """
    from matplotlib import colors as mcolors

    # Sort colors by hue, saturation, value and name.
    by_hsv = sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgba(color)[:3])), name) for name, color in colors.items())
    sorted_names = [name for hsv, name in by_hsv]

    n = len(sorted_names)
    nrows = n // ncols + 1

    fig, ax = plt.subplots(figsize=figsize)

    # Get height and width
    x, y = fig.get_dpi() * fig.get_size_inches()
    h = y / (nrows + 1)
    w = x / ncols

    for i, name in enumerate(sorted_names):
        col = i % ncols
        row = i // ncols
        y = y - (row * h) - h

        xi_line = w * (col + 0.05)
        xf_line = w * (col + 0.25)
        xi_text = w * (col + 0.3)

        ax.text(
            xi_text,
            y,
            "%s %s" % (name, colors[name]),
            fontsize=(h * 0.4),
            horizontalalignment="left",
            verticalalignment="center",
        )
        ax.hlines(y + h * 0.1, xi_line, xf_line, color=colors[name], linewidth=(h * 0.6))

    ax.set_xlim(0, x)
    ax.set_ylim(0, y)
    ax.set_axis_off()

    fig.subplots_adjust(left=0, right=1, top=1, bottom=0, hspace=0, wspace=0)
    if save:
        print(f"Saving to {save}")
        plt.savefig(save, bbox_to_anchor="tight")
    plt.show()


def relative_frequencies_lineplot(
    relative_frequencies: pd.DataFrame,
    order,
    cluster,
    xlabel: str = "days",
    ylabel: str = "relative frequency",
    hue: str = None,
    smooth: bool = None,
    cols=None,
    title: str = None,
    rotation: int = None,
    figsize: Tuple[int, int] = (15, 5),
    tick_size: int = None,
    label_size: int = None,
    order_smooth: int = 3,
    conf_int=None,
    scatter=None,
    save: str = None,
):
    """Plot relative frequencies as a line plot.

    Args:
        relative_frequencies: Data frame containing relative Frequencies as calculated by calc_relFreq()
        order: Order of x-axis labels from left to right
        cluster: Which cluster to plot
        xlabel: x-axis label
        ylabel: y-axis label
        hue: Value to color by
        smooth: Whether to smoothen the plot
        cols: List of colors to use for line plot
        title: Title of the plot
        rotation: Rotation of the x-axis labels
        figsize: Size of the figure as specified in matplotlib
        tick_size: Size of the ticks as specified in matplotlib
        label_size: Size of the labels as specified in matplotlib
        order_smooth: If greater than 1, numpy.polyfit is used to estimate a polynomial regression
        conf_int: Size of the confidence interval for the regression estimate
        scatter: Set to True to add average expression values per sample ID as dots
        save: Path to save the plot to

    Example:
        .. image:: /_images/relative_frequencies_lineplots.png
    """
    if hue:
        sub_freqs = relative_frequencies.loc[:, [cluster] + [xlabel, hue]]
        sub_freqs = pd.melt(sub_freqs, id_vars=[xlabel, hue])
    else:
        sub_freqs = relative_frequencies.loc[:, [cluster] + [xlabel]]
        sub_freqs = pd.melt(sub_freqs, id_vars=[xlabel])

    standard_lineplot(
        sub_freqs,
        order=order,
        xlabel=xlabel,
        ylabel=ylabel,
        hue=hue,
        gene="value",
        smooth=smooth,
        palette=cols,
        title=title,
        rotation=rotation,
        figsize=figsize,
        tick_size=tick_size,
        label_size=label_size,
        order_smooth=order_smooth,
        confidence_interval=conf_int,
        scatter=scatter,
        save=save,
    )


def annotated_cell_type_umap(
    adata,
    primary_color: Union[str, Sequence[str]],
    cell_type_color: str,
    legend_loc: str = "on data",
    legend_fontsize: int = 8,
    title: str = "Plot title",
    palette=None,
    cmap=None,
    figsize=(8, 6),
    save=None,
):
    """Plots a UMAP which is colored by the primary_color, but also draws all labels on top of all clusters.

    Args:
        adata: AnnData object
        primary_color: Primary color to color all cells by, e.g. 'genotype'
        cell_type_color: Key containing all cell types, e.g. 'cell_type'
        legend_loc: Location of the legend (default: 'on data')
        legend_fontsize: Font size of the legend (default: 8)
        title: Title of the plot
        palette: Color
        cmap: Color map of the UMAP
        figsize: Size of the figure
        save: Path to save the plot to

    Returns:
        fig and axs Matplotlib objects

    Example:
        .. image:: /_images/annotated_cell_type_umap.png
    """
    fig, axs = plt.subplots(figsize=figsize)
    sc.pl.umap(adata, color=primary_color, show=False, palette=palette, cmap=cmap, ax=axs)
    sc.pl.umap(
        adata,
        color=cell_type_color,
        alpha=0,
        legend_loc=legend_loc,
        legend_fontsize=legend_fontsize,
        title=title,
        show=False,
        ax=axs,
    )

    if save:
        fig.savefig(save, dpi=1200, format="pdf", bbox_inches="tight")

    return fig, axs


def genotype_vs_genotype_umaps(
    adata,
    genotype_key: str,
    genotype_label_1: str,
    genotype_label_2: str,
    color: str,
    hide_one_legend: bool = True,
    figsize: Tuple[int, int] = (12, 6),
):
    """Plots a two UMAPs of genotypes next to each other displaying only the colors of the second UMAP.

    Args:
        adata: AnnData object
        genotype_key: Key of the genotypes
        genotype_label_1: Name of the first genotype; Must be contained in the genotypes
        genotype_label_2: Name of the second genotype; Must be contained in the genotypes
        color: Key to color by
        hide_one_legend: Whether to hide the legend of the genotype_label_1
        figsize: Size of the figure

    Example:
        .. image:: /_images/genotype_vs_genotype_umaps.png
    """
    genotype_data_1 = adata[adata.obs[genotype_key].isin([genotype_label_1])].copy()
    genotype_data_2 = adata[adata.obs[genotype_key].isin([genotype_label_2])].copy()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    sc.pl.umap(
        genotype_data_1,
        color=color,
        ax=ax1,
        palette=sc.pl.palettes.default_20,
        legend_fontsize="xx-small",
        size=40,
        show=False,
    )
    if hide_one_legend:
        ax1.get_legend().remove()
    ax1.set_title(genotype_label_1)
    sc.pl.umap(
        genotype_data_2,
        color=color,
        ax=ax2,
        palette=sc.pl.palettes.default_20,
        legend_fontsize="xx-small",
        size=40,
        show=False,
    )
    _ = ax2.set_title(genotype_label_2)
