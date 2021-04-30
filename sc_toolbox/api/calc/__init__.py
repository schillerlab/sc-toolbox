import os
from typing import List
from typing import Literal
from typing import Optional

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from rich import print

WORKING_DIRECTORY = os.path.dirname(__file__)


def generate_expression_table(
    adata,
    cluster: str = "all",
    subset_by: str = "cell_type",
    xlabel: str = "days",
    hue: str = None,
    use_raw: bool = None,
):
    """

    Args:
        adata: Anndata object
        cluster: Which label of the subsets to generate the table for. Use 'all' if for all subsets.
        subset_by: Which label to subset the clusters by
        xlabel: x-axis
        hue: Value to color by
        use_raw: Whether to use adata.raw.X for the calculations

    Returns:
         Gene expression table
    """
    if cluster == "all":
        cells = adata.obs_names
    else:
        cells = [True if val in cluster else False for val in adata.obs[subset_by]]
    if use_raw:
        gen_expression_table = pd.DataFrame(
            adata[cells].raw.X.todense(), index=adata[cells].obs_names, columns=adata[cells].raw.var_names
        )
    else:
        gen_expression_table = pd.DataFrame(
            adata[cells].X, index=adata[cells].obs_names, columns=adata[cells].var_names
        )

    gen_expression_table["identifier"] = adata[cells].obs["identifier"]
    gen_expression_table[xlabel] = adata[cells].obs[xlabel]
    if hue:
        # For multiple cluster, split internally per condition
        if isinstance(cluster, list) and len(cluster) > 1 and subset_by != hue:
            gen_expression_table[hue] = [f"{t}_{c}" for t, c in zip(adata[cells].obs[hue], adata[cells].obs[subset_by])]

        else:
            gen_expression_table[hue] = adata[cells].obs[hue]

    return gen_expression_table


def relative_frequencies(adata, group_by: str = "cell_type", xlabel: str = "days", condition: str = "batch"):
    """
    Calculates the relative frequencies of conditions grouped by an observation.

    Args:
        adata: AnnData Objet containing the data
        group_by:
        xlabel: x-axis label
        condition:

    Returns:
        Relative frequencies in a Pandas DataFrame

    """
    freqs = adata.obs.groupby(["identifier", group_by]).size()
    samples = np.unique(adata.obs["identifier"])
    ind = adata.obs[group_by].cat.categories

    relative_frequencies = [freqs[ident] / sum(freqs[ident]) for ident in samples]
    relative_frequencies = pd.DataFrame(relative_frequencies, columns=ind, index=samples).fillna(0)

    # relFreqs[xlabel] = grouping.loc[samples, xlabel]  ## when using Grouping Table
    cell_types = {}
    combis = adata.obs.groupby(["identifier", xlabel]).groups.keys()

    for c in combis:
        cell_types[c[0]] = c[1]
    relative_frequencies[xlabel] = [cell_types[label] for label in relative_frequencies.index]  # type: ignore

    # Todo, add for condition
    if condition:
        combis = adata.obs.groupby(["identifier", condition]).groups.keys()
        for c in combis:
            cell_types[c[0]] = c[1]
        relative_frequencies[condition] = [cell_types[label] for label in relative_frequencies.index]  # type: ignore

    return relative_frequencies


def relative_frequency_per_cluster(adata, group_by: str = "cell_type", xlabel: str = "days", condition=None):
    """
    Calculates relative frequencies per cluster

    Args:
        adata: AnnData object containing the data
        group_by: The label to group by for the clusters
        xlabel: x-axis label
        condition: condition to combine by

    Returns:
        Pandas DataFrame of relative frequencies
    """
    frequencies = adata.obs.groupby([group_by, xlabel]).size()
    celltypes = np.unique(adata.obs[group_by])
    ind = adata.obs[xlabel].cat.categories

    relative_frequencies = [frequencies[ident] / sum(frequencies[ident]) for ident in celltypes]
    relative_frequencies = pd.DataFrame(relative_frequencies, columns=ind, index=celltypes).fillna(0)

    cell_types = {}
    combinations = adata.obs.groupby([group_by, xlabel]).groups.keys()

    for combination in combinations:
        cell_types[combination[0]] = combination[1]
    relative_frequencies[group_by] = relative_frequencies.index  # type: ignore

    # Todo, add for condition
    if condition:
        combinations = adata.obs.groupby([group_by, condition]).groups.keys()
        for combination in combinations:
            cell_types[combination[0]] = combination[1]
        relative_frequencies[condition] = [cell_types[label] for label in relative_frequencies.index]  # type: ignore

    return relative_frequencies


def correlate_to_signature(
    adata,
    marker: pd.DataFrame,
    log_fc_threshold: float = 0.7,
    cell_type: str = "AT2 cells",
    cell_type_label: str = "cell_type",
    log_fc_label: str = "logfoldchange",
    gene_label: str = "gene",
    use_raw: bool = True,
):
    """
    Correlations Score (based on cell type signature (logFC)) - alternative to sc.tl.score

    Args:
        adata: AnnData object containing the data
        marker: Pandas DataFrame containing marker genes
        log_fc_threshold: Log fold change label
        cell_type: Cell type to calculate the correlation for
        cell_type_label: Label of all cell types in the AnnData object
        log_fc_label: Label of fold change in the AnnData object
        gene_label: Label of genes in the AnnData object
        use_raw: Whether to use adata.raw.X

    Returns:
        List of correlations
    """
    from scipy.sparse import issparse

    topmarker = marker[marker.loc[:, cell_type_label] == cell_type]
    topmarker = topmarker.loc[topmarker.loc[:, log_fc_label] > log_fc_threshold, [gene_label, log_fc_label]]

    gene_names = list(np.intersect1d(adata.var_names, topmarker.loc[:, gene_label].astype(str)))
    topmarker = topmarker[topmarker.loc[:, gene_label].isin(gene_names)]
    print(f"[bold blue]{len(gene_names)} genes used for correlation score to {cell_type}")

    if use_raw:
        if issparse(adata.raw.X):
            gene_expression = adata.raw[:, gene_names].X.todense()
        else:
            gene_expression = adata.raw[:, gene_names].X
    else:
        if issparse(adata.X):
            gene_expression = adata[:, gene_names].X.todense()
        else:
            gene_expression = adata[:, gene_names].X
    gene_expression = pd.DataFrame(gene_expression.T, index=gene_names)

    # For each cell separately
    gene_expression = pd.DataFrame.fillna(gene_expression, value=0)
    res = [
        np.correlate(topmarker.loc[:, log_fc_label], gene_expression.iloc[:, c])[0]
        for c in range(gene_expression.shape[1])
    ]

    return res


def remove_outliers(cords, eps: int = 1, min_samples: int = 2):
    """
    Remove outlying cells based on UMAP embeddings with DBScan (density based clustering)
    Call as: sub.obs["d_cluster"] = remove_outliers(sub.obsm["X_umap"], min_samples = 10)

    Args:
        cords: adata UMAP coordinates, typically adata.obsm["X_umap"]
        eps: Maximum distance between two clusters to still be considered neighbors
        min_samples: Minimum samples of a cluster

    Returns:
        Pandas DataFrame of clusters

    """
    from sklearn.cluster import DBSCAN
    from natsort import natsorted

    clustering = DBSCAN(eps=eps, min_samples=min_samples).fit(cords)
    cluster = clustering.labels_.astype("U")

    return pd.Categorical(cluster, categories=natsorted(np.unique(cluster)))


def add_percentages(adata, table, ids, group_by: str, threshold: int = 0, gene_label: str = "gene"):
    """
    Add columns to existing diffxpy table specifying percentage of expressing cells

    Args:
        adata: AnnData object containing the data
        table: Table as generated by diffxpy
        ids:
        group_by: Label to group by
        threshold:
        gene_label: Label of the genes

    Returns:
        Table containing percentage of expressing cells
    """
    for ident in ids:
        cells = adata.obs_names[adata.obs[group_by] == ident]
        data_temp = pd.DataFrame(
            ((adata[cells].layers["counts"] > threshold).sum(0) / adata[cells].layers["counts"].shape[0]).T,
            index=adata.var_names,
        )

        if gene_label == "index":
            table["pct.%s" % ident] = data_temp.loc[table.index.values].values
        else:
            table["pct.%s" % ident] = data_temp.loc[table.loc[:, gene_label]].values
    return table


def ranksums_between_groups(
    table, id1: str = "bystander", id2: str = "infected", xlabel: str = "condition", cells=None, score: str = "Axin2"
):
    """
    Perform Wilcoxon Rank-sum test between two groups.

    Args:
        table:
        id1:
        id2:
        xlabel: x-axis label
        cells:
        score:

    Returns:
        Pandas DataFrame containing test statistic and p-value
    """
    from scipy import stats

    if cells is not None:
        table = table.loc[cells].copy()
    group1 = table[table.loc[:, xlabel] == id1].copy()
    group2 = table[table.loc[:, xlabel] == id2].copy()

    t, p = stats.ranksums(group1.loc[:, score], group2.loc[:, score])
    result = pd.DataFrame(columns=["wilcoxon_ranksum", "pval"])
    result.loc[0] = [t, p]

    return result


def generate_count_object(
    adata,
    hue: str = "disease",
    cell_type_label: str = "cell_type",
    cell_type: List[str] = None,
    min_samples: int = 2,
    min_cells: int = 5,
    ref: str = "healthy",
    subset: List[str] = None,
    layer: str = "counts",
    outliers_removal: bool = False,
):
    """
    @Meshal what is this really supposed to do?

    Args:
        adata: AnnData object
        hue: Value to color by
        cell_type_label: Label containing cell types
        cell_type: Cells type to generate counts for
        min_samples: Minimum samples for outlier removal with DBScan
        min_cells: Minimal number of cells
        ref:
        subset:
        layer:
        outliers_removal: Whether to remove outliers or not

    Returns:
        AnnData object containing counts

    Example Call:
    subset = ['3d PI-KO', '3d PI-WT']

    raw_counts = generate_count_object(adata,
                                       condition = "grouping",
                                       cell_type_label = "celltype_refined", cell_type = ["AT2"],
                                       ref = "3d PI-WT",
                                       subset = subset)
    """
    adata_subset = adata[adata.obs.grouping.isin(subset)]
    cells = [
        True if (adata_subset.obs[cell_type_label][i] in cell_type) else False for i in range(adata_subset.n_obs)  # type: ignore
    ]

    # Raw count data for diffxpy
    obs = adata_subset[cells].obs.copy()
    var = adata_subset.var_names.copy()

    adata_raw = sc.AnnData(X=adata_subset[cells].layers[layer].copy())
    adata_raw.obs = obs
    adata_raw.var.index = var
    adata_raw.obsm = adata_subset[cells].obsm.copy()

    # Also automate tidy up with DBScan :)
    if outliers_removal:
        adata_raw.obs["dcluster"] = remove_outliers(adata_raw.obsm["X_umap"], min_samples=min_samples)
        sc.pl.umap(adata_raw, color=[hue, "dcluster"])
        adata_raw = adata_raw[adata_raw.obs.dcluster == "0"].copy()

    sc.pp.filter_genes(adata_raw, min_cells=min_cells)

    # Set reference as first column
    adata_raw.obs.loc[:, hue].cat.reorder_categories([ref, np.setdiff1d(subset, ref)[0]], inplace=True)
    pal = adata_subset.uns[f"{hue}_colors"]
    sc.pl.umap(adata_raw, color=[hue], palette=list(pal))

    return adata_raw


def tidy_de_table(de_test, adata, cells, ids=None, qval_thresh: float = 0.9, group_by: str = "treatment", cols=None):
    """
    Sorts diffxpy de table and adds percentages of expression per group

    Args:
        de_test: diffxpy de test
        adata: AnnData object
        cells:
        ids:
        qval_thresh:
        group_by:
        cols:

    Returns:
        Pandas Dataframe of diffxpy table with percentages
    """
    result = de_test.summary().sort_values(by=["qval"], ascending=True)
    result = result[result.qval < qval_thresh].loc[:, cols].copy()

    # Add percentages
    result = add_percentages(adata[cells], result, ids=ids, group_by=group_by)

    return result


def correlate_means_to_gene(means: pd.DataFrame, corr_gene: str = "EOMES"):
    """
    Calculate gene to gene correlation based on a mean expression table

    Args:
        means:
        corr_gene:

    Returns:
        Pandas DataFrame of correlations
    """
    import scipy.stats

    genes = means.columns.values
    cors = pd.DataFrame(index=genes, columns=["spearman_corr", "pvalue"])
    # tab = sc.get.obs_df(sub, keys = [corr_gene], layer = None, use_raw = True)
    table = means.loc[:, [corr_gene]].values

    # Loop over all genes.
    for gene in genes:
        tmp = scipy.stats.spearmanr(table, means.loc[:, [gene]])  # Spearman's rho
        cors.loc[gene, :] = tmp[0:2]

    cors.dropna(axis=0, inplace=True)
    cors.sort_values("spearman_corr", ascending=False, inplace=True)

    return cors


def automated_marker_annotation(
    adata: AnnData,
    organism: str,
    tissue: str,
    marker_file: str,
    key: str = "rank_genes_groups",
    normalize: Optional[Literal["reference", "data"]] = "reference",
    p_value: float = 0.05,
    log_fold_change: float = 2,
):
    """
    Calculates a marker gene overlap based on pre-existing annotations.

    Args:
        adata: AnnData object containing ranked genes
        organism: Currently supported: 'mouse'
        tissue: Currently supported: 'lung'
        marker_file: Currently supported: 'lung_particle_markers.txt'
        key: Key of ranked genes in adata (default: 'rank_genes_groups')
        normalize: Normalization option for the marker gene overlap output (default: 'reference')
        p_value: p-value threshold for existing marker genes (default: 0.05)
        log_fold_change: log fold change threshold for existing marker genes (default: 2)

    Returns:
        Pandas DataFrame of overlapping genes. Visualize with a Seaborn Heatmap
    """
    supported_organisms = {"mouse"}
    supported_tissues = {"lung"}
    supported_marker_files = {"lung_particle_markers.txt"}

    if organism not in supported_organisms:
        print(f"[bold red]Unfortunately organism {organism} is not yet supported.")
        return

    if tissue not in supported_tissues:
        print(f"[bold red]Unfortunately tissue {tissue} is not yet supported.")
        return

    if marker_file not in supported_marker_files:
        print(f"[bold red]Unfortunately marker file {marker_file} could not be found. Please check your spelling.")
        return

    marker_table = pd.read_csv(f"{WORKING_DIRECTORY}/markers/{marker_file}", sep="\t", index_col=None)
    marker_table = marker_table[
        (marker_table.logfoldchange > log_fold_change) & (marker_table.pval_adj < p_value)
    ].copy()

    marker = dict()
    for ct in marker_table["cluster"].unique():
        tmp = marker_table[marker_table["cluster"] == ct]
        marker[ct] = tmp.gene.values

    return sc.tl.marker_gene_overlap(adata, marker, key=key, normalize=normalize)
