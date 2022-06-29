from __future__ import annotations

import os
from typing import List

from pandas import Categorical
from statsmodels.stats.multitest import fdrcorrection

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal  # type: ignore

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
    xlabel: str = None,
    condition: str = None,
    use_raw: bool = None,
):
    """Generates a table of cells by genes of expression values as a Pandas DataFrame.

    Args:
        adata: Anndata object
        cluster: Which label of the subsets to generate the table for. Use 'all' if for all subsets.
        subset_by: Which label to subset the clusters by
        xlabel: Label that will be used for subsequent line plots as x-axis label. Typically a time series such as "days".
        condition: Column name of the condition to include.
        use_raw: Whether to use adata.raw.X for the calculations

    Returns:
         Gene expression table.
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
    if xlabel:
        gen_expression_table[xlabel] = adata[cells].obs[xlabel]
    if condition:
        # For multiple cluster, split internally per condition
        if isinstance(cluster, list) and len(cluster) > 1 and subset_by != condition:
            gen_expression_table[condition] = [
                f"{t}_{c}" for t, c in zip(adata[cells].obs[condition], adata[cells].obs[subset_by])
            ]

        else:
            gen_expression_table[condition] = adata[cells].obs[condition]

    return gen_expression_table


def relative_frequencies(adata, group_by: str = "cell_type", xlabel: str = "days", condition: str = "batch"):
    """Calculates the relative frequencies of conditions grouped by an observation.

    Args:
        adata: AnnData Objet containing the data
        group_by: Column name to group by
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


def remove_outliers(cords, eps: int = 1, min_samples: int = 2) -> Categorical:
    """Remove outlying cells based on UMAP embeddings with DBScan (density based clustering).

    Call as: sub.obs["d_cluster"] = remove_outliers(sub.obsm["X_umap"], min_samples = 10)

    Args:
        cords: adata UMAP coordinates, typically adata.obsm["X_umap"]
        eps: Maximum distance between two clusters to still be considered neighbors
        min_samples: Minimum samples of a cluster

    Returns:
        Pandas Categorical of clusters
    """
    from natsort import natsorted
    from sklearn.cluster import DBSCAN

    clustering = DBSCAN(eps=eps, min_samples=min_samples).fit(cords)
    cluster = clustering.labels_.astype("U")

    return pd.Categorical(cluster, categories=natsorted(np.unique(cluster)))


def add_percentages(adata, table, ids, group_by: str, threshold: int = 0, gene_label: str = "gene"):
    """Add columns to existing diffxpy table specifying percentage of expressing cells.

    Args:
        adata: AnnData object containing the data
        table: Table as generated by diffxpy
        ids: Identifiers to add percentages for.
        group_by: Label to group by
        threshold: Cell count threshold.
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
            table[f"pct.{ident}s"] = data_temp.reindex(table.index.values).values
        else:
            table[f"pct.{ident}s"] = data_temp.reindex(table.loc[:, gene_label]).values
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
        True if (adata_subset.obs[cell_type_label][i] in cell_type) else False  # type: ignore
        for i in range(adata_subset.n_obs)
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


def extended_marker_table(
    adata: AnnData,
    qval_thresh: float = 0.05,
    cell_type_label: str = "cell_type",
    gene_ranks_key: str = "rank_genes_groups",
):
    """
    Generates an extended marker table with cell types and percentages of expressed cell types per cluster.

    Run scanpy.tl.rank_genes_groups before using this function.

    Args:
        adata: AnnData object containing ranked genes
        qval_thresh: Threshold to filter the log fold change for
        cell_type_label: Label containing all cell types
        gene_ranks_key: Key for the ranked gene groups (generated by sc.tl.rank_genes_groups)

    Returns:
        A Pandas DataFrame
    """
    result = adata.uns[gene_ranks_key]
    all_markers = []

    for cluster in result["names"].dtype.names:
        current = pd.DataFrame(
            {
                "gene": result["names"][cluster],
                "score": result["scores"][cluster],
                "log_FC": result["logfoldchanges"][cluster],
                "pval": result["pvals"][cluster],
                "pval_adj": result["pvals_adj"][cluster],
                "cell_type": cluster,
            }
        )

        # Add percentage expressed per cell type
        adata.obs["group"] = ["within" if ct == cluster else "outside" for ct in adata.obs.loc[:, cell_type_label]]
        current = add_percentages(adata, table=current, group_by="group", gene_label="gene", ids=["within", "outside"])

        all_markers.append(current)

    all_markers_df = pd.concat(all_markers)
    all_markers_df = all_markers_df[all_markers_df.pval_adj < qval_thresh].copy()

    return all_markers_df


def generate_pseudobulk(adata: AnnData, group_key: str = "identifier", sep="\t", save: str = None) -> pd.DataFrame:
    """
    Generates a pseudobulk for a given key of groups in the AnnData object.

    Looks like:

    +------------+------------------+------------------+
    | Genes      | Group Member 1   | Group Member 2   |
    +============+==================+==================+
    | Gene 1     | Value 1          | Value 2          |
    +------------+------------------+------------------+
    | Gene 2     | Value 2          | Value 3          |
    +------------+------------------+------------------+

    Args:
        adata: AnnData object
        group_key: The key to group by. E.g. by mice, by condition, ... (default: 'identifier')
        sep: Separator to use when saving the pseudobulk table (default: '\t')
        save: Path to save the pseudobulk table to (default: None)

    Returns:
        A Pandas DataFrame containing the pseudobulk table
    """
    pseudobulk = pd.DataFrame(data=adata.var_names.values, columns=["Genes"])

    for i in adata.obs.loc[:, group_key].cat.categories:
        temp = adata.obs.loc[:, group_key] == i
        pseudobulk[i] = adata[temp].X.sum(0, dtype=int)  # column sums (genes)

    if save:
        pseudobulk.to_csv(save, sep=sep, index=False)

    return pseudobulk


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
    """Calculates a marker gene overlap based on pre-existing annotations.

    Currently supported marker files:

    +------------+------------+------------------------------+
    | Organism   | Tissue     | Marker File                  |
    +============+============+==============================+
    | Mouse      | Lung       | lung_particle_markers.txt    |
    +------------+------------+------------------------------+
    | Human      | NA         |                              |
    +------------+------------+------------------------------+

    Args:
        adata: AnnData object containing ranked genes
        organism: Currently supported: 'mouse'
        tissue: Currently supported: 'lung'
        marker_file: Name of the marker file to be used - refer to table
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
    for ct in marker_table["cell_type"].unique():
        tmp = marker_table[marker_table["cell_type"] == ct]
        marker[ct] = tmp.gene.values

    return sc.tl.marker_gene_overlap(adata, marker, key=key, normalize=normalize)


def de_res_to_anndata(
    adata: AnnData,
    de_res: pd.DataFrame,
    *,
    groupby: str,
    gene_id_col: str = "gene_symbol",
    score_col: str = "score",
    pval_col: str = "pvalue",
    pval_adj_col: Optional[str] = None,
    lfc_col: str = "lfc",
    key_added: str = "rank_genes_groups",
) -> None:
    """Add a tabular differential expression result to AnnData as if it was produced by scanpy.tl.rank_genes_groups.

    Args:
        adata: Annotated data matrix
        de_res: Tablular DE result as Pandas DataFrame
        groupby: Column in `de_res` that indicates the group. This column must also exist in `adata.obs`.
        gene_id_col: Column in `de_res` that holds the gene identifiers
        score_col: Column in `de_res` that holds the score (results will be ordered by score).
        pval_col: Column in `de_res` that holds the unadjusted pvalue
        pval_adj_col: Column in `de_res` that holds the adjusted pvalue.
                      If not specified, the unadjusted p values will be FDR-adjusted.
        lfc_col: Column in `de_res` that holds the log fold change
        key_added: Key under which the results will be stored in adata.uns
    """
    if groupby not in adata.obs.columns or groupby not in de_res.columns:
        raise ValueError("groupby column must exist in both adata and de_res. ")
    res_dict = {
        "params": {
            "groupby": groupby,
            "reference": "rest",
            "method": "other",
            "use_raw": True,
            "layer": None,
            "corr_method": "other",
        },
        "names": [],
        "scores": [],
        "pvals": [],
        "pvals_adj": [],
        "logfoldchanges": [],
    }
    df_groupby = de_res.groupby(groupby)
    for _, tmp_df in df_groupby:
        tmp_df = tmp_df.sort_values(score_col, ascending=False)
        res_dict["names"].append(tmp_df[gene_id_col].values)  # type: ignore
        res_dict["scores"].append(tmp_df[score_col].values)  # type: ignore
        res_dict["pvals"].append(tmp_df[pval_col].values)  # type: ignore
        if pval_adj_col is not None:
            res_dict["pvals_adj"].append(tmp_df[pval_adj_col].values)  # type: ignore
        else:
            res_dict["pvals_adj"].append(fdrcorrection(tmp_df[pval_col].values)[1])  # type: ignore
        res_dict["logfoldchanges"].append(tmp_df[lfc_col].values)  # type: ignore

    for key in ["names", "scores", "pvals", "pvals_adj", "logfoldchanges"]:
        res_dict[key] = pd.DataFrame(
            np.vstack(res_dict[key]).T,  # type: ignore
            columns=list(df_groupby.groups.keys()),
        ).to_records(index=False, column_dtypes="O")
    adata.uns[key_added] = res_dict
