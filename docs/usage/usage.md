# Usage

Import the sc-toolbox API as follows:

```python
import sc_toolbox as sct
```

You can then access the respective modules like:

```python
sct.pl.cool_fancy_plot()
```

```{eval-rst}
.. currentmodule:: sc_toolbox
```

## Preprocessing

```{eval-rst}
.. autosummary::
    :toctree: preprocessing
```

## Tools

```{eval-rst}
.. autosummary::
    :toctree: tools
    :nosignatures:

    tools.generate_expression_table
    tools.relative_frequencies
    tools.relative_frequency_per_cluster
    tools.correlate_to_signature
    tools.remove_outliers
    tools.add_percentages
    tools.ranksums_between_groups
    tools.generate_count_object
    tools.tidy_de_table
    tools.correlate_means_to_gene
    tools.extended_marker_table
    tools.generate_pseudobulk
    tools.automated_marker_annotation
    tools.de_res_to_anndata
```

## Plots

```{eval-rst}
.. autosummary::
    :toctree: plot

    plot.Colormaps
    plot.custom_plot_size
    plot.standard_lineplot
    plot.average_expression
    plot.average_expression_per_cluster
    plot.average_expression_split_cluster
    plot.average_expression_per_cell
    plot.gene_expression_dpt_ordered
    plot.relative_frequencies_boxplots
    plot.split_boxplot
    plot.marker_dendrogram
    plot.volcano_plot
    plot.cluster_composition_stacked_barplot
    plot.gene_boxplot
    plot.colors_overview
    plot.relative_frequencies_lineplot
    plot.annotated_cell_type_umap
    plot.genotype_vs_genotype_umaps
```
