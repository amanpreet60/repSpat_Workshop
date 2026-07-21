# plot_cluster_feature_presence()

    plot_cluster_feature_presence(adata, top_n=5, feature_columns=None, figsize=(8, 3), point_size=2, alpha=1.0, label_key="labels")

Plot each region alongside its most defining features.

## Description

For every region, draws a spatial highlight next to a chart of its most
defining features. Binary features are shown as presence rates within
the region; if no binary features are available, numeric features are
shown as cluster-vs-rest standardized enrichment scores.

## Arguments

| Argument | Description |
|----|----|
| `adata` | AnnData object with region labels assigned. |
| `top_n` | Number of top features to display per region. Default `5`. |
| `feature_columns` | Optional list restricting which feature columns are considered. Defaults to all binary columns, or all numeric columns if none are binary. |
| `figsize` | Figure size per region, in inches. Default `(8, 3)`. |
| `point_size` | Marker size for the spatial highlight. Default `2`. |
| `alpha` | Marker transparency. Default `1.0`. |
| `label_key` | Column in `adata.obs` containing the region labels. Default `"labels"`. |

## Value

A dict mapping each region’s cluster ID to its `(figure, axes)` pair.

## Example

``` python
cluster_figures = plot_cluster_feature_presence(adata, top_n=10)
```

[← Back to reference
index](https://amanpreet60.github.io/repSpat_Workshop/articles/function_reference/index.md)
