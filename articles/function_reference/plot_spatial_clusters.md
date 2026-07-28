# plot_spatial_clusters()

## Description

Draws a spatial scatter plot of all cells, colored by their assigned
cluster label.

## Arguments

| Argument | Description |
|----|----|
| `adata` | AnnData object with cluster labels assigned (see [`spatial_constrained_hac()`](https://amanpreet60.github.io/repSpat_Workshop/articles/function_reference/spatial_constrained_hac.md)). |
| `label_key` | Column in `adata.obs` containing the cluster labels to plot. Default `"labels"`. |
| `figsize` | Figure size in inches. Default `(4, 4)`. |
| `point_size` | Marker size for each cell. Default `10`. |
| `alpha` | Marker transparency. Default `1.0`. |
| `title` | Plot title. |
| `show_legend` | Whether to draw the cluster legend. Default `True`. |

## Value

A `(figure, axes)` tuple from `matplotlib`.

## Example

``` python
plot_spatial_clusters(adata, label_key="labels", point_size=4)
```
