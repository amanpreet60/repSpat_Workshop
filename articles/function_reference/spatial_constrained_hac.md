# spatial_constrained_hac()

## Description

Hierarchical agglomerative clustering with a spatial connectivity
constraint: cells can only be grouped together when they are connected
through nearby cells in the tissue, keeping clusters spatially coherent.
Stores labels in `adata.obs["labels"]`, prints the storage location, and
returns the same `adata`. Clustering uses the active feature layer
recorded by
[`compute_distances()`](https://amanpreet60.github.io/repSpat_Workshop/articles/function_reference/compute_distances.md).

## Arguments

| Argument | Description |
|----|----|
| `adata` | AnnData object with feature distances already computed. |
| `n_clusters` | Number of clusters to construct. Default `7`. |
| `n_neighs` | Number of nearby cells used to define each spatial neighborhood. Default `8`. |
| `coord_type` | Coordinate type passed to `squidpy.gr.spatial_neighbors`. Default `"generic"`. |
| `delaunay` | Whether to build the spatial graph using Delaunay triangulation. Default `False`. |
| `linkage` | Linkage method: `"ward"`, `"single"`, `"complete"`, or `"average"`. Default `"ward"`. |

## Value

The same AnnData object, updated in place, with the assigned cluster
label for every cell stored in `adata.obs["labels"]`.

## Example

``` python
adata = spatial_constrained_hac(adata, n_clusters=7, n_neighs=8, linkage="ward")
```
