# spatial_silhouette_analysis()

## Description

Runs a spatially-aware silhouette analysis over combinations of
k-nearest spatial neighbors and cluster counts. For each combination,
cells are clustered with spatially constrained hierarchical clustering,
then a silhouette-like score is computed comparing each cell to its own
cluster and to spatially adjacent clusters only. Stores a DataFrame of
scores in `adata.uns["silhouette_scores"]` and returns the same `adata`.

## Arguments

| Argument | Description |
|----|----|
| `adata` | AnnData object with feature distances already computed (see [`compute_distances()`](https://amanpreet60.github.io/repSpat_Workshop/articles/function_reference/compute_distances.md)). |
| `n_neighbors_list` | List of neighbor counts to test, e.g. `[6, 8]`. |
| `n_clusters_range` | Range of cluster counts to test, e.g. `range(4, 9)`. |
| `linkage` | Linkage method: `"ward"`, `"single"`, `"complete"`, or `"average"`. Default `"ward"`. |

## Value

The same AnnData object, updated in place, with a DataFrame of scores
for every `(n_neighbors, n_clusters)` combination stored in
`adata.uns["silhouette_scores"]`.

## Example

``` python
adata = spatial_silhouette_analysis(
    adata,
    n_neighbors_list=[6, 8],
    n_clusters_range=range(4, 9),
)
adata.uns["silhouette_scores"]
```

[← Back to reference
index](https://amanpreet60.github.io/repSpat_Workshop/articles/function_reference/index.md)
