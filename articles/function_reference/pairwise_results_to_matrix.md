# pairwise_results_to_matrix()

## Description

Builds a cluster-by-cluster adjacency matrix from
`adata.uns["repspat_mmd_results"]`, connecting clusters whose adjusted
p-value is at least `alpha` (i.e. clusters that could not be
statistically distinguished). Optionally draws the corresponding network
graph.

## Arguments

| Argument | Description |
|----|----|
| `adata` | AnnData object containing MMD^2 comparison results (see [`multiple_comparison()`](https://amanpreet60.github.io/repSpat_Workshop/articles/function_reference/multiple_comparison.md)), or a DataFrame with the same columns. |
| `plot` | Whether to draw the network graph in addition to returning the matrix. Default `True`. |
| `alpha` | Significance threshold used to decide which clusters are linked. Default `0.05`. |

## Value

A cluster-by-cluster adjacency `DataFrame`. Clusters considered similar
(repeated patterns) are connected with their observed MMD^2 as the edge
weight; clusters found to be significantly different are left
unconnected (`0`).

## Example

``` python
matrix = pairwise_results_to_matrix(adata)
```

[← Back to reference
index](https://amanpreet60.github.io/repSpat_Workshop/articles/function_reference/index.md)
