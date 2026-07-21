# compute_distances()

    compute_distances(adata, layer, metric, distance_key="repspat_distances")

Compute pairwise feature distances from a required AnnData layer.

## Description

Computes pairwise distances from the required AnnData layer, stores them
in `adata.obsp[distance_key]`, records the active feature layer in
`adata.uns["repspat"]["layer"]`, and returns the same `adata`.

## Arguments

| Argument | Description |
|----|----|
| `adata` | AnnData object with the feature matrix of interest stored in `adata.layers`. |
| `layer` | Name of the `adata.layers` entry to use as the feature matrix. Required. |
| `metric` | Distance metric passed to `scipy.spatial.distance.pdist` (e.g. `"euclidean"`, `"jaccard"`). |
| `distance_key` | Name under which the resulting distance matrix is stored in `adata.obsp`. Defaults to `"repspat_distances"`. |

## Value

The same AnnData object, updated in place, with the distance matrix
stored in `adata.obsp[distance_key]`.

## Example

``` python
from repspat import compute_distances

adata = compute_distances(adata, layer="layer 1", metric="euclidean")
```

[← Back to reference
index](https://amanpreet60.github.io/repSpat_Workshop/articles/function_reference/index.md)
