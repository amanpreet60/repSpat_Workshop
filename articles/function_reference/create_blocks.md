# create_blocks()

    create_blocks(adata, knn)

Create permutation blocks within each region for block-based testing.

## Description

Splits each region into smaller KMeans-based blocks so that cells within
the same tissue region are not treated as fully independent observations
during permutation testing. Stores block IDs in
`adata.obs["repspat_block_id"]`, prints the storage location, and
returns the same `adata`.

## Arguments

| Argument | Description |
|----|----|
| `adata` | AnnData object with region labels already assigned (see [`spatial_constrained_hac()`](https://amanpreet60.github.io/repSpat_Workshop/articles/function_reference/spatial_constrained_hac.md)). |
| `knn` | Approximate number of cells per block. |
| `region_key` | Column in `adata.obs` containing the region labels. Default `"labels"`. |
| `block_key` | Column in `adata.obs` to store the resulting block IDs. Default `"repspat_block_id"`. |

## Value

The same AnnData object, updated in place, with the block assignment for
every cell stored in `adata.obs[block_key]`.

## Example

``` python
adata = create_blocks(adata, knn=8)
```

[← Back to reference
index](https://amanpreet60.github.io/repSpat_Workshop/articles/function_reference/index.md)
