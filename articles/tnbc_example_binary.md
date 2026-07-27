# TNBC Example (Binary): Detecting Repeated Spatial Patterns

## Introduction

This is the binary version of the TNBC workflow. It follows the same
steps as the main example, but instead of using the continuous marker
measurements, we threshold each marker into a present (`1`) or absent
(`0`) value and compare cells with the Jaccard metric. See the main TNBC
vignette for the detailed explanation of each step.

## Load libraries

``` r

library(BiocStyle)
library(SpatialExperiment)
library(reticulate)
```

``` r

anndata <- import("anndata", convert = FALSE)
repspat <- import("repspat", convert = FALSE)
```

## Prepare the Sample

We convert the `.rds` file into an AnnData object, threshold each marker
into a present/absent value, and store the result in a `"binary"` layer
alongside the original measurements in `X`. A marker is present when its
value is equal to or greater than its threshold.

``` r

rds_path <- "../inst/extdata/03_TNBC_2018_spe.rds"
h5ad_path <- "../inst/extdata/03_TNBC_2018_spe_binary.h5ad"

load(rds_path)

expression <- t(as.matrix(assay(spe, "exprs")))

threshold_values <- unlist(thresholds, use.names = TRUE)
binary_matrix <- sweep(expression,2,threshold_values[colnames(expression)],FUN = ">=") + 0L

cell_metadata <- as.data.frame(colData(spe), optional = TRUE)
rownames(cell_metadata) <- colnames(spe)

feature_metadata <- data.frame(
    feature_name = rownames(spe),
    row.names = rownames(spe)
)

adata <- anndata$AnnData(
    X = expression,
    obs = cell_metadata,
    var = feature_metadata
)

py_set_item(adata$layers, "binary", binary_matrix)
py_set_item(adata$obsm, "spatial", spatialCoords(spe)[, 1:2, drop = FALSE])

adata$write_h5ad(h5ad_path)
adata
```

## Compute Feature Distances

We select `Sample_04` and compute pairwise distances from the `"binary"`
layer using the `"jaccard"` metric, which suits binary present/absent
data.

``` r

py$adata <- adata
adata <- py_eval("adata[adata.obs['sample_id'] == 'Sample_04'].copy()", convert = FALSE)

adata <- repspat$compute_distances(
    adata,
    layer = "binary",
    metric = "jaccard"
)
adata$obs$sample_id$unique()
```

``` r

cat("Feature matrix: first 5 cells\n")
print(py_to_r(py_get_item(adata$layers, "binary"))[1:5, 1:11])

cat("\nSpatial coordinates: first 5 cells\n")
print(head(py_to_r(py_get_item(adata$obsm, "spatial")), 5))

cat("\nDistance matrix: first 5 cells x first 5 cells\n")
print(py_to_r(py_get_item(adata$obsp, "repspat_distances"))[1:5, 1:5])
```

## Select Clustering Parameters

We test combinations of neighborhood size and cluster count and score
each with the spatial silhouette score.

``` r

data <- repspat$spatial_silhouette_analysis(
    adata = adata,
    n_neighbors_list = list(as.integer(6), as.integer(8)),
    n_clusters_range = as.integer(4:8)
)
cat("\n--------------------\n\n")
py_get_item(data$uns, "silhouette_scores")
```

We proceed with `n_neighbors = 8` and `n_clusters = 7`.

## Construct Spatially Constrained Clusters

We apply constrained agglomerative hierarchical clustering (CAHC) with
the chosen parameters.

``` r

adata <- repspat$spatial_constrained_hac(
    adata = adata,
    n_clusters = as.integer(7),
    n_neighs = as.integer(8)
)
cat("\n--------------------\n\n")
py_get_item(data$obs, "labels")
```

## Visualize the Spatially Constrained Clusters

``` r

plt <- reticulate::import("matplotlib.pyplot")
plot_result <- repspat$plot_spatial_clusters(
    adata = adata,
    label_key = "labels",
    figsize = reticulate::tuple(10, 7),
    point_size = 4L
)
plt$subplots_adjust(right = 0.72)
plt$show()
```

## Visualize Features Across Clusters

For binary data, this shows the features most frequently present in each
cluster.

``` r

plt <- reticulate::import("matplotlib.pyplot")
feature_plots <- repspat$plot_cluster_feature_presence(
    adata = data,
    label_key = "labels",
    top_n = as.integer(10),
    figsize = reticulate::tuple(8, 3),
    point_size = 2L,
    alpha = 1.0
)

for (fig_axes in reticulate::iterate(feature_plots$values())) {
    fig <- py_get_item(fig_axes, 0L)
    plt$figure(fig$number)
    plt$show()
}
```

## Create blocks for permutation

We group cells into blocks within each cluster so permutation happens at
the block level.

``` r

adata <- repspat$create_blocks(
    adata = adata,
    knn = as.integer(8)
)

cat("\n--------------------\n\n")
py_get_item(data$obs, "repspat_block_id")
```

## Test Spatial Invariance Between Clusters

We compare every pair of clusters with the MMD^2 statistic, using the
`"IMQ"` kernel, 200 block permutations, and Benjamini–Hochberg
correction.

``` r

adata <- repspat$multiple_comparison(
    adata = adata,
    kernel = "IMQ",
    nperm = as.integer(200),
    adj_p = "BH"
)

cat("\n--------------------\n\n")
py_to_r(py_get_item(adata$uns, "repspat_mmd_results")$head(5L)$drop(columns = list("null_dist")))
```

## Visualize Repeated Spatial Patterns

Cluster pairs with `adj_p >= 0.05` are repeated spatial patterns and are
drawn as a network.

``` r

similarity_matrix <- py_to_r(repspat$pairwise_results_to_matrix(
    adata,
    plot = TRUE
))

similarity_matrix
```

## Relabel Repeated Spatial Patterns

We relabel the clusters that could not be statistically distinguished
with a single shared label.

``` r

clusters_to_merge <- c(3, 5, 6)   # edit: select labels to reassign
new_label <- min(clusters_to_merge)
labels <- py_to_r(py_get_item(adata$obs, "labels"))
labels[labels %in% clusters_to_merge] <- new_label
py_set_item(adata$obs, "updated_labels", as.integer(labels))
```

## Visualize the Relabeled Repeated Spatial Patterns

``` r

plt <- reticulate::import("matplotlib.pyplot")
plot_result <- repspat$plot_spatial_clusters(
    adata = adata,
    label_key = "updated_labels",
    figsize = reticulate::tuple(10, 7),
    point_size = 4L
)
plt$subplots_adjust(right = 0.72)
plt$show()
```

## Conclusion

Applied to the TNBC sample using binary marker profiles and the Jaccard
metric, `repSpat` recovered repeated patterns among the clusters. As
with the continuous version, repeated patterns are best interpreted
alongside cluster size, location, and biological context rather than the
test result in isolation.
