# TNBC Example (Binary)

## Introduction

This is the binary version of the TNBC workflow. It follows the same
steps as the main example, except that each marker measurement is
thresholded into a present (`1`) or absent (`0`) value and cells are
compared using the Jaccard metric, which is suited to binary
present/absent data. For a detailed explanation of each step, see the
main TNBC vignette.

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

In this step, the `.rds` file is converted into an AnnData object and
each marker is thresholded into a present/absent value. The thresholded
values are stored in a `"binary"` layer alongside the original
measurements in `X`. A marker is considered present when its value is
greater than or equal to its threshold.

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

The first step of the workflow selects `Sample_04` and computes pairwise
feature distances between its cells from the `"binary"` layer, using the
`"jaccard"` metric suited to binary present/absent data. These distances
are used in the subsequent steps to construct spatially constrained
clusters.

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

Before clustering, combinations of neighborhood size and cluster count
are evaluated and scored with the spatial silhouette score, which
measures how well the resulting clusters are spatially separated.

``` r

adata <- repspat$spatial_silhouette_analysis(
    adata = adata,
    n_neighbors_list = list(as.integer(6), as.integer(8)),
    n_clusters_range = as.integer(4:8)
)
cat("\n--------------------\n\n")
py_get_item(adata$uns, "silhouette_scores")
```

We proceed with `n_neighbors = 8` and `n_clusters = 7` for the
clustering step.

## Construct Spatially Constrained Clusters

Using the selected parameters, the tissue is partitioned into spatially
constrained clusters with constrained agglomerative hierarchical
clustering (CAHC), which groups cells by feature similarity while
enforcing spatial connectivity. The resulting cluster labels are stored
in `adata.obs["labels"]`.

``` r

adata <- repspat$spatial_constrained_hac(
    adata = adata,
    n_clusters = as.integer(7),
    n_neighs = as.integer(8)
)
cat("\n--------------------\n\n")
py_get_item(adata$obs, "labels")
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

For binary features, this displays the features most frequently present
in each cluster.

``` r

plt <- reticulate::import("matplotlib.pyplot")
feature_plots <- repspat$plot_cluster_feature_presence(
    adata = adata,
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

To account for local spatial dependence, cells within each cluster are
grouped into smaller blocks so that permutation is performed at the
block level rather than the individual-cell level.

``` r

adata <- repspat$create_blocks(
    adata = adata,
    knn = as.integer(8)
)

cat("\n--------------------\n\n")
py_get_item(adata$obs, "repspat_block_id")
```

## Test Spatial Invariance Between Clusters

Every pair of clusters is compared using the MMD² statistic to test
whether their feature distributions differ. The test uses the `"IMQ"`
kernel, 200 block permutations to approximate the null distribution, and
the Benjamini–Hochberg correction for multiple comparisons.

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

Cluster pairs that are not significantly different (`adj_p >= 0.05`) are
treated as repeated spatial patterns and drawn as a network.

``` r

similarity_matrix <- py_to_r(repspat$pairwise_results_to_matrix(
    adata,
    plot = TRUE
))

similarity_matrix
```

## Relabel Repeated Spatial Patterns

Based on the test results, the clusters that could not be statistically
distinguished are assigned a single shared label. This step changes only
the cluster labels.

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

## Running Multiple Samples

Here we worked with just one sample, `Sample_04`, because every sample
is its own piece of tissue with its own layout. If you have more
samples, you can run the same steps on each one, one at a time, and then
compare the repeated spatial patterns you find across them.

## Conclusion

Applied to the TNBC sample using binary marker profiles and the Jaccard
metric, `repSpat` recovered repeated patterns among the clusters. As
with the continuous version, repeated patterns are best interpreted
alongside cluster size, location, and biological context rather than the
test result in isolation.
