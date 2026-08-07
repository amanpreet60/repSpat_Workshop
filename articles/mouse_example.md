# Mouse CosMx Example

## Introduction

This is the mouse CosMx version of the workflow. It follows the same
steps as the main TNBC example, using continuous features and the
Euclidean metric. The one difference is in preprocessing: because the
full object is large (~48,556 cells), only a small contiguous crop of
the tissue (~5,000 cells) is retained for the analysis. For a detailed
explanation of each step, see the main TNBC vignette.

## Load libraries

``` r

library(BiocStyle)
library(reticulate)
```

``` r

anndata <- import("anndata", convert = FALSE)
repspat <- import("repspat", convert = FALSE)
```

## Prepare the Sample

Because the mouse object is a single tissue, rather than selecting one
sample the spatial coordinates are cropped to a rectangular window
centered on the tissue and sized to retain roughly 5,000 cells. Since
the crop is contiguous, local neighborhoods and cell density are
preserved. The feature matrix is stored in a `"layer 1"` layer and a
single `sample_id` label is added. Edit `target_cells` to change the
crop size.

``` r

adata_full <- anndata$read_h5ad(brain_path)

coords <- py_to_r(py_get_item(adata_full$obsm, "spatial"))[, 1:2, drop = FALSE]
x <- coords[, 1]
y <- coords[, 2]

# Size a square window, centered on the tissue, to capture ~target_cells.
n_total <- nrow(coords)
frac    <- min(1, target_cells / n_total)
cx      <- median(x)
cy      <- median(y)
half_x  <- (diff(range(x)) / 2) * sqrt(frac)
half_y  <- (diff(range(y)) / 2) * sqrt(frac)

in_window <- x >= (cx - half_x) & x <= (cx + half_x) &
             y >= (cy - half_y) & y <= (cy + half_y)

cat("Cells kept:", sum(in_window), "of", n_total, "\n")

py$adata_full <- adata_full
py$in_window <- r_to_py(which(in_window) - 1L)   # 0-based indices for Python
adata <- py_eval("adata_full[in_window].copy()", convert = FALSE)

# Store the feature matrix as a named layer and add a sample label.
py_set_item(adata$layers, "layer 1", adata$X)
n_cells <- as.integer(py_to_r(adata$n_obs))
py_set_item(adata$obs, "sample_id", r_to_py(rep("mouse_brain", n_cells)))

adata$write_h5ad(h5ad_path, compression = "gzip")
adata
```

## Compute Feature Distances

The first step of the workflow selects the cropped tissue and computes
pairwise feature distances between its cells from the `"layer 1"` layer,
using the `"euclidean"` metric. These distances are used in the
subsequent steps to construct spatially constrained clusters.

``` r

py$adata <- adata
adata <- py_eval("adata[adata.obs['sample_id'] == 'mouse_brain'].copy()", convert = FALSE)

adata <- repspat$compute_distances(
    adata,
    layer = "layer 1",
    metric = "euclidean"
)
adata$obs$sample_id$unique()
```

``` r

cat("Feature matrix: first 5 cells\n")
print(py_to_r(py_get_item(adata$layers, "layer 1"))[1:5, 1:11])

cat("\nSpatial coordinates: first 5 cells\n")
print(head(py_to_r(py_get_item(adata$obsm, "spatial")), 5))

cat("\nDistance matrix: first 5 cells x first 5 cells\n")
print(py_to_r(py_get_item(adata$obsp, "repspat_distances"))[1:5, 1:5])
```

## Select Clustering Parameters

Before clustering, combinations of neighborhood size and number of
clusters are evaluated and scored based on the spatial silhouette score.

``` r

adata <- repspat$spatial_silhouette_analysis(
    adata = adata,
    n_neighbors_list = list(as.integer(6), as.integer(8)),
    n_clusters_range = as.integer(4:8)
)
cat("\n--------------------\n\n")
py_get_item(adata$uns, "silhouette_scores")
```

For this crop the scores are highest at a small number of clusters, so
we proceed with `n_neighbors = 8` and `n_clusters = 7` following
clustering step.

## Construct Spatially Constrained Clusters

Using the selected parameters, the cropped tissue is partitioned into
spatially constrained clusters with constrained agglomerative
hierarchical clustering (CAHC), which groups cells by feature similarity
while enforcing spatial connectivity. The resulting cluster labels are
stored in `adata.obs["labels"]`.

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

For continuous features, it displays those with the largest standardized
mean difference between the cluster and all remaining cells.

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
grouped into smaller blocks.

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
the Benjamini–Hochberg procedure for multiple comparisons.

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
treated as repeated spatial patterns.

``` r

similarity_matrix <- py_to_r(repspat$pairwise_results_to_matrix(
    adata,
    plot = TRUE
))

similarity_matrix
```

For this dataset, every cluster pair was found to be significantly
different (`adj_p < 0.05`). No repeated spatial patterns were detected.
Thus, we do not relabel the clusters.

## Conclusion

The application of repSpat to a cropped region of the mouse CosMx sample
showed statistically significant differences in feature distributions
among all seven spatially constrained clusters. Thus, no repeated
spatial patterns were identified within this region. The absence of
repeated spatial patterns is also informative, indicating that none of
the identified clusters could be relabelled together based on the
repSpat pairwise tests.

For the whole tissue, the tissue can be divided into contiguous regions,
and CAHC can be applied to these regions in parallel. The resulting
clusters can then be collected, and the repSpat pairwise tests can be
performed across all clusters to identify repeated spatial patterns in
the entire tissue.
