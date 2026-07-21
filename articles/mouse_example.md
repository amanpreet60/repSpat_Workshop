# Mouse Visium Example: Detecting Repeated Spatial Patterns

## Introduction

This is the mouse Visium version of the workflow. It follows the same
steps as the main TNBC example, using continuous features and the
Euclidean metric. The one difference is preprocessing: the full object
is large (~48,556 cells), so we keep only a small contiguous crop of the
tissue (~5,000 cells). See the main TNBC vignette for the detailed
explanation of each step.

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

The mouse object is a single tissue, so instead of selecting one sample
we crop the spatial coordinates to a rectangular window centered on the
tissue, sized to keep roughly 5,000 cells. Because it is a contiguous
crop, local neighborhoods and density are preserved. The feature matrix
is stored in a `"layer 1"` layer and a single `sample_id` label is
added. Edit `target_cells` to change the crop size.

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

adata$write_h5ad(h5ad_path)
adata
```

## Compute Feature Distances

We select the cropped tissue and compute pairwise distances from the
`"layer 1"` layer using the `"euclidean"` metric.

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

## Select Spatial Region Settings

We test combinations of neighborhood size and region count and score
each with the custom silhouette score.

``` r

data <- repspat$spatial_silhouette_analysis(
    adata = adata,
    n_neighbors_list = list(as.integer(6), as.integer(8)),
    n_clusters_range = as.integer(4:8)
)
cat("\n--------------------\n\n")
py_get_item(data$uns, "silhouette_scores")
```

For this crop the scores are highest at a small number of clusters, so
we proceed with `n_neighbors = 8` and `n_clusters = 7`.

## Divide the Tissue Into Spatial Regions

We apply spatially constrained hierarchical clustering with the chosen
settings.

``` r

adata <- repspat$spatial_constrained_hac(
    adata = adata,
    n_clusters = as.integer(7),
    n_neighs = as.integer(8)
)
cat("\n--------------------\n\n")
py_get_item(data$obs, "labels")
```

The crop is divided into 7 spatially coherent regions, each occupying a
connected part of the tissue.

## Visualize the Spatial Regions

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

## Visualize Features Across Regions

This shows the numeric features most enriched in each region compared to
the rest of the tissue.

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

## Create Permutation Blocks

We group cells into blocks within each region so permutation happens at
the block level.

``` r

adata <- repspat$create_blocks(
    adata = adata,
    knn = as.integer(8)
)

cat("\n--------------------\n\n")
py_get_item(data$obs, "repspat_block_id")
```

## Compare Region Pairs

We compare every pair of regions with MMD, using the `"IMQ"` kernel, 200
block permutations, and Benjamini–Hochberg correction.

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

Region pairs with a small `obs_mmd_sq` and a large `adj_p` have similar
feature profiles; those with a large MMD and small `adj_p` are the
clearly distinct regions.

## Visualize Which Regions Are Not Significantly Different

Region pairs with `adj_p >= 0.05` are repeated spatial patterns and are
drawn as a network.

``` r

similarity_matrix <- py_to_r(repspat$pairwise_results_to_matrix(
    adata,
    plot = TRUE
))

similarity_matrix
```

For this dataset, every region pair was significantly different
(`adj_p < 0.05`), so the network has no connections. No repeated spatial
patterns were detected, and there is nothing to merge.

## Conclusion

Applied to a cropped portion of the mouse Visium sample, `repSpat` found
all seven regions to be statistically distinct, so no repeated spatial
patterns were recovered here. This is itself an informative result:
within this small contiguous crop, the regions carry genuinely different
feature profiles. As with the TNBC example, results are best interpreted
alongside region size, location, and biological context rather than the
test result in isolation.
