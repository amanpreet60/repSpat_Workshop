# TNBC Example: Detecting Repeated Spatial Patterns

## Introduction

Spatial omics data tells us both **what is present in a tissue** and
**where it is located**. Sometimes, two areas that are far apart contain
similar cells or biological features. These areas may represent a
repeated spatial pattern. In this notebook, we will use `repSpat` to
find and compare these repeated spatial patterns in a TNBC dataset.

The workflow uses five main package components:

- `compute_distances()` computes pairwise feature distances from a
  selected AnnData layer.
- `spatial_silhouette_analysis()` compares combinations of k-nearest
  neighbors and cluster counts.
- `spatial_constrained_hac()` performs hierarchical aglomerative
  clustering based on the k and number of clusters.
- `create_blocks()` creates smaller spatial groups within each region
  for block-based permutation testing.
- `multiple_comparison()` compares pairs of regions to identify similar
  or distinct spatial patterns.

At the end, we identify region pairs for which the analysis does not
find strong evidence of a difference. These pairs are repeated spatial
patterns.

We will introduce one component at a time, explain its inputs and
outputs, run it, and inspect the result before continuing.

## Load libraries

We first load the R packages needed to format the notebook, work with
spatial omics data, and call the Python `repSpat` package from R.

``` r

library(BiocStyle)
library(SpatialExperiment)
library(reticulate)
```

This code imports the Python packages needed for the analysis: `anndata`
for working with AnnData objects and `repspat` for detecting repeated
spatial patterns.

``` r

anndata <- import("anndata", convert = FALSE)
repspat <- import("repspat", convert = FALSE)
```

## Prepare the Sample

Before using `repSpat`, we first preprocess the TNBC dataset. We start
with the `.rds` file, convert it into an AnnData object, and then select
one tissue sample for analysis.

The original dataset contains multiple samples, but each sample
represents a separate tissue layout. We therefore analyze one sample at
a time.

In this step, we:

- import the TNBC `.rds` file,
- extract the marker measurements, cell metadata, and spatial
  coordinates,
- create an AnnData object from these data,
- select `Sample_04` for the analysis.

``` r

rds_path <- "../inst/extdata/03_TNBC_2018_spe.rds"
h5ad_path <- "../inst/extdata/03_TNBC_2018_spe.h5ad"

load(rds_path)

expression <- t(as.matrix(assay(spe, "exprs")))
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

py_set_item(adata$layers, "layer 1", expression)
py_set_item(adata$obsm, "spatial", spatialCoords(spe)[, 1:2, drop = FALSE])

adata$write_h5ad(h5ad_path)
adata
```

The result of this preprocessing step is an AnnData object containing
only `Sample_04`. This prepared object is then used as the input for the
`repSpat` workflow.

## Compute Feature Distances

Before dividing the tissue into regions, we first select `Sample_04`
from the full dataset and compute pairwise feature distances between its
cells.

This step measures how similar or different cells are based on their
marker profiles, storing the result as a distance matrix that later
steps use to build spatial regions.

The inputs are:

- `adata`: the AnnData object, filtered to `Sample_04`.
- `layer = "layer 1"`: the feature matrix used to compute distances.
- `metric = "euclidean"`: the distance metric applied between cells.

``` r

py$adata <- adata
adata <- py_eval("adata[adata.obs['sample_id'] == 'Sample_04'].copy()", convert = FALSE)

adata <- repspat$compute_distances(
    adata,
    layer = "layer 1",
    metric = "euclidean"
)
adata$obs$sample_id$unique()
```

The resulting `adata` object contains the selected `Sample_04` data in
AnnData format.

The result includes:

- `layers["layer 1"]`: the main feature matrix.
- `obsm["spatial"]`: the spatial coordinates.
- `obs`: the cell metadata, including the sample labels and cell-type
  labels.
- `uns["repspat_distances"]`: the feature-based distance matrix computed
  by `repspat$compute_distances`.

We can check the size of each part:

``` r

cat("Feature matrix: first 5 cells\n")
print(py_to_r(py_get_item(adata$layers, "layer 1"))[1:5, 1:11])

cat("\nSpatial coordinates: first 5 cells\n")
print(head(py_to_r(py_get_item(adata$obsm, "spatial")), 5))

cat("\nDistance matrix: first 5 cells x first 5 cells\n")
print(py_to_r(py_get_item(adata$obsp, "repspat_distances"))[1:5, 1:5])
```

This confirms that the sample was loaded and its main inputs were
prepared successfully. The distances are now stored in `adata`, ready to
guide the choice of spatial region settings.

## Select Spatial Region Settings

Before dividing the tissue into regions, we need to choose two settings
and test which combination works best.

The analysis tests different combinations of neighborhood size and
region count and gives each combination a silhouette score. Higher
scores indicate better separation between spatial regions. The custom
silhouette score compares each cell with its own cluster and with
spatially neighboring clusters only. This makes the score focus on how
well nearby tissue regions are separated.

The ranges tested here can be adjusted based on the dataset and in
collaboration with domain experts who understand the tissue structure
and expected biological patterns.

The inputs are:

- `adata`: the AnnData object with feature distances already computed.
- `n_neighbors_list = [6, 8]`: how many nearby cells are used to define
  local spatial neighborhoods.
- `n_clusters_range = 4:8`: how many spatial regions we want to create.

``` r

data <- repspat$spatial_silhouette_analysis(
    adata = adata,
    n_neighbors_list = list(as.integer(6), as.integer(8)),
    n_clusters_range = as.integer(4:8)
)
cat("\n--------------------\n\n")
py_get_item(data$uns, "silhouette_scores")
```

The highest silhouette score is observed for `n_neighbors = 8` and
`n_clusters = 4`, suggesting this setting gives the clearest separation
among the tested options.

The setting `n_neighbors = 8` and `n_clusters = 7` is also close, so the
final choice can consider both the score and biological
interpretability. With `n_neighbors = 8` and `n_clusters = 7` selected,
the data is now ready for spatial region clustering.

## Divide the Tissue Into Spatial Regions

After choosing the clustering settings from previous function, we
finally divide the tissue into spatial regions.

This step uses spatially constrained hierarchical clustering to create
tissue regions. Cells are grouped based on their marker profiles, but
the clustering is restricted by the spatial neighborhood graph. This
means cells can be grouped together only when they are connected through
nearby cells in the tissue, helping the regions remain spatially
coherent.

The inputs are:

- `adata`: the AnnData object.
- `n_clusters = 7`: the number of spatial regions to create.
- `n_neighs = 8`: the number of nearby cells used to define local
  neighborhoods.

``` r

adata <- repspat$spatial_constrained_hac(
    adata = adata,
    n_clusters = as.integer(7),
    n_neighs = as.integer(8)
)
cat("\n--------------------\n\n")
py_get_item(data$obs, "labels")
```

The function returns the updated AnnData object, with the assigned
region number for every cell stored in `adata.obs["labels"]`. The data
is now ready to be visualized and grouped into blocks for comparison.

## Visualize the Spatial Regions

Now that each cell has an assigned region, we can plot the tissue
colored by region to check that the regions look spatially coherent. The
coordinates and region labels are stored in `adata`, so they can also be
pulled out directly to build your own custom visualizations or analyses.

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

This gives us a simple check of how the sample was divided before
continuing. The data is now ready for block based comparison, once
grouped into blocks in the next step.

## Visualize Features Across Regions

The region labels tell us how the tissue was partitioned, but not what
makes each region distinct. This function plots each region alongside
its most defining features: if the data has binary features, it shows
the features most frequently present in that region; otherwise, it shows
the numeric features most enriched in that region compared to the rest
of the tissue.

The inputs are:

- `adata`: the sample containing the region labels.
- `label_key = "labels"`: the column identifying each cell’s assigned
  region.
- `top_n = 10`: the number of top features to display per region.
- `figsize`: the size of each region’s plot.
- `point_size`: the size of the spatial scatter points.
- `alpha`: the transparency of the spatial scatter points.

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

The function returns a dictionary mapping each region’s cluster ID to
its `(figure, axes)` pair, which we loop through below to display every
region’s plot.

## Create Permutation Blocks

Cells within the same spatial region tend to look similar to their
neighbors, so they shouldn’t be treated as fully independent
observations during permutation testing. To account for this, we group
cells into smaller blocks within each region and permute at the block
level instead of the individual-cell level.

The inputs are:

- `adata`: the sample with region labels assigned in the previous step.
- `knn = 8`: the approximate number of cells per block, matching the
  neighborhood size used earlier. This keeps nearby cells grouped
  together while still splitting larger regions into several blocks.

``` r

adata <- repspat$create_blocks(
    adata = adata,
    knn = as.integer(8)
)

cat("\n--------------------\n\n")
py_get_item(data$obs, "repspat_block_id")
```

The function returns the updated AnnData object, with the block
assignment for every cell stored in `adata.obs["repspat_block_id"]` The
data is now ready for block-based comparison between spatial regions.

## Compare Region Pairs

The tissue has now been divided into regions and blocks, but we still
need to determine which regions have similar biological profiles.

- **Null hypothesis (H0):** the two regions come from the same
  underlying feature distribution. Any observed difference is due to
  chance, so the regions are spatially invariant and can be treated as a
  repeated spatial pattern.
- **Alternative hypothesis (H1):** the two regions come from different
  feature distributions, meaning they are biologically distinct.

To check whether two regions are similar, we need a way to compare their
feature profiles. A simple option would be to compare their means or
variances, but MMD does this more reliably: instead of reducing each
region to a single summary number, it compares the full shape of the two
distributions, so it can pick up differences in shape, spread, and
structure that a mean or variance comparison would miss.

MMD tells us how different two regions’ feature profiles are, as a
single number:

- Low MMD: the regions look similar.
- High MMD: the regions look different.

The inputs are:

- `adata`: the selected sample containing distances, regions, and
  blocks.
- `kernel`: the rule used to convert distances into similarities.
  `repSpat` supports `"IMQ"` and `"Gaussian"`. Both give higher
  similarity to cells with more similar feature profiles — `"IMQ"`
  decreases similarity gradually as distance increases, while
  `"Gaussian"` decreases it more quickly and emphasizes very close
  profiles. Here, we use `"IMQ"`.
- `nperm`: the number of block permutations used to build the
  chance-based reference distribution. We use 200; more permutations
  give a more stable estimate but take longer to run.
- `adj_p`: the method used to correct p-values when testing many region
  pairs. We use the Benjamini–Hochberg method (`"BH"`), which limits
  false discoveries across all the pairwise tests.

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

The results includes:

- `region_1` and `region_2`: the regions being compared.
- `obs_mmd_sq`: the observed MMD difference.
- `p_value`: the unadjusted permutation p-value.
- `adj_p`: the p-value after correcting for multiple comparisons (here,
  using the Benjamini–Hochberg false discovery rate method).

## Visualize Which Regions Are Not Significantly Different

Region pairs where we couldn’t reject the null hypothesis
(`adj_p >= 0.05`) are our repeated spatial patterns, regions that can’t
be clearly distinguished from each other. This function draws those
pairs as a network, connecting regions that are similar and labeling
each connection with its observed MMD.

The inputs are:

- `adata`: the sample containing the MMD comparison results
  (`adata.uns["repspat_mmd_results"]`).
- `plot = TRUE`: whether to draw the network graph in addition to
  returning the matrix.
- `alpha`: the significance threshold used to decide which regions are
  linked (region pairs with `adj_p >= alpha` are considered similar).

``` r

similarity_matrix <- py_to_r(repspat$pairwise_results_to_matrix(
    adata,
    plot = TRUE
))

similarity_matrix
```

A small adjusted p-value provides evidence that two regions differ. A
larger adjusted p-value means that the test did not find strong evidence
of a difference; it does not prove that the regions are identical.

The function returns a region-by-region adjacency matrix: regions
considered similar (repeated patterns) are connected with their MMD
value as the edge weight, while regions that were significantly
different are left unconnected (0).

## Merge Regions Into Repeated Spatial Patterns

The MMD comparison and network plot showed which region pairs could not
be statistically distinguished from each other, these are the repeated
spatial patterns. Here, we manually specify which region labels to merge
based on that result, combining them into a single shared label.

``` r

clusters_to_merge <- c(2, 3, 7)   # edit: region labels you want merged
new_label <- min(clusters_to_merge)
labels <- py_to_r(py_get_item(adata$obs, "labels"))
labels[labels %in% clusters_to_merge] <- new_label
py_set_item(adata$obs, "updated_labels", as.integer(labels))
```

The merged region assignment is now ready to be visualized alongside the
original regions.

## Visualize the Merged Regions

Plotting the tissue by `updated_labels` shows how the regions look once
statistically indistinguishable ones have been merged into shared
groups.

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

`repSpat` is a nonparametric framework for detecting repeated spatial
patterns in spatial data. It combines constrained hierarchical
clustering with a kernel based maximum mean discrepancy test and a block
permutation procedure. Together, these correct for the over segmentation
that constrained clustering methods often produce and allow spatially
separate regions with similar feature distributions to be identified as
the same underlying pattern.

Because the framework relies on distance based kernels, it generalizes
across a variety of attribute types and controls the false discovery
rate effectively even under weak to moderate spatial dependence. Applied
here to the TNBC sample, it successfully recovered repeated patterns
among the clustered regions, demonstrating its broader value for spatial
omics data, where the same biological neighborhood can recur at multiple
locations within a tissue.

One limitation worth noting is that the test shows limited sensitivity
to small shifts in spatial dependence alone, such as changes in
autocorrelation strength. Repeated patterns are therefore best
interpreted alongside region size, location, and biological context
rather than the test result in isolation.
