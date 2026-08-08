# TNBC Example: Detecting Repeated Spatial Patterns

## Introduction

MIBI-TOF (Multiplexed Ion Beam Imaging by Time-of-Flight) is a spatial
proteomics technology that measures the expression of dozens of protein
markers while preserving their spatial locations within tissue. In
tumors such as triple-negative breast cancer (TNBC), the tumor
microenvironment consists of diverse cell populations that interact to
form distinct spatial organizations. Similar cellular neighborhoods or
microenvironments may occur in multiple, spatially separated regions of
the same tissue, giving rise to repeated spatial patterns. In this
vignette, we use repSpat to identify and compare repeated spatial
patterns in a TNBC MIBI-TOF dataset.

The workflow uses five main package components:

- `compute_distances()` computes pairwise feature distances from a
  selected AnnData layer.
- `spatial_silhouette_analysis()` evaluates combinations of the number
  of nearest neighbors and the number of clusters to guide the selection
  of constrained clustering parameters.
- `spatial_constrained_hac()` performs constrained agglomerative
  hierarchical clustering using the selected parameters.
- `create_blocks()` partitions each cluster into smaller blocks for
  block-permutation testing
- `multiple_comparison()` performs pairwise hypothesis tests between
  clusters to identify repeated spatial patterns using the maximum mean
  discrepancy squared (MMD²) statistic and a block-permutation procedure

At the end of the workflow, clusters that are not found to differ
significantly are grouped together as repeated spatial patterns.

In this vignette, we introduce each function in turn, describe its
inputs and outputs, demonstrate its use, and examine the resulting
output before proceeding to the next step.

## Load libraries

We first load the R packages required for the workshop. These packages
provide tools for working with spatial omics data and for calling the
Python implementation of repSpat from within R using the reticulate
package.

``` r

library(BiocStyle)
library(SpatialExperiment)
library(reticulate)
```

Next, we import the Python packages used in the analysis: anndata for
working with AnnData objects and repspat for detecting repeated spatial
patterns.

``` r

anndata <- import("anndata", convert = FALSE)
repspat <- import("repspat", convert = FALSE)
```

## Prepare the Sample

Before running `repSpat`, we preprocess the TNBC MIBI-TOF dataset. The
original dataset is stored as an R SpatialExperiment object (.rds)
containing multiple tissue samples. Since each sample represents an
independent tissue section with its own spatial layout, the analysis is
performed one sample at a time.

In this section, we:

- load the SpatialExperiment object from the `.rds` file,
- extract the marker measurements, cell information, and spatial
  coordinates,
- convert the selected sample into an AnnData object,
- save the AnnData object as an .h5ad file for use with the Python
  implementation of repSpat.

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

The resulting AnnData object contains the selected tissue sample
(`Sample_04`) and serves as the input to the repSpat workflow.

## Compute Feature Distances

The first step of the repSpat workflow computes pairwise feature
distances between cells based on their marker profiles. These distances
quantify the similarity between cells and are used in subsequent steps
to construct spatially constrained clusters.

The function takes the following inputs:

- `adata`: the AnnData object containing `Sample_04`.
- `layer = "layer 1"`: the feature matrix used to compute distances.
- `metric = "euclidean"`: the distance metric used to compare cells.

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

After computing the distances, the resulting distance matrix is stored
in adata\$obsp\[“repspat_distances”\] for use in the clustering step.
The feature matrix (layers\[“layer 1”\]), spatial coordinates
(obsm\[“spatial”\]), and cell information (obs) remain unchanged.

We can verify that the data have been loaded correctly and inspect the
main components of the AnnData object before proceeding.

``` r

cat("Feature matrix: first 5 cells\n")
print(py_to_r(py_get_item(adata$layers, "layer 1"))[1:5, 1:11])

cat("\nSpatial coordinates: first 5 cells\n")
print(head(py_to_r(py_get_item(adata$obsm, "spatial")), 5))

cat("\nDistance matrix: first 5 cells x first 5 cells\n")
print(py_to_r(py_get_item(adata$obsp, "repspat_distances"))[1:5, 1:5])
```

This confirms that the data have been prepared successfully and that the
pairwise feature distances have been computed. The resulting distance
matrix will be used to construct spatially constrained clusters in the
next step.

## Select Clustering Parameters

Before partitioning the tissue into spatially constrained clusters, we
first select appropriate clustering parameters. Specifically, we choose
the number of nearest neighbors used to define local spatial
connectivity and the number of clusters used to partition the tissue.

To guide this selection, the workflow evaluates combinations of these
parameters using a spatial silhouette score. Unlike the classical
silhouette score, the spatial silhouette score compares each cell with
its own cluster and only neighboring clusters, providing a measure of
cluster separation that accounts for the spatial organization of the
tissue. Higher scores indicate better-separated spatially constrained
clusters.

The candidate values can be adjusted according to the characteristics of
the dataset and informed by domain knowledge regarding the tissue
architecture.

The function takes the following inputs:

- `adata`: the AnnData object containing the precomputed pairwise
  feature distances.
- `n_neighbors_list = [6, 8]`: candidate numbers of nearest neighbors
  used to define the spatial connectivity graph.
- `n_clusters_range = 4:8`: candidate numbers of clusters.

``` r

adata <- repspat$spatial_silhouette_analysis(
    adata = adata,
    n_neighbors_list = list(as.integer(6), as.integer(8)),
    n_clusters_range = as.integer(4:8)
)
cat("\n--------------------\n\n")
py_get_item(adata$uns, "silhouette_scores")
```

The highest spatial silhouette score is obtained for n_neighbors = 8 and
n_clusters = 4, indicating the best cluster separation among the
parameter combinations evaluated. The combination n_neighbors = 8 and
n_clusters = 7 also achieves a comparable score.

Therefore, the final choice should consider both the spatial silhouette
score and the biological interpretability of the resulting clusters. In
this workshop, we proceed with n_neighbors = 8 and n_clusters = 7 for
the subsequent clustering analysis.

## Construct Spatially Constrained Clusters

After selecting the clustering parameters, we partition the tissue into
spatially constrained clusters using constrained agglomerative
hierarchical clustering (CAHC). Cells are grouped according to the
similarity of their feature profiles while enforcing spatial
connectivity through the neighborhood graph. As a result, cells are
merged only if they are connected through nearby cells, producing
spatially contiguous clusters.

The function takes the following inputs:

- `adata`: the AnnData object.
- `n_clusters = 7`: the number of clusters to construct.
- `n_neighs = 8`: the number of nearest neighbors used to define the
  spatial connectivity graph.

``` r

adata <- repspat$spatial_constrained_hac(
    adata = adata,
    n_clusters = as.integer(7),
    n_neighs = as.integer(8)
)
cat("\n--------------------\n\n")
py_get_item(adata$obs, "labels")
```

The function returns an updated AnnData object with the cluster
assignment for each cell stored in adata.obs\[“labels”\]. These
spatially constrained clusters are used in the subsequent steps for
block construction and hypothesis testing.

## Visualize the Spatially Constrained Clusters

After assigning each cell to a spatially constrained cluster, we
visualize the tissue using the cluster labels. This provides a simple
way to assess whether the resulting clusters are spatially contiguous
and consistent with the tissue structure.

The spatial coordinates and cluster labels are stored in the AnnData
object and can also be extracted to create customized visualizations or
conduct additional analyses.

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

The resulting plot shows the spatial distribution of the seven clusters
across the tissue. After inspecting the clustering result, we proceed to
partition each cluster into blocks for block-permutation testing.

## Visualize Features Across Clusters

The spatially constrained clusters show how the tissue has been
partitioned, but they do not explain what distinguishes each cluster.
This function visualizes the features that are most characteristic of
each cluster. For binary features, it displays the features with the
highest prevalence within each cluster. For continuous features, it
displays those with the largest standardized mean difference between the
cluster and all remaining cells.

The function takes the following inputs:

- `adata`: the AnnData object containing the cluster labels.
- `label_key = "labels"`: the column containing the cluster assignment
  for each cell.
- `top_n = 10`: the number of top features to display for each cluster.
- `figsize`: the size of each cluster in the plot
- `point_size`: the size of the points in the plot.
- `alpha`: the transparency of the plotted points.

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

The function returns a dictionary mapping each cluster to its
corresponding (figure, axes) pair. The code below iterates through these
figures and displays the feature plots for all clusters.

## Create blocks for permutation

Permutation tests assume that observations are exchangeable under the
null hypothesis. However, cells within a spatially constrained cluster
often exhibit local spatial dependence and therefore cannot be treated
as independent observations. To account for this dependence, repSpat
partitions each cluster into smaller blocks and performs permutations at
the block level rather than the individual-cell level.

The function takes the following inputs:

- `adata`: the AnnData object containing the cluster assignments.
- `knn = 8`: the approximate number of cells per block, matching the
  neighborhood size used during clustering. This keeps nearby cells
  grouped together while allowing larger clusters to be partitioned into
  multiple blocks.

``` r

adata <- repspat$create_blocks(
    adata = adata,
    knn = as.integer(8)
)

cat("\n--------------------\n\n")
py_get_item(adata$obs, "repspat_block_id")
```

The function returns an updated AnnData object with the block assignment
for each cell stored in adata.obs\[“repspat_block_id”\]. These block
assignments are then used in the block-permutation procedure to perform
hypothesis tests between clusters.

## Test Spatial Invariance Between Clusters

The tissue has now been partitioned into spatially constrained clusters,
and each cluster has been divided into permutation blocks. The next step
is to compare the feature distributions of all cluster pairs and
determine which pairs show evidence of difference in distribution

- **Null hypothesis (H0):** the two clusters have the same underlying
  feature distribution.
- **Alternative hypothesis (H1):** the two clusters have different
  feature distributions.

Failure to reject the null hypothesis indicates that the analysis does
not find sufficient evidence of a difference between the two cluster
distributions. Such cluster pairs are treated as candidate repeated
spatial patterns. In contrast, rejection of the null hypothesis provides
evidence that the clusters have different feature distributions.

The observed MMD² value measures the discrepancy between the feature
distributions of two clusters:

- A smaller MMD² value indicates greater similarity between the
  distributions.
- A larger MMD² value indicates greater dissimilarity between the
  distributions.

The function takes the following inputs:

- `adata`: the AnnData object containing the feature data, cluster
  assignments, and permutation-block assignments.
- `kernel`: the kernel used to compare feature distributions. repSpat
  supports “IMQ” and “Gaussian” kernels. Both assign greater similarity
  to cells with similar feature profiles. The inverse multiquadric
  (“IMQ”) kernel decreases more gradually as feature distance increases,
  whereas the Gaussian kernel decreases more rapidly. In this example,
  we use “IMQ”.
- `nperm`: the number of block permutations used to approximate the null
  distribution of the MMD² statistic.
- `adj_p = "BH"`: the method used to adjust p-values for the multiple
  pairwise tests. The Benjamini–Hochberg (“BH”) procedure is used here
  to control the false discovery rate.

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

The results are stored in adata\$uns\[“repspat_mmd_results”\] and
include:

- `region_1` and `region_2`: the pair of clusters being compared.
- `obs_mmd_sq`: the observed MMD² statistic.
- `p_value`: the unadjusted p-value.
- `adj_p`: the p-value adjusted for multiple comparisons using the
  selected procedure.

In the future versions, false discovery rate threshold (`alpha`) for
`multiple_comparison`and kernel parameter (`lambda`) can be changed by
users.

## Visualize Repeated Spatial Patterns

The pairwise hypothesis tests identify cluster pairs that do not show
evidence of different feature distributions. Using the
Benjamini–Hochberg procedure with a false discovery rate (FDR) threshold
of 0.05, cluster pairs with adj_p ≥ 0.05 are treated as candidate
repeated spatial patterns.

This function visualizes these relationships as a network. Each node
represents a cluster, and an edge connects two clusters that are not
significantly different. The edge is labeled with the corresponding
observed MMD² statistic.

The function takes the following inputs:

- `adata`: the AnnData object containing the pairwise comparison
  results.
- `plot = TRUE`: whether to display the network in addition to returning
  the similarity matrix.

``` r

similarity_matrix <- py_to_r(repspat$pairwise_results_to_matrix(
    adata,
    plot = TRUE
))

similarity_matrix
```

An adjusted p-value smaller than the FDR threshold (adj_p \< 0.05)
provides evidence that the two clusters have different feature
distributions. An adjusted p-value greater than or equal to the
threshold (adj_p ≥ 0.05) indicates that the analysis did not find
sufficient evidence to distinguish their feature distributions. These
cluster pairs are therefore treated as candidate repeated spatial
patterns.

The function returns a cluster-by-cluster adjacency matrix. Clusters
identified as repeated spatial patterns are connected by edges weighted
by their observed MMD² statistic, whereas significantly different
cluster pairs are represented by 0.

## Relabel Repeated Spatial Patterns

The pairwise hypothesis tests identify repeated spatial patterns. Based
on these results, we can assign a common label to clusters belonging to
the same repeated spatial pattern. This step changes only the cluster
labels.

``` r

clusters_to_merge <- c(2, 3, 7)   # edit: select labels to reassign
new_label <- min(clusters_to_merge)
labels <- py_to_r(py_get_item(adata$obs, "labels"))
labels[labels %in% clusters_to_merge] <- new_label
py_set_item(adata$obs, "updated_labels", as.integer(labels))
```

The updated labels identify clusters belonging to the same repeated
spatial pattern and can be visualized alongside the original cluster
assignments.

## Visualize the Relabeled Repeated Spatial Patterns

Plotting the tissue using updated_labels visualizes the repeated spatial
patterns identified by the hypothesis tests.

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
patterns in spatial omics data. The framework consists of four main
steps: constrained agglomerative hierarchical clustering, pairwise
hypothesis testing using the MMD² statistic, block permutation to
account for spatial dependence, and cluster relabeling based on the
hypothesis testing results.

By comparing the full feature distributions of spatially constrained
clusters, repSpat identifies spatially separated tissue regions that
exhibit similar feature distributions and therefore represent repeated
spatial patterns. Because the framework is based on distance-based
kernels, it can be applied to a wide range of feature types, including
continuous, binary, categorical, and compositional measurements.

In this vignette, we applied repSpat to a triple-negative breast cancer
(TNBC) MIBI-TOF dataset and identified repeated spatial patterns across
spatially separated regions of the tissue. These results illustrate how
repSpat can help characterize recurring tissue organization in spatial
omics data.

An important direction for future development is feature selection. The
current framework compares clusters using all available features,
whereas many applications would benefit from identifying the subset of
features that most strongly characterize repeated spatial patterns. We
are actively developing feature selection methods for integration into
repSpat, enabling more interpretable repeated spatial pattern detection
and improved biological insight.
