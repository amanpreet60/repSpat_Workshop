# Reference

| Function | Description |
|----|----|
| [`compute_distances()`](https://amanpreet60.github.io/repSpat_Workshop/articles/function_reference/compute_distances.md) | Compute pairwise feature distances from a selected AnnData layer. |
| [`spatial_silhouette_analysis()`](https://amanpreet60.github.io/repSpat_Workshop/articles/function_reference/spatial_silhouette_analysis.md) | Score combinations of neighborhood size and number of clusters using a spatial silhouette score. |
| [`spatial_constrained_hac()`](https://amanpreet60.github.io/repSpat_Workshop/articles/function_reference/spatial_constrained_hac.md) | Partition a tissue sample into spatially constrained clusters with constrained agglomerative hierarchical clustering. |
| [`create_blocks()`](https://amanpreet60.github.io/repSpat_Workshop/articles/function_reference/create_blocks.md) | Create permutation blocks within each cluster for block-based testing. |
| [`multiple_comparison()`](https://amanpreet60.github.io/repSpat_Workshop/articles/function_reference/multiple_comparison.md) | Pairwise MMD² comparison between all clusters, with multiple-testing correction. |
| [`plot_spatial_clusters()`](https://amanpreet60.github.io/repSpat_Workshop/articles/function_reference/plot_spatial_clusters.md) | Plot the tissue colored by assigned cluster labels. |
| [`plot_cluster_feature_presence()`](https://amanpreet60.github.io/repSpat_Workshop/articles/function_reference/plot_cluster_feature_presence.md) | Plot each cluster alongside its most defining features. |
| [`pairwise_results_to_matrix()`](https://amanpreet60.github.io/repSpat_Workshop/articles/function_reference/pairwise_results_to_matrix.md) | Build a cluster similarity matrix and network plot from MMD² results. |

All user-facing `repspat` functions, called from R via `reticulate` in
the workshop vignettes.

Source: [`amanpreet60/repspat`](https://github.com/amanpreet60/repspat)
