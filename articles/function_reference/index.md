# Reference

| Function | Description |
|----|----|
| [`compute_distances()`](https://amanpreet60.github.io/repSpat_Workshop/articles/function_reference/compute_distances.md) | Compute pairwise feature distances from a required AnnData layer. |
| [`spatial_silhouette_analysis()`](https://amanpreet60.github.io/repSpat_Workshop/articles/function_reference/spatial_silhouette_analysis.md) | Score combinations of neighborhood size and cluster count using a spatially-aware silhouette score. |
| [`spatial_constrained_hac()`](https://amanpreet60.github.io/repSpat_Workshop/articles/function_reference/spatial_constrained_hac.md) | Divide a tissue sample into spatial regions with spatially constrained hierarchical clustering. |
| [`create_blocks()`](https://amanpreet60.github.io/repSpat_Workshop/articles/function_reference/create_blocks.md) | Create permutation blocks within each region for block-based testing. |
| [`multiple_comparison()`](https://amanpreet60.github.io/repSpat_Workshop/articles/function_reference/multiple_comparison.md) | Pairwise MMD comparison between all regions, with multiple-testing correction. |
| [`plot_spatial_clusters()`](https://amanpreet60.github.io/repSpat_Workshop/articles/function_reference/plot_spatial_clusters.md) | Plot the tissue colored by assigned region labels. |
| [`plot_cluster_feature_presence()`](https://amanpreet60.github.io/repSpat_Workshop/articles/function_reference/plot_cluster_feature_presence.md) | Plot each region alongside its most defining features. |
| [`pairwise_results_to_matrix()`](https://amanpreet60.github.io/repSpat_Workshop/articles/function_reference/pairwise_results_to_matrix.md) | Build a region similarity matrix and network plot from MMD results. |

All user-facing `repspat` functions, called from R via `reticulate` in
the workshop vignettes.

Source: [`amanpreet60/repspat`](https://github.com/amanpreet60/repspat)
