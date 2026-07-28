# repSpat

This workshop demonstrates how to use the Python implementation of repSpat within R through the reticulate package. RepSpat detects repeated spatial patterns in spatial omics data, defined as spatially separated tissue regions with similar distributions of features, such as gene expression, cell types, binary markers, and other molecular or cellular measurements.
The repSpat provides a nonparametric statistical inference framework based on spatially constrained clustering followed by multiple hypothesis testing between clusters using the maximum mean discrepancy squared (MMD²) statistic and a block-permutation procedure. This enables formal hypothesis testing for repeated spatial patterns. 
The Python implementation of repSpat accepts input data in the [AnnData](https://anndata.readthedocs.io/en/latest/index.html) (.h5ad) format, a widely used data structure for single-cell and spatial omics analyses. For Bioconductor users, SpatialExperiment objects can be readily converted to AnnData format within R before analysis.

The Python package is available on GitHub, and this workshop demonstrates how to use its functionality from R through the reticulate package.

```r
library(reticulate)
reticulate::py_install("git+https://github.com/amanpreet60/repspat.git")
repspat <- import("repspat", convert = FALSE)
```

## Workshop Goals

By the end of this workshop, participants will be able to:

- Understand the concept of repeated spatial patterns and their importance in spatial omics data.
- Apply constrained agglomerative hierarchical clustering (CAHC) to obtain partition tissue into spatially contiguous clusters.
- Test for spatial invariance between clusters using the MMD² statistic
- Approximate the null distribution of the MMD² statistic using block permutation that preserves spatial dependence.
- Reassign cluster labels based on pairwise similarity to identify repeated spatial patterns.
- Apply the complete repSpat workflow to analyze spatial omics datasets and interpret the resulting repeated spatial patterns.

## Workshop Length

90 minutes.

## Required R Packages

- BiocStyle (Bioconductor)
- SpatialExperiment (Bioconductor)
- SingleCellExperiment (Bioconductor)
- SummarizedExperiment (Bioconductor)
- reticulate (CRAN)
- knitr (CRAN) 
- rmarkdown (CRAN)

## Prerequisites

The following background is preferred but not required:

- Basic familiarity with R.
- Familiarity with spatial omics data eg. Visium, MIBI-TOF.
- Basic experience working with SpatialExperiment object.

## Docker

This workshop can be run using the provided Docker image.

Steps to run the workshop using Docker:

1. Open a terminal and navigate to the directory where you would like to run the workshop.

```terminal
cd /path/to/your/working/directory
```

2. Pull the workshop Docker image.

```text
docker pull ghcr.io/amanpreet60/repspat_workshop:f433744
```

3. Start the Docker container.

```terminal
docker run -e PASSWORD=bioc -p 8787:8787 ghcr.io/amanpreet60/repspat_workshop:f433744
```

4. Open your web browser and visit:

```text
http://localhost:8787
```

5. Log in to RStudio using

```text
username: rstudio
password: bioc
```

## Vignettes

The workshop includes the following example analyses:

MIBI-TOF spatial proteomics TNBC (Triple-Negative Breast Cancer)

```text
vignettes/tnbc_example.Rmd
```

MIBI-TOF spatial proteomics TNBC with binary markers

```text
vignettes/tnbc_example_binary.Rmd
```

10x Visium spatial transcriptomics of mouse brain

```text
vignettes/mouse_example.Rmd
```

## References

1. Senanayake, R. & Jeganathan, P. (2026). A Robust Nonparametric Framework for Detecting Repeated Spatial Patterns. Spatial Statistics, 101025. https://doi.org/10.1016/j.spasta.2026.101025

2. González-Almagro, G., Peralta, D., De Poorter, E., Cano, J.-R., & García, S. (2023). Semi-Supervised Constrained Clustering: An In-Depth Overview, Ranked Taxonomy and Future Research Directions. arXiv:2303.00522. https://arxiv.org/abs/2303.00522

3. Keren, L., Bosse, M., Marquez, D., et al. (2018). A Structured Tumor-Immune Microenvironment in Triple Negative Breast Cancer (TNBC) Revealed by Multiplexed Ion Beam Imaging. *Cell*, 174(6), 1373-1387. https://www.cell.com/fulltext/S0092-8674(18)31100-0

4. 10x Genomics. Visium spatial gene expression, mouse coronal brain section. Distributed via the Bioconductor `STexampleData` package. https://bioconductor.org/packages/STexampleData

5. Virshup, I., Rybakov, S., Theis, F. J., Angerer, P., & Wolf, F. A. (2024). anndata: Access and store annotated data matrices. *Journal of Open Source Software*, 9(101), 4371. https://doi.org/10.21105/joss.04371
