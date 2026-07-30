FROM bioconductor/bioconductor_docker:RELEASE_3_23-R-4.6.1

# install.packages() and BiocManager::install() only warn when a package fails to
# install, and `R -e` exits 0 either way, so each install is followed by a check
# that turns a missing package into a failed build rather than a broken image.
RUN R -q -e "pkgs <- c('remotes', 'reticulate', 'knitr', 'rmarkdown'); install.packages(pkgs); missing <- setdiff(pkgs, rownames(installed.packages())); if (length(missing)) stop('failed to install: ', paste(missing, collapse = ', '))"

RUN R -q -e "pkgs <- c('BiocStyle', 'SpatialExperiment', 'SingleCellExperiment', 'SummarizedExperiment'); BiocManager::install(pkgs, ask = FALSE, update = FALSE); missing <- setdiff(pkgs, rownames(installed.packages())); if (length(missing)) stop('failed to install: ', paste(missing, collapse = ', '))"

RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-venv \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m venv /opt/repspat-venv

ENV PATH="/opt/repspat-venv/bin:${PATH}"
ENV RETICULATE_PYTHON="/opt/repspat-venv/bin/python"

RUN python -m pip install --upgrade pip setuptools

RUN python -m pip install git+https://github.com/amanpreet60/repspat.git

COPY --chown=rstudio:rstudio . /home/rstudio