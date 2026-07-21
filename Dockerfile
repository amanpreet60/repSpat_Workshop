FROM bioconductor/bioconductor_docker:devel

RUN R -e "install.packages(c('remotes', 'reticulate', 'knitr', 'rmarkdown'))"

RUN R -e "BiocManager::install(c( \
    'BiocStyle', \
    'SpatialExperiment', \
    'SingleCellExperiment', \
    'SummarizedExperiment' \
))"

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

COPY --chown=rstudio:rstudio . /home/rstudio/repSpat_Workshop