# =============================================================================
#  ChAMP2 container
#
#  docker build -f docker/Dockerfile --build-arg R_VERSION=4.6.1 -t champ2:4.6.1 .
#  docker build -f docker/Dockerfile --build-arg R_VERSION=4.2.3 -t champ2:4.2.3 .
#
#  The build installs every required and optional dependency with the bundled
#  installer, then ChAMP itself, then runs the test suite.
# =============================================================================
ARG R_VERSION=4.6.1
FROM rocker/r-ver:${R_VERSION}

LABEL maintainer="champ450k@gmail.com" \
      org.opencontainers.image.title="ChAMP2" \
      org.opencontainers.image.description="Chip Analysis Methylation Pipeline for Illumina 450K/EPIC/EPICv2"

ENV DEBIAN_FRONTEND=noninteractive
ENV R_LIBS_USER=/usr/local/lib/R/site-library

RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev libssl-dev libxml2-dev zlib1g-dev \
        libpng-dev libjpeg-dev libtiff-dev libbz2-dev liblzma-dev \
        libglpk-dev libgsl-dev pandoc \
    && rm -rf /var/lib/apt/lists/*

COPY inst/scripts/install_deps.R /opt/champ2/install_deps.R
RUN Rscript /opt/champ2/install_deps.R --quiet || \
    echo "Some optional packages failed to install; ChAMP still works."

COPY . /opt/champ2/pkg
RUN R CMD INSTALL --no-multiarch /opt/champ2/pkg

RUN Rscript -e "library(ChAMP); champ_check_env(); testthat::test_dir('/opt/champ2/pkg/tests/testthat', package = 'ChAMP', stop_on_failure = TRUE)"

WORKDIR /data
CMD ["R"]
