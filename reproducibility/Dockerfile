# Dockerfile for exact reproduction of the v2.1 manuscript analyses.
# Base: rocker/r-ver:4.5.1 (matches the local R 4.5.1 snapshot, 2025-06-13).
#
# Build context = repository root:
#   docker build -f reproducibility/Dockerfile -t nhanes-opioid-disparities .
#
# Run with your local NHANES + MEPS data directory mounted:
#   docker run --rm -it \
#     -e STUDY_ROOT=/work \
#     -v $(pwd):/work \
#     nhanes-opioid-disparities

FROM rocker/r-ver:4.5.1

# System dependencies for arrow, survey, haven, and downstream packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    libxml2-dev libcurl4-openssl-dev libssl-dev \
    libsodium-dev libicu-dev libgit2-dev \
    pandoc git wget \
    && rm -rf /var/lib/apt/lists/*

# Install remotes first, then pinned packages via requirements.R
RUN R -e "install.packages('remotes', repos='https://cloud.r-project.org')"

WORKDIR /work
COPY reproducibility/requirements.R /tmp/requirements.R
RUN Rscript /tmp/requirements.R

# Default entrypoint: drop into R at the work directory
CMD ["R", "--no-save"]
