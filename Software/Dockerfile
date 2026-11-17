FROM bioconductor/bioconductor_docker:RELEASE_3_18

# Install system libraries required by tidyverse (ragg/textshaping/xml2) and Bioconductor packages
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff-dev \
    libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# Install required CRAN and Bioconductor packages using a Posit Package Manager
# binary snapshot aligned with Bioconductor 3.18 / R 4.3 (Ubuntu 22.04 Jammy)
RUN Rscript -e "options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/__linux__/jammy/2024-04-30')); BiocManager::install(c('psych', 'reshape2', 'tidyverse', 'optparse', 'DESeq2', 'Biostrings'), update = FALSE, ask = FALSE)"

# Verify all required R packages installed and load cleanly
RUN Rscript -e "pkgs <- c('psych', 'reshape2', 'tidyverse', 'optparse', 'DESeq2', 'Biostrings', 'MASS'); missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]; if (length(missing)) stop('Missing R packages: ', paste(missing, collapse = ', ')); invisible(lapply(pkgs, library, character.only = TRUE))"

# Create a directory for the scripts
RUN mkdir -p /opt/Y2H-SCORES

# Copy only the necessary R scripts to the container
COPY *.R /opt/Y2H-SCORES/

# Create a wrapper script so the container executes "run_scores.R" without the user typing it
RUN echo '#!/bin/bash' > /usr/local/bin/run_scores && \
    echo 'exec Rscript /opt/Y2H-SCORES/run_scores.R "$@"' >> /usr/local/bin/run_scores && \
    chmod +x /usr/local/bin/run_scores

# Default working directory for Docker (run_scores.R also resolves its own script directory
# so sourcing sibling scripts works from any working directory, including under Apptainer)
WORKDIR /opt/Y2H-SCORES

# Define the container entrypoint to automatically trigger the wrapper
ENTRYPOINT ["run_scores"]
