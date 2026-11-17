FROM bioconductor/bioconductor_docker:RELEASE_3_18

# Install required CRAN packages
RUN R -e "install.packages(c('psych', 'reshape2', 'tidyverse', 'optparse'), repos='https://cloud.r-project.org')"

# Install required Bioconductor packages
RUN R -e "BiocManager::install(c('DESeq2', 'Biostrings'))"

# Verify all required R packages installed and load cleanly
RUN Rscript -e "invisible(lapply(c('psych', 'reshape2', 'tidyverse', 'optparse', 'DESeq2', 'Biostrings', 'MASS'), library, character.only = TRUE))"

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
