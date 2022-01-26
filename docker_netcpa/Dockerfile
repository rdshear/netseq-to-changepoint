# /bin/bash -c docker build . -t rdshear/seqcpa
#############################################################
# Changepoint analysis of NET-sdeq data sets
#############################################################

# TODO Add user and working directory

# TODO dclare version
FROM bioconductor/bioconductor_docker
LABEL description="Docker image with components for netseq-to-changepoints pipeline"

WORKDIR $HOME

RUN Rscript --vanilla -e 'BiocManager::install(c("plyranges", "Biostrings", \
    "GenomicRanges", "rtracklayer", \
    "breakpoint", "doMC", "foreach", "foreach", "kableExtra"), \
     quietly = TRUE, update = TRUE, ask = FALSE)'

ADD scripts/ /scripts/

CMD ["/bin/bash"]
