# DiscoverBreakpointsScatter.R
# This step distributes the reference gene ranges to
# an array of files, each of which will have ranges for the 
# worker task to execute in gff3 format
# 
# Usage:
# Rscript --vanilla scripts/DiscoverBreakpointsScatter.R \
#   {reference.gene.list.filename} \
#   {occupancy.bedgraph.pos.gz} 
#   {occupancy.bedgraph.neg.gz}
#   {max.gene.length}     # regions of interest longer than this will only 
#                         # be searched to this lengh
#   {n.genes} \
#   {n.shards} 
#
# Example:
# Rscript --vanilla scripts/DiscoverBreakpointsScatter.R \
    # /n/groups/churchman/rds19/data/S005/genelist.gff \
    # 800
    # 12 \
    # 2

# To run with embedded parameters, set DEBUG.TEST <- TRUE

suppressPackageStartupMessages({
  library(rtracklayer)
  library(tidyverse)
  library(plyranges)
  library(glue)
})
set.seed(20190416)

# DEBUG ONLY FROM HERE.....
DEBUG.TEST <- TRUE
if (interactive() && exists("DEBUG.TEST")) {
  print("DEBUG IS ON -- COMMAND LINE PARAMETERS IGNORED")
  setwd("~/temp/")
  commandArgs <- function(trailingOnly) {
    c("/n/groups/churchman/rds19/data/S005/genelist.gff",
      "/n/groups/churchman/rds19/data/S005/wt-1.pos.bedgraph.gz",
      "/n/groups/churchman/rds19/data/S005/wt-1.neg.bedgraph.gz",
      "800", # Maximum gene body length
      "12", # n.genes
      "1") # n.shards
  }
}
#  .............TO HERE

args <- commandArgs(trailingOnly = TRUE)

subject_genes.filename <- args[1]
bedgraph.filename.pos <- args[2]
bedgraph.filename.neg <- args[3]
GeneMaxLength <- as.numeric(args[4]) # Truncate gene to this length (or inf if 0)
maxGenes <- as.numeric(args[5]) # if > 0 sample this number of genes
n.shards <- as.numeric(args[6]) # shards


sprintf("Starting at %s.  shards = %s. Max genes = %d", 
        Sys.time(), n.shards, maxGenes)

g <- import.gff3(subject_genes.filename, genome = "sacCer3", 
                  feature.type = "gene", colnames = "ID")

if (maxGenes > 0 & maxGenes < length(g)) {
  g <- g[sort(sample(length(g), maxGenes))]
}

if (GeneMaxLength > 0) {
  g <- resize(g, fix = "start", 
              ifelse(width(g) > GeneMaxLength, 
                     GeneMaxLength, width(g)))
}

map2(c(bedgraph.filename.pos, bedgraph.filename.neg), c('+', '-'), 
     function(filename, strand) {
       result <- import(filename, genome = "sacCer3")
       strand(result) <- strand
       result
     }) %>% 
  GRangesList() %>% # Can't unlist a list of GRanges (as of 9/2021)
  unlist() %>%
  join_overlap_inner_directed(g) %>%
  as_tibble() %>%
  nest_by(ID) %>%
  inner_join(g, copy = TRUE) %>%
  ungroup() %>%
  mutate(width = NULL, shard = row_number() %% n.shards) %>%
  group_by(shard) %>%
  group_map(function(u, v)  saveRDS(u, glue("shard_{v}.rds")))

print(sessionInfo())
