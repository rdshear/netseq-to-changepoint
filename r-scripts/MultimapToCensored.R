# MultimapToCensored.R
library(glue)
library(tidyverse)
library(Rsamtools)
library(GenomicAlignments)
library(rtracklayer)

GRangesToGPos <- function(z) GPos(z, score = rep(z$score, width(z)))

BedgraphToGranges <- function(path, sample) {
  x <- mapply(function(strand_sym, infilename) {
    x <- import(infilename, genome = "sacCer3")
    x <- x[seqnames(x) != "chrM"]
    strand(x) <- strand_sym
    x
  }, list('+', '-'),  paste0(path,sample, ".", c("neg","pos"),".bedgraph.gz"), SIMPLIFY = FALSE)
  x <- GRangesToGPos(sort(c(x[[1]], x[[2]])))
}

SamToScore <- function(x) {
    x %>%
    resize(., 1, fix="start") %>%
    split(., strand(.)) %>%
    as.list(.) %>% .[1:2] %>%
    map2(names(.), function(u, s) {
      coverage(u) %>% 
        bindAsGRanges() %>% 
        .[.$V1 > 0] %>%
        {strand(.) <- s; .}
    }) %>% 
    GRangesList %>% 
    unlist %>% 
    {colnames(mcols(.)) <- "score";.} %>%
    sort.GenomicRanges %>%
    (function(x) {
      y <- GPos(x)
      y$score <- rep(x$score, width(x))
      y
    })
}

AddScores <- function(a, b) {
  u <- findOverlaps(a,b)
  x <- a[queryHits(u)]
  g <- gaps(GRanges(x))
  g <- g[strand(g) != "*"]
  v <- a[subjectHits(GenomicRanges::findOverlaps(g, a))]
  w <- b[subjectHits(findOverlaps(g, b))]
  
  x$score <- x$score + b[subjectHits(u)]$score
  y <- c(x, v, w)
  sort(y)
}

jaccard_similarity <- function(a,b) sum(width(GenomicRanges::intersect(a, b))) / sum(width(GenomicRanges::union(a, b)))

f <- "/n/groups/churchman/rds19/data/S005/mm-to-censor/SRR12840066.bam"

reads <- f %>%
  readGAlignments(f, param=ScanBamParam(tag = c("NH", "HI"), which = GRanges("chrI:1-200000"))) %>%
  GRanges(.) %>%
  split(., .$HI) %>% as.list %>%
  map(SamToScore) 

unique_scores <- reads[[1]]
dup_list <- reads[2]
reads <- reads[-1:-2]
for (r in reads) {
  dup_list <- base::append(dup_list, AddScores(dup_list[[length(dup_list)]], r))
}
names(dup_list) <- paste0("mask_", seq(length(dup_list)))
  
gene_list <- import("/Users/robertshear/Documents/n/groups/churchman/rds19/data/S005/genelist.gff")

dup_scores <- dup_list[[length(dup_list)]]
# For entire genome
jsim_all <- jaccard_similarity(unique_scores, dup_scores)

print(glue("Total unique scores = {sum(unique_scores$score)}"))
print(glue("Total dups scores = {sum(dup_scores$score)}"))

# For regions of interest 
z <- lapply(list(u = unique_scores, d = dup_scores), function(x) {
  ov <- findOverlaps(gene_list, x)
  x[subjectHits(ov)]
  })

jsim_row <- jaccard_similarity(z[[1]], z[[2]])

print(glue("Jaccard sim (genome) = {round(jsim_all,5)}"))
print(glue("Jaccard sim (row) = {round(jsim_row,5)}"))

ov <- findOverlaps(z[[1]], z[[2]])
diffs <- z[[1]][queryHits(ov)]$score -z[[2]][subjectHits(ov)]$score
u2 <- z[[1]][queryHits(ov)]
u2$dups_score <- diffs
# Investigate wide overlaps between unique and dup
ovx <- GenomicRanges::reduce(GRanges(u2), min.gapwidth = 3)
ovx_review <- ovx[width(ovx) > 4]
ovx_review





