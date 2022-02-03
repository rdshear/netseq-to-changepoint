# MultimapToCensored.R
library(glue)
library(tidyverse)
library(Rsamtools)
library(GenomicAlignments)
library(rtracklayer)
library(corrplot)

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
  if (is.na(a)) {
    return(b)
  }
  if (is.na(b)) {
    return(a)
  }
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

gene_list <- import("/Users/robertshear/Documents/n/groups/churchman/rds19/data/S005/genelist.gff", genome = "sacCer3")

n_genes <- 10
gene_list <- sort(sample(gene_list, n_genes))
bam_directory <- "/n/groups/churchman/rds19/data/S005/mm-to-censor/"

f <- tibble(sample_id = paste0("SRR1284006", 6:9),
       bam_file = paste0(bam_directory, sample_id, ".bam"))

f <- f[1,] # TODO
scores <- f %>% 
  mutate(mask_list = map(bam_file, function(bf) {
    readGAlignments(bf, param=ScanBamParam(tag = c("NH", "HI"), which = gene_list)) %>%
      GRanges(.) %>%
      split(., .$HI) %>% as.list %>%
      map(SamToScore) -> r
      gmask <- r[-1]
      for (i in 2:length(gmask)) {
        gmask[[i]] <- AddScores(gmask[[i - 1]], gmask[[i]])
      }
      tibble(signal = list(r[[1]]), tibble(gmask, n_multi = as.integer(names(gmask))))
    }),
    bam_file = NULL) %>%
  unnest(mask_list)
  

# check similarity between mask and signal
scores %>% 
  mutate(bam_file = NULL, jsim = map2_dbl(signal, gmask, jaccard_similarity)) %>%
  select(sample_id, n_multi, jsim) ->z

z %>% pivot_wider(names_from = n_multi, values_from = jsim) -> zmat


z %>% ggplot(aes(group = n_multi, x = n_multi, y = jsim)) + geom_col() + facet_wrap(vars(sample_id))

scores %>% filter(n_multi == 4) %>% select(sample_id, signal) -> signals
seq(nrow(signals)) %>% expand.grid(A=., B=.) %>% filter(A < B) %>%
  mutate(jsim = map2_dbl(A, B, function(a, b) jaccard_similarity(signals$signal[[a]], signals$signal[[b]]))) %>% 
  rbind(data.frame(A = 1:4, B = 1:4, jsim = 1.0)) -> zcor

n <- nrow(signals)
z_mat <- matrix(ncol = n, nrow = n)
z_mat[as.matrix(as.matrix(zcor[,1:2]))] = zcor$jsim
rownames(z_mat) <- colnames(z_mat) <- signals$sample_id
corrplot(z_mat, is.corr = FALSE, diag = FALSE, type = "upper", method = "pie")
z_mat
