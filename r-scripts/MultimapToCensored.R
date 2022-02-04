# MultimapToCensored.R
library(glue)
library(magrittr)
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

SamToScore <- function(u) {
  split(u, strand(u))[1:2] %>% as.list %>%
     map2(names(.), function(u, s) {
       coverage(u) %>% 
        bindAsGRanges() %>%
        .[.$V1 > 0] %T>%
        {strand(.) <- s}
     }) %>%
    GRangesList %>% unlist %T>%
    sort.GenomicRanges %>%
    GPos(., score = rep(.$V1, width(.)))
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

gene_list <- import("/Users/robertshear/Documents/n/groups/churchman/rds19/data/S005/genelist.gff", genome = "sacCer3")
n_genes <- 10
gene_list <- sort(sample(gene_list, n_genes))
bam_directory <- "/n/groups/churchman/rds19/data/S005/mm-to-censor/"

f <- tibble(sample_id = paste0("SRR1284006", 6:9),
       bam_file = paste0(bam_directory, sample_id, ".bam"))

scores <- f %>%
  mutate(gmask = map(bam_file, function(u)
    readGAlignments(u, param=
                ScanBamParam(tag = c("NH", "HI"))))) %>%
  mutate(gmask = map(gmask, GRanges)) %>%
  mutate(gmask = map(gmask, GenomicRanges::resize, width = 1, fix = "start")) %>%
  mutate(gmask = map(gmask, function(u) 
        u[subjectHits(findOverlaps(gene_list, u))]), 
    gmask = map(gmask, function(u) as.list(split(u, u$HI))),
    gmask = map(gmask, function(u) map(u,SamToScore)),
    gsignal = map(gmask, function(u) u[[1]]),
    gmask = map(gmask, function(u) {
      Tot <<- GRanges()
      map(u[-1], function(w) Tot <<- AddScores(w, Tot))
    })
  )

# check similarity between mask and osignal
scores %>% 
  mutate(bam_file = NULL, jsim = map2_dbl(osignal, gmask, jaccard_similarity)) %>%
  select(sample_id, n_multi, jsim) ->z

z %>% pivot_wider(names_from = n_multi, values_from = jsim) -> zmat


z %>% ggplot(aes(group = n_multi, x = n_multi, y = jsim)) + geom_col() + facet_wrap(vars(sample_id))

scores %>% filter(n_multi == 4) %>% select(sample_id, osignal) -> signals
seq(nrow(signals)) %>% expand.grid(A=., B=.) %>% filter(A < B) %>%
  mutate(jsim = map2_dbl(A, B, function(a, b) jaccard_similarity(signals$osignal[[a]], signals$osignal[[b]]))) %>% 
  rbind(data.frame(A = 1:4, B = 1:4, jsim = 1.0)) -> zcor

n <- nrow(signals)
z_mat <- matrix(ncol = n, nrow = n)
z_mat[as.matrix(as.matrix(zcor[,1:2]))] = zcor$jsim
rownames(z_mat) <- colnames(z_mat) <- signals$sample_id
corrplot(z_mat, is.corr = FALSE, diag = FALSE, type = "upper", method = "pie")
z_mat
