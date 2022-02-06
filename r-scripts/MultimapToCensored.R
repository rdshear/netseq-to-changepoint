# MultimapToCensored.R
library(glue)
library(magrittr)
library(tidyverse)
library(Rsamtools)
library(GenomicAlignments)
library(rtracklayer)
library(corrplot)
library(fitdistrplus)

set.seed(20200101)
gene_list <- import("/Users/robertshear/Documents/n/groups/churchman/rds19/data/S005/genelist.gff", genome = "sacCer3")
names(gene_list) <- gene_list$ID
bam_directory <- "/n/groups/churchman/rds19/data/S005/mm-to-censor/"
f <- tibble(sample_id = paste0("SRR1284006", 6:9),
            bam_file = paste0(bam_directory, sample_id, ".bam"))
n_genes <- 20

OverlappedRanges <- function(q, s) s[subjectHits(findOverlaps(q, s))]

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

complementStrand <- function(u) c(`+`="-", `-`="+")[as.character(strand(u))]

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
  OverlappedRanges
  v1 <- a[subjectHits(GenomicRanges::findOverlaps(g, a))]
  v <- OverlappedRanges(g, a)
  w1 <- b[subjectHits(findOverlaps(g, b))]
  w <- OverlappedRanges(g, b)
  browser()
  
  x$score <- x$score + b[subjectHits(u)]$score
  y <- c(x, v, w)
  sort(y)
}

jaccard_similarity <- function(a,b) sum(width(GenomicRanges::intersect(a, b))) / 
  sum(width(GenomicRanges::union(a, b)))

ApplyMask <- function(s, m) {
  s <- scores$gsignal[1][[1]]
  m <- scores$gmask[1][[1]]
  
  msk <- gaps(GRanges(m)) %>% .[strand(.) != "*"]
  r1 <- s[subjectHits(findOverlaps(msk, s))]
  r <- OverlappedRanges(msk, s)
  browser()
  #TODO: add NA's and ZEROS
  
}

if (n_genes > 0) {
  gene_list <- sort(sample(gene_list, n_genes))
}

# This is a 'rough' filter. Gets rid of anti-sense reads and reads far away from target gene
# Just for efficiency. There is no provision for strand-aware filtering in readGAlignments
bam_read_mask <- gene_list
strand(bam_read_mask) <- "*"
bam_read_mask <- GenomicRanges::reduce(bam_read_mask)

  # The reads are from cDNA, therefore the 5'-end is the 3'-end of the nascent RNA
  # which is the last base exposed from the elongation complex.
  # THerefore, we will declare the occupancy to be the 5'-end. And we will
  # swap the strand information when we actually du=o the counts.
scores <- f  %>%
  mutate(bamreads = map(bam_file, function(u) {
    u <- GRanges(readGAlignments(u, param=
                ScanBamParam(tag = c("NH", "HI"),
                  what = c("qname", "cigar", "qwidth"), which = bam_read_mask)), 
              seqinfo = seqinfo(gene_list))
    glst <- gene_list
    strand(glst) <- complementStrand(glst)
    u <- OverlappedRanges(glst, u)
    u$cigstart <- explodeCigarOps(u$cigar) %>% 
      map(paste0, collapse="") %>% 
      unlist %>%
      stringi::stri_sub(., if_else(as.character(strand(u)) == "+", 1, -1), length = 1)
    u <- u[u$cigstart == "M"]
    strand(u) <- complementStrand(u)
    u <- split(u, u$HI)
    tibble(n_multi = as.integer(names(u)), gmask = as.list(u))
  })) %>%
  mutate(gsignal = map(bamreads, function(u) {
      GenomicRanges::resize(u$gmask[[1]], width = 1, fix = "end")  %>%
      sort
    })) %>%
  mutate(gsignal = map(gsignal, function(u) 
    SamToScore(OverlappedRanges(gene_list, u)))) %>% 
  mutate(gmask = map(bamreads, function(masktab) {
    cum <- GRanges(seqinfo = seqinfo(masktab$gmask[[1]]))
    result <- list()
    masktab <- masktab[-1, ]
    for (i in masktab$gmask) {
      cum <- OverlappedRanges(gene_list, c(cum, i))
      result <- append(result, list(cum))
    }
    masktab$gmask <- result
    masktab
  }))

# mask x sample
scores %>% 
  unnest(gmask) %>%
  dplyr::filter(n_multi == 4) -> w

jaccard_similarity(w$gmask[[1]], w$gmask[[4]])
# check similarity between mask and gsignal



# TODO Construct the appropriate censor mask
  
scores %>%
  unnest(gmask) %>% 
  mutate(jsim = map2_dbl(gsignal, gmask, jaccard_similarity),
        npos_signal = map_int(gsignal, length), 
        npos_mask =  map_int(gmask, length), 
        mean_signal = map_dbl(gsignal, function(u) mean(u$score))) %>%
  dplyr::select(!(bam_file:gsignal)) %>%
  dplyr::select(!gmask) -> z

z %>% ggplot(aes(group = n_multi, x = n_multi, y = jsim)) +
  geom_col() + facet_wrap(vars(sample_id))

scores %>%
  unnest(gmask)%>% filter(n_multi == 4) %>% dplyr::select(sample_id, gmask) -> signals
seq(nrow(signals)) %>% expand.grid(A=., B=.) %>% filter(A < B) %>%
  mutate(jsim = map2_dbl(A, B, function(a, b) 
    jaccard_similarity(signals$gmask[[a]], signals$gmask[[b]]))) %>%
  rbind(data.frame(A = 1:4, B = 1:4, jsim = 1.0)) -> zcor

pivot_wider(zcor, names_from = B, values_from = jsim)
# TODO: Exclude wt-4?
# TODO: look at mask jitter?
n <- nrow(signals)
z_mat <- matrix(ncol = n, nrow = n)
z_mat[as.matrix(as.matrix(zcor[,1:2]))] = zcor$jsim
rownames(z_mat) <- colnames(z_mat) <- signals$sample_id
corrplot(z_mat, is.corr = FALSE, diag = FALSE, type = "upper", method = "pie")
z_mat
