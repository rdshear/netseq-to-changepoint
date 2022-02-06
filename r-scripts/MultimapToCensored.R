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
bam_directory <- "/n/groups/churchman/rds19/data/S005/mm-to-censor/"
f <- tibble(sample_id = paste0("SRR1284006", 6:9),
            bam_file = paste0(bam_directory, sample_id, ".bam"))
n_genes <- 100

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
  v <- a[subjectHits(GenomicRanges::findOverlaps(g, a))]
  w <- b[subjectHits(findOverlaps(g, b))]
  
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
  r <- s[subjectHits(findOverlaps(msk, s))]
  
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

reads <- f %>%
  # The reads are from cDNA, therefore the 5'-end is the 3'-end of the nascent RNA
  # which is the last base exposed from the elongation complex.
  # THerefore, we will declare the occupancy to be the 5'-end. And we will
  # swap the strand information when we actually du=o the counts.
  mutate(bamreads = map(bam_file, function(u) {
    x <- GRanges(readGAlignments(u, param=
                ScanBamParam(tag = c("NH", "HI"),
                  what = c("qname", "cigar", "qwidth"), which = bam_read_mask)), 
              seqinfo = seqinfo(gene_list))
    glst <- gene_list
    strand(glst) <- complementStrand(glst)
    x <- x[subjectHits(findOverlaps(glst, x))]
    x$cigstart <- explodeCigarOps(x$cigar) %>% 
      map(paste0, collapse="") %>% 
      unlist %>%  
      stringi::stri_sub(., if_else(as.character(strand(x)) == "+", 1, -1), length = 1)
    x
    }))


# convert to occupancy table
scores <- reads %>%
  mutate(bamreads = map(bamreads, 
                     GenomicRanges::resize, width = 1, fix = "start")) %>%
  mutate(bamreads = map(bamreads, function(u) u[u$cigstart == "M"])) %>%
  mutate(bamreads = map(bamreads, function(u) {
    strand(u) <- complementStrand(u)
    u
  })) %>% 
  mutate(bamreads = map(bamreads, function(u) {
    u[subjectHits(findOverlaps(gene_list, u))]
  })) %>%
  mutate(gmask = map(bamreads, function(u) as.list(split(u, u$HI))),
    gmask = map(gmask, function(u) map(u,SamToScore))) %>%
  mutate(gsignal = map(gmask, function(u) u[[1]])) %>%
  # TODO Outer join on depth (HI)
  mutate(gmask = map(gmask, function(u) {
      Tot <- GRanges()
      r <- list()
      for (w in u[-1])  {
        Tot <- AddScores(w, Tot)
        r <- append(r, list(Tot))
      }
      tibble(n_multi = names(u[-1]),gmask = r)
    })
  ) %>%

# check similarity between mask and gsignal

# TODO Construct the appropriate censor mask

  mutate(roi_reads = map(bamreads, function(u)
         (function(v) split(u[subjectHits(v)], 
                  gene_list$ID[queryHits(v)]))(findOverlaps(gene_list, u)))) %>%
  mutate(roi_score = map(gsignal, function(u)
    (function(v) split(SamToScore(u[subjectHits(v)]), 
                       gene_list$ID[queryHits(v)]))(findOverlaps(gene_list, u))))
  
scores %>%
  unnest(gmask) %>% 
  mutate(jsim = map2_dbl(gsignal, gmask, jaccard_similarity),
        npos_signal = map_int(gsignal, length), 
        npos_mask =  map_int(gmask, length), 
        mean_signal = map_dbl(gsignal, function(u) mean(u$score)), 
        mean_mask = map_dbl(gmask, function(u) mean(u$score))) %>%
  dplyr::select(!(bam_file:bamreads)) %>% dplyr::select(!(gmask:gsignal)) -> z

z %>% pivot_wider(names_from = n_multi, values_from = c(jsim)) -> zmat


z %>% ggplot(aes(group = n_multi, x = n_multi, y = jsim)) +
  geom_col() + facet_wrap(vars(sample_id))

scores %>%
  unnest(gmask)%>% filter(n_multi == 4) %>% dplyr::select(sample_id, gsignal) -> signals
seq(nrow(signals)) %>% expand.grid(A=., B=.) %>% filter(A < B) %>%
  mutate(jsim = map2_dbl(A, B, function(a, b) 
    jaccard_similarity(signals$gsignal[[a]], signals$gsignal[[b]]))) %>% 
  rbind(data.frame(A = 1:4, B = 1:4, jsim = 1.0)) -> zcor

n <- nrow(signals)
z_mat <- matrix(ncol = n, nrow = n)
z_mat[as.matrix(as.matrix(zcor[,1:2]))] = zcor$jsim
rownames(z_mat) <- colnames(z_mat) <- signals$sample_id
corrplot(z_mat, is.corr = FALSE, diag = FALSE, type = "upper", method = "pie")
z_mat
