# MultimapToCensored.R
library(glue)
library(magrittr)
library(tidyverse)
library(Rsamtools)
library(GenomicAlignments)
library(rtracklayer)
library(corrplot)
library(fitdistrplus)

print(glue("Start time = {Sys.time()}"))
start_time <- proc.time()
set.seed(20200101)
gene_list <- import("/Users/robertshear/Documents/n/groups/churchman/rds19/data/S005/genelist.gff", genome = "sacCer3")
n_genes <- 20

names(gene_list) <- gene_list$ID 

z <- disjoin(gene_list, with.revmap = TRUE, ignore.strand = FALSE)$revmap
w <- unique(unlist(z[which(sapply(z, function(u) length(u) > 1))]))
if (length(w) > 0) {
  gene_list <- gene_list[-w]
}

if (n_genes > 0 && n_genes < length(gene_list)) {
  gene_list <- sample(gene_list, n_genes)
  }

gene_list <- GenomicRanges::sort(gene_list)

bam_directory <- "/n/groups/churchman/rds19/data/S006/"
f <- tibble(sample_id = paste0("SRR1284006", 6:9),
            bam_file = paste0(bam_directory, sample_id, ".bam"))
n_genes <- 0

IncludedRanges <- function(q, s) {
  x <- c(q, s)
  y <- disjoin(x, with.revmap = TRUE, ignore.strand = FALSE)
  s1 <- y[unique(which((sapply(as.vector(y$revmap), max) > length(q))))]
  OverlappedRanges(q, s1)
}

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
  #TODO: add NA's and ZEROS
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
      cum <- GenomicRanges::reduce(OverlappedRanges(gene_list, c(cum, i)))
      result <- append(result, list(cum))
    }
    masktab$gmask <- result
    masktab
  })) %>%
  mutate(topmask = map(gmask, function(u) u$gmask[which.max(u$n_multi)][[1]]))  %>%
  mutate(g_scores = map2(gsignal, topmask, function(s, m){
    masked_scores <- OverlappedRanges(gaps(m) %>% 
        .[as.character(strand(.)) != "*"], s)
    ov <- findOverlaps(gene_list, masked_scores)
    result <- split(masked_scores[subjectHits(ov)], names(gene_list)[queryHits(ov)])
  }))
  # calculate mask losses per gene


export(scores$g_scores[[1]] %>% unlist, "~/Downloads/tmp.bedgraph", format = "bedGraph", index = TRUE)
export(gene_list, "~/Downloads/genes_20.gff3", format = "gff3", index = TRUE)
export(scores$topmask[[1]], "~/Downloads/tmp_mask.gff3", format = "gff3", index = TRUE)

# check similarity between mask and gsignal

# TODO parameterize save location
#saveRDS(scores, "~/Downloads/MultimapToCensored.rds")

elapsed <-  proc.time() - start_time
print("elapsed time")
print(elapsed)
print(glue("time / gene  (n = {n_genes})"))
print(elapsed / n_genes)

scores %>%
  mutate(mask_losses = map(topmask, function(u) {
  u <-  IncludedRanges(gene_list, u)
  ov <-  findOverlaps(gene_list, u)
  s <- split(u[subjectHits(ov)], names(gene_list[queryHits(ov)]))
  u <- rbind(tibble(gene_id = names(s), mask_width = sum(width(s)), mask_segments = sapply(s, length), gene_width = width(gene_list[names(s)])),
             setdiff(names(gene_list), names(s)) %>% tibble(gene_id = ., mask_width = 0, mask_segments = 0,  gene_width = width(gene_list[.])))
  u$base_loss = u$mask_width / u$gene_width
  u
  })
) %>%
  transmute(sample_id, mask_losses) %>%
  unnest(mask_losses) %>%
  arrange(base_loss, sample_id)-> mask_coverage_summary

mask_coverage_summary %>%
  kableExtra::kable(format = "pipe")

mask_coverage_summary %>%
filter(base_loss > 0) %>%
ggplot(aes(x = base_loss)) + geom_histogram()

s <- mask_coverage_summary$sample_id[1]
g <- mask_coverage_summary$gene_id[1]
scores %>% 
  filter(sample_id == s) -> target
sq <- target$g_scores[[1]][g] %>% as.list %>% unlist %>% .[[1]]
gn <- gene_list[g]
# quick & dirty for testing only
x <- GPos(gn)
x$score <- 0
y <- AddScores(x, sq)
z <- y$score
plot(z)
library(changepoint.np)
a <- cpt.np(z)

library(HresSE)
qq <- readRDS("/Users/robertshear/Documents/n/groups/churchman/rds19/data/S004/hresse/HresSE_screen.rds")
b <- assay(qq, "cp")[g, 1][[1]]
cumsum(width(b))
width(gn)
# scores %>%
#   unnest(gmask) %>% 
#   mutate(jsim = map2_dbl(gsignal, gmask, jaccard_similarity),
#         npos_signal = map_int(gsignal, length), 
#         npos_mask =  map_int(gmask, length), 
#         mean_signal = map_dbl(gsignal, function(u) mean(u$score))) %>%
#   dplyr::select(!(bam_file:gsignal)) %>%
#   dplyr::select(!gmask) -> z
# 
# z %>% ggplot(aes(group = n_multi, x = n_multi, y = jsim)) +
#   geom_col() + facet_wrap(vars(sample_id))
# 
# scores %>%
#   unnest(gmask)%>% filter(n_multi == 4) %>% dplyr::select(sample_id, gmask) -> signals
# seq(nrow(signals)) %>% expand.grid(A=., B=.) %>% filter(A < B) %>%
#   mutate(jsim = map2_dbl(A, B, function(a, b)
#     jaccard_similarity(signals$gmask[[a]], signals$gmask[[b]]))) %>%
#   rbind(data.frame(A = 1:4, B = 1:4, jsim = 1.0)) -> zcor
# 
# pivot_wider(zcor, names_from = B, values_from = jsim)
# # TODO: Exclude wt-4?
# # TODO: look at mask jitter?
# n <- nrow(signals)
# z_mat <- matrix(ncol = n, nrow = n)
# z_mat[as.matrix(as.matrix(zcor[,1:2]))] = zcor$jsim
# rownames(z_mat) <- colnames(z_mat) <- signals$sample_id
# corrplot(z_mat, is.corr = FALSE, diag = FALSE, type = "upper", method = "pie")
# z_mat
