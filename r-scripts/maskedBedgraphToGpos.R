library(GenomicRanges)
library(rtracklayer)
library(GPosExperiment)
library(tidyverse)
library(breakpoint)

set.seed(20220311)

matrix_apply <- function(x, f, ...) {
  result <- sapply(x, f, ...)
  dim(result) <- dim(x)
  dimnames(result) <- dimnames(x)
  result
}

data_root <- "/n/groups/churchman/rds19/data/S006"
# rds_filename <- file.path(data_root, "GPosExp-Uzun.rds")
rds_filename <- file.path(data_root, "GPosExp-100Rand.rds")
n_genes <- 100

which <- import(file.path(data_root, "genelist.gff"), genome = "sacCer3")
names(which) <- which$ID

# gx <- c("YDR152W", "YDR311W", "YDR381W") # examples  (Uzun et al. 2021)
# which <- which[gx]
n_genes <- 100
which <- sort(sample(which, n_genes))

begraphpos_fnames <- list.files(file.path(data_root),
                         pattern = "*[.]pos[.]bedgraph[.]gz$", full.names = TRUE)
begraphneg_fnames <- list.files(file.path(data_root),
                                pattern = "*[.]neg[.]bedgraph[.]gz$", full.names = TRUE)

mask_pos <- list.files(file.path(data_root),
                       pattern = "*[.]mask_pos[.]bedgraph[.]gz$", full.names = TRUE)

mask_neg <- list.files(file.path(data_root),
                       pattern = "*[.]mask_neg[.]bedgraph[.]gz$", full.names = TRUE)

sample_ids <- basename(sapply(strsplit(begraphpos_fnames, ".", fixed = TRUE), function(u) u[1]))

nsd <- NETseqDataFromBedgraph(sample_ids,
                            begraphpos_fnames, begraphneg_fnames,
                            seqinfo = seqinfo(which))

# TODO need "which" parameter for GPosExperiment
# TODO add masks to the constructor
for (i in seq_along(nsd)) {
  m1 <- import(mask_pos[i])
  m2 <- import(mask_neg[i])
  strand(m1) <- "+"
  strand(m2) <- "-"
  mask(nsd[[i]]) <- sort(c(m1, m2))
}
e <- GPosExperiment(nsd, rowRanges = which)
s <- scores(e, apply_mask = TRUE, zero_fill = TRUE)
n_mask <- matrix_apply(s, function(u) sum(is.na(u$score)))

n_s <- matrix_apply(s, function(u) sum(u$score, na.rm = TRUE))
n_mask
n_s

# consider s[2,1], which has a single lost region
# subject <- s[2,1][[1]]
# s_na_idx <- which(is.na(subject$score))
# subject_censored <- subject[-s_na_idx]
# sum(subject_censored$score)
# plot(pos(subject_censored), subject_censored$score)
# abline(v = pos(subject[s_na_idx]), col = "red")

calc_cp <- function(subject) {
  trace <<- trace + 1
  print(trace)
  s_na_idx <- which(is.na(subject$score))
  if (length(s_na_idx) == 0) {
    subject_censored <- subject
  } else {
    subject_censored <- subject[-s_na_idx]
  }
  
  seg <- try(CE.NB(data = data.frame(subject_censored$score),
                   Nmax = Kmax, parallel = FALSE), silent = FALSE)
  if (inherits(seg, c("character", "try-error"))) {
    seg <- list(bpts = GRanges())
  } else {
    list(bpts = subject_censored[seg$BP.Loc], BIC = seg$BIC, ll = seg$ll)
  }
}

x.bar <- matrix_apply(s, function(u) mean(u$score, na.rm = TRUE))
assay(e, "x.bar") <- x.bar
prop.unmasked <- apply(GPosExperiment::mask(e), 2, function(u) sum(width(GRangesList(u))) / width(rowRanges(e)))
assay(e, "prop.unmasked") <- prop.unmasked

saveRDS(e, rds_filename)

library(breakpoint)
Kmax <- 10

# TODO breakout to sepearte source file
trace <- 0
result <- lapply(s, calc_cp)
bpts <- matrix_apply(result, function(u) start(u$bpts))
dim(bpts) <- dim(s)
dimnames(bpts) <- dimnames(s)
assay(e, "bpts") <- bpts
bic <-  matrix_apply(result, function(u) u$BIC)
dim(bic) <- dim(s)
dimnames(bic) <- dimnames(s)
assay(e, "bic") <- bic
assay(e, "k") <- matrix_apply(bpts, length)

saveRDS(e, rds_filename)

bp <- assay(e, "bpts")

par(mfcol = c(ncol(e),nrow(e)), mar = c(2, 2, 1, 1) + 0.1)
for (i in 1:nrow(e))
  for (j in 1:ncol(e)) {
    u <- s[i,j][[1]]
    plot(pos(u), log2(u$score), type = "h", col = "pink")
    censored <- pos(u[is.na(u$score)])
    abline(v = censored, col = "green")
    v <- bp[i,j][[1]]
    abline(v = v, col = "blue")
  }

#e <- readRDS(rds_filename)

# result <- matrix_apply(s, calc_cp)


# abline(v = pos(result$bpts), col = "purple")
# bp.zi <- CE.ZINB(data = data.frame(subject_censored$score), Nmax = Kmax, parallel = FALSE)
# s_bp.zi <- subject_censored[bp$BP.Loc]
# abline(v = pos(s_bp), col = "blue")
# 
# bp <- CE.NB(data = data.frame(subject_censored$score), Nmax = Kmax, parallel = FALSE)
# s_bp <- subject_censored[bp$BP.Loc]
# abline(v = pos(s_bp), col = "blue")
# 
# library(fitdistrplus)
# fd <- fitdist(subject_censored$score, distr = "nbinom")
# 
# cd <- subject$score
# cd[is.na(cd)] <- 0
# fdc <- fitdistcens(cd, distr = "nbinom")
