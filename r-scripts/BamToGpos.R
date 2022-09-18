library(rtracklayer)
library(GPosExperiment)
library(tidyverse)
library(breakpoint)

matrix_apply <- function(x, f, ...) {
  result <- sapply(x, f, ...)
  dim(result) <- dim(x)
  dimnames(result) <- dimnames(x)
  result
}

data_root <- "/n/groups/churchman/rds19/data/S006"

which <- import(file.path(data_root, "genelist.gff"), genome = "sacCer3")
names(which) <- which$ID
gx <- c("YDR152W", "YDR311W", "YDR381W") # examples  (Uzun et al. 2021)
which <- which[gx]

bam_fnames <- list.files(file.path(data_root),
                         pattern = "*[.]bam$", full.names = TRUE)

nsd <- mapply(function(id, fn) {
  print(id)
  print(fn)
  NETseqDataFromBAM(id, 
                         fn,
                         which, seqinfo(which))
  },
             id = basename(bam_fnames),
             fn = bam_fnames)

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

library(breakpoint)
Kmax <- 7

calc_cp <- function(subject) {
  s_na_idx <- which(is.na(subject$score))
  if (length(s_na_idx) == 0) {
    subject_censored <- subject
  } else {
    subject_censored <- subject[-s_na_idx]
  }
  bp <- CE.NB(data = data.frame(subject_censored$score), Nmax = Kmax, parallel = FALSE)
  c(bpts = subject_censored[bp$BP.Loc], calc = bp)
}

result <- lapply(s, calc_cp)
bpts <- matrix_apply(result, function(u) pos(u["bpts"][[1]]))
dim(bpts) <- dim(s)
dimnames(bpts) <- dimnames(s)
assay(e, "bpts") <- bpts
assay(e, "k") <- matrix_apply(bpts, length)

rds_filename <- file.path(data_root, "GPosExp-Uzun.rds")
#saveRDS(e, rds_filename)
e <- readRDS(rds_filename)
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
