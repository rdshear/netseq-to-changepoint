library(rtracklayer)
library(GPosExperiment)
library(tidyverse)
library(breakpoint)

data_root <- "~/n/groups/churchman/rds19/data/S006"

which <- import(file.path(data_root, "genelist.gff"), genome = "sacCer3")
names(which) <- which$ID
gx <- c("YDR152W", "YDR311W", "YDR381W") # examples  (Uzun et al. 2021)
which <- which[gx]

bam_fnames <- list.files(file.path(data_root),
                         pattern = "^.*wt.*bam$", full.names = TRUE)
bam_fnames

nsd <- mapply(function(id, fn) {
  print(id)
  print(fn)
  NETseqDataFromBAM(id, 
                         fn,
                         which, seqinfo(which))
  },
             id = basename(bam_fnames),
             fn = bam_fnames)
str(nsd, max.level = 3)
e <- GPosExperiment(nsd, rowRanges = which)
s <- scores(e, apply_mask = TRUE, zero_fill = TRUE)
n_mask <- apply(s, 1:2, function(u) sum(is.na(u[[1]]$score)))

n_s <- apply(s, 1:2, function(u) sum(u[[1]]$score, na.rm = TRUE))
n_mask
n_s
assay(e, "count") <- n_s

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
# saveRDS(e, rds_filename)
# e <- readRDS(rds_filename)
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
