library(GPosExperiment)
library(rtracklayer)
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
# count of censored positions
assay(e, "n_censored") <- apply(s, 1:2, function(u) sum(is.na(u[[1]]$score)))

assay(e, "count") <- n_s <- apply(s, 1:2, function(u) sum(u[[1]]$score, na.rm = TRUE))

assay(e, "n_censored")
assay(e, "count")

