## ---- include=FALSE-----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----libraries----------------------------------------------------------------
suppressWarnings({
  library(HRseqExperiment, quietly = TRUE)
  library(rtracklayer, quietly = TRUE)
})

## ----row-definitions----------------------------------------------------------
data_root <- file.path(system.file(package = "HRseqExperiment"), "extdata")
genelist = rtracklayer::import(file.path(data_root, "genelist.gff3.bgz"), genome = "sacCer3")
names(genelist) <- genelist$ID
# Keep only the gene name, no other columns are needed
mcols(genelist) <- list(gene = genelist$gene)

# The regions of interest. For demonstration purposes, take 3 genes
genelist <- genelist[c("YAL062W", "YAL060W", "YAL059W")]
genelist

## ----BAM-to-samples-----------------------------------------------------------
bam_files = c("wt-1" = "SRR12840066.bam", 
              "wt-2" = "SRR12840067.bam")


# Create a list of HRseqData objests
nsd <- mapply(function(id, fn) {
  HRseqDataFromBAM(id, 
                     fn,
                     genelist, 
                    seqinfo(genelist))
  },
             id = names(bam_files),
             fn = file.path(data_root, bam_files))

# and now create

## ----create-Experiment--------------------------------------------------------
e <- HRseqExperiment(nsd, rowRanges = genelist)
e

## ----investigate-samples------------------------------------------------------
s <- scores(e, apply_mask = TRUE, zero_fill = TRUE)
n_mask <- apply(s, 1:2, function(u) sum(is.na(u[[1]]$score)))

n_s <- apply(s, 1:2, function(u) sum(u[[1]]$score, na.rm = TRUE))
n_mask
n_s

## ----plots--------------------------------------------------------------------
plot(e)

## ----SessionInfo--------------------------------------------------------------
sessionInfo()

