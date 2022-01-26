# DiscoverBreakpoints.R
# bedgraph files --> gff3 files
# disregarding multi-mapped regions
# 
# Usage:
# Rscript --vanilla scripts/DiscoverBreakpoints.R \
#   {reference.gene.list.filename} \
#   {max.gene.length} \
#   {K.max} \
#   {max.genes} \
#   {algorithm} \
#   {pos.begraph.filename} \
#   {neg.begraph.filename} \
#   {output.file name}

# To run with embedded parameters, set DEBUG.TEST <- TRUE

suppressPackageStartupMessages({
  library(parallel)
  library(GenomicRanges)
  library(rtracklayer)
  library(breakpoint)
})
set.seed(20190416)

# DEBUG ONLY FROM HERE.....
DEBUG.TEST <- TRUE
if (interactive() && exists("DEBUG.TEST")) {
  print("DEBUG IS ON -- COMMAND LINE PARAMETERS IGNORED")
  commandArgs <- function(trailingOnly) {
    c("/n/groups/churchman/rds19/data/S005/genelist.gff",
      "800", # Maximum gene body length
      "8", # Kmax (maximum number of segments)
      "256", # Sample size
      "CEZINB",
      "/n/groups/churchman/rds19/data/S005/wt-2.pos.bedgraph.gz",
      "/n/groups/churchman/rds19/data/S005/wt-2.neg.bedgraph.gz",
      "/n/groups/churchman/rds19/data/S005/wt-2.CEZINB.gff3")

  }
}
#  .............TO HERE

args <- commandArgs(trailingOnly = TRUE)

subject_genes.filename <- args[1]
GeneMaxLength <- as.numeric(args[2]) # Truncate gene to this length (or inf if 0)
Kmax <- as.numeric(args[3])
maxGenes <- as.numeric(args[4]) # if > 0 sample this number of genes
sample.name <- args[5]
bedgraph.pos.filename <- args[6]
bedgraph.neg.filename <- args[7]
out.filename <- args[8]

start.time <- Sys.time()
cat(sprintf("Starting at %s. Sample name = %s. Max genes = %d", 
        start.time, sample.name, maxGenes))

algorithm <- "CEZINB"

cores <- detectCores()
options(mc.cores = max(cores - 1, 1))
cat(sprintf("Number of cores detected = %d. Cores to use = %d.", cores, getOption("mc.cores")))


# read the two bedgraph files into a single GRanges object
scores <- mapply(function(strand_sym, infilename) {
  x <- import(infilename, genome = "sacCer3")
  strand(x) <- strand_sym
  x
}, list('+', '-'), list(bedgraph.pos.filename, bedgraph.neg.filename), SIMPLIFY = FALSE)
scores <- unlist(GRangesList(scores))

# now make the GRanges object dense (zeros instead of gaps)
x <- gaps(scores)
x$score <- 0
scores <- sort(c(scores, x))
scores <- scores[strand(scores) != "*"]

# the genes to query
g <- import(subject_genes.filename, genome = 'sacCer3')

# truncate gene lengths if so desired
if (GeneMaxLength > 0) {
  g <- resize(g, fix = "start", ifelse(width(g) > GeneMaxLength, GeneMaxLength, width(g)))
}

# subset the genes if so desired
if (maxGenes > 0 & maxGenes < length(g)) {
  g <- g[sort(sample(length(g), maxGenes))]
}


h <- findOverlaps(g, scores)
h <- split(subjectHits(h), queryHits(h))

result <- mclapply(seq_along(g), function(i) {
    u <- g[i]
    x <- scores[h[[i]]]
    start(x[1]) <- start(u)
    end(x[length(x)]) <- end(u)
    s <- rep(x$score, width(x))
    
    mu <- mean(s)
    v <- var(s)
    # heuristic. For the CE.ZINB for mean occupancy < 0.05 reads / NT
    if (mu < 0.05) {
      tau <- numeric(0)
    } else
    {
      # NOTE: This is *always* + direction. We assume no bias in this algorithm by strand
      seg <- try(CE.ZINB(data = data.frame(s), Nmax = Kmax, parallel = FALSE), silent = TRUE)
      if (inherits(seg, c("character", "try-error"))) {
        tau <- numeric(0)
      } else {
        tau <- seg$BP.Loc
      }
    }
    n <- length(tau) + 1
    
    tr <- cbind(c(1, tau), c(tau - 1, length(s)))
    
    stats <- apply(tr, 1, function(w) {
      v <- s[w[1]:w[2]]
      c(m = mean(v), v = var(v))
    })
    if (as.vector(strand(u)) == "+") {
      sqi <- 1:n
    } else
    {
      sqi <- n:1
    }
    result <- GRanges(seqnames = seqnames(u),
                strand = strand(u),
                ranges = IRanges(start = start(u) - 1 + tr[, 1], 
                  end = start(u) - 1 + tr[, 2]),
                seq_index = sqi,
                type = "seq_index",
                source = "DiscoverBreakpoints2",
               algorithm = algorithm,
               tx_name = as.character(u$Name),
               m = stats["m", ],
               v = stats["v", ])

    # TODO LOG
    # writeLines(kableExtra::kable(result, format = "pipe"), stderr())
    # writeLines(as.character(as.numeric(Sys.time() - time.in, unit = "secs")), stderr())
    # writeLines("---", stderr())
    # 
    result
  })

result <- unlist(GRangesList(result))
result$m <- round(result$m, 3)
result$v <- round(result$v, 3)
result <- sort(result)

  # # TODO: Carry BIC and logLikelihood
  # # # TODO: report multi-map removal areas

export(result, con = out.filename)

end.time <- Sys.time()
run.time <-  as.numeric(end.time - start.time, units = "secs")
n.genes <- length(g)
n.segments <- length(result)
n.bases <- sum(width(g))

cat(sprintf("Completed at %s\nRun time: %.3f sec\ngenes: %.0f \nbases: %.0f\nsec / gene: %.3f\nmsec/base: %.3f", 
            as.character(end.time), run.time, n.genes, n.bases, run.time / n.genes, run.time * 1000 / n.bases))


if (!interactive()) {
  print(sessionInfo())
}
