# DiscoverBreakpointsWorker.R
# bedgraph files --> gff3 files
# appropriate for scatter/gather model
# 
# Usage:
# Rscript --vanilla scripts/DiscoverBreakpointsWorker.R \
#   {shard_n.rds}    # changepoint input
#   {changepoints.gff}    # changepoint output
#   {K.max}   # maximum number of changepoints
#   {algorithm} # name of algorithm to run
#   {threads} # if 1, then no multiprocessing, otheerwise # of threads, 0 for all - 1
#

# To run with embedded parameters, set DEBUG.TEST <- TRUE

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(tidyverse)
  library(plyranges)
  library(parallel)
  library(foreach)
  library(doMC)
})
set.seed(20210915)

# DEBUG ONLY FROM HERE.....
DEBUG.TEST <- TRUE
if (interactive() && exists("DEBUG.TEST")) {
  setwd("~/temp")
  print("DEBUG IS ON -- COMMAND LINE PARAMETERS IGNORED")
  commandArgs <- function(trailingOnly) {
    c("~/temp/shard_0.rds",
      "~/temp/cp_shard_0.gff",
      "12", # Kmax (maximum number of segments)
      "CEZINB",
      "1"
      )
  }
}
#  .............TO HERE

args <- commandArgs(trailingOnly = TRUE)

input.filename <- args[1]
output.filename <- args[2]
Kmax <- as.numeric(args[3])
algorithm <- args[4]
threads <- as.numeric(args[5])

# Select algorithm
# Each algorithm has a function at the getbp entry:
# f(s, kMax), where s is an integer array of scores and kMax is an integer, maximum number of breakpoints to consider
# the result is a named list with the vanue named tau an array of integers representing the breakpoints 
selected_algorithm <- switch(algorithm,
           CEZINB = list(lib_name = "breakpoint", 
                   getbp = function(s, Kmax) {
                     u <- data.frame(score = s)
                     # NOTE: This is *always* + direction. We assume no bias in this algorithm by strand
                     seg <- try(CE.ZINB(data = u, Nmax = Kmax, parallel = FALSE), silent = FALSE)
                     if (inherits(seg, c("character", "try-error"))) {
                       seg <- list(tau = numeric(0))
                     }
                     list(tau = seg$BP.Loc, BIC = seg$BIC, ll = seg$ll)
                   }
                 )
          )

if (is.null(selected_algorithm)) {
  stop(sprintf("algorithm '%s' not available", algorithm))
}

if (!require(selected_algorithm$lib_name, character.only = TRUE))
{
  stop(sprintf("Could not load package '%s' for algorithm '%s'", selected_algorithm$lib_name, algorithm))
}

get_change_points <- function(s) {
  start.time <- Sys.time()
  props <- selected_algorithm$getbp(s, Kmax)
  elapsed.time <- as.numeric(difftime(Sys.time(), start.time, units = "secs"))
  tau <- props$tau
  props$elapsed_time <- elapsed.time
  result <- tibble(seq_index = seq(length(tau) + 1),
         s.start = c(1, tau), s.end = c(tau - 1, length(s))) %>% 
    mutate(m = map2_dbl(s.start, s.end,
                        function(a, b, c) mean(c[a:b]), s),
                      var = map2_dbl(s.start, s.end,
                          function(a, b, c) var(c[a:b]), s))
  list(result = result, props = props)
}

# Select correct number of threads and register local multicore backend
if (threads == 0) {
  threads <- max(detectCores(), 2)
}
registerDoMC(cores = threads - 1)

start.time <- Sys.time()


print(sprintf("Starting at %s.\n  Shard name=%s.  Kmax=%d.   \nAlgorithm=%s.  Threads=%d", 
              Sys.time(), input.filename, Kmax, algorithm, threads))


result <- readRDS(input.filename) %>% head(3) %>%
  mutate(scores = pmap(list(start, end, data), function(s, e, d) {
      locs <- map2(d$start - s + 1, d$end - d$start + 1, function(u, v) u + seq(0, v-1))
      x <- rep(0, e - s + 1)
      # convert GRanges scores to an integer array. A loop sometimes beats obscurity
      for (i in seq_along(locs)) {
        x <- modify_at(x, locs[[i]], function(a,b) a + b, d$score[i])
      }
      x
      })) %>%
  mutate(mu = map_dbl(scores, mean), v = map_dbl(scores, var)) %>% 
  mutate(tmp = foreach(s = .$scores) %dopar% get_change_points(s)) %>% 
  mutate(segments = lapply(tmp, function(u) u[[1]]),
                  props = lapply(tmp, function(u) u[[2]]), tmp = NULL)

r <- result %>% unnest(segments)
gr.out <- GRanges(seqnames = r$seqnames,
                  strand = r$strand,
                  ranges = IRanges(start = r$s.start + r$start - 1, 
                                   end = r$s.end + r$start - 1),
                  tx_name = r$ID,
                  seq_index = r$seq_index,
                  type = "seq_index",
                  algorithm = algorithm,
                  m = round(r$m,4),
                  v = round(r$var,4))

end.time <- Sys.time()
elapsed.time <- difftime(end.time, start.time, units = "secs")
gene.count <- length(unique(result$ID))
time.per.gene <- elapsed.time / gene.count
total_width <- sum(result$end - result$start)
time.per.nt <- elapsed.time / total_width

  
  # # # TODO: reverse negative strand 
  # # IRanges(start = c(1, seg + 1), end = c(seg, length(s)))
  # # # TODO: report multi-map removal areas
  # c(tx_name = u$tx_name, mu = mu, v = v, 
  #   mzl = head(sort(runLength(s), decreasing = TRUE)))

export.gff3(gr.out, con = output.filename)

cat(sprintf("Elapsed time: %.0f sec,  genes: %.0f,   bases: %.0f \n  %0.2f sec/gene,   %0.1f msec/base \n Completed at %s\n",
        elapsed.time, gene.count, total_width, time.per.gene, time.per.nt * 1000, Sys.time()))

print(sessionInfo())
