# DiscoverBreakpointsWorker.R
# 
# Usage:
# Rscript --vanilla scripts/DiscoverBreakpointsGather.R \
#   {changepoints.gff}    # combined output
#   {shard1.gff} {shard2.gff} ... {shardn.gff}
#

# To run with embedded parameters, set DEBUG.TEST <- TRUE

suppressPackageStartupMessages({
  library(rtracklayer)
})

# DEBUG ONLY FROM HERE.....
DEBUG.TEST <- TRUE
if (interactive() && exists("DEBUG.TEST")) {
  print("DEBUG IS ON -- COMMAND LINE PARAMETERS IGNORED")
  commandArgs <- function(trailingOnly) {
    c("~/temp/small_wt-1.cp.gff",
      "~/temp/cp_shard_0.gff"
      )
  }
}
#  .............TO HERE

args <- commandArgs(trailingOnly = TRUE)

output.filename <- args[1]
shards <- args[-1]

result <- unlist(GRangesList(lapply(shards, function(u) {
    v <- import.gff3(u, genome = "sacCer3")
    # HACK Because no good way to close the connection
    # See https://github.com/lawremi/rtracklayer/issues/14
    gc()
    v
  })))
result <- sort(result)
export.gff3(result, con = output.filename)

sprintf("Completed at %s\n", Sys.time())
print(sessionInfo())

