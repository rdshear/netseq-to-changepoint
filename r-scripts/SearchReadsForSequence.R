library(Biostrings)
library(rtracklayer)

# chrIV:437860-437879, … TTATATTCCCAATATTTTT - 1 match
# chrIV:437831-437850 … AGGTTCGAGTCCTGCAGTT - 11 matches

targets = c("TTATATTCCCAATATTTTT", "AGGTTCGAGTCCTGCAGTT")
targets <- DNAStringSet(targets)
targets <- c(targets, reverseComplement(targets))
infile <- "/n/groups/churchman/GSE159603/wt-1.fastq.gz"

reads <- readLines(infile)
seqs <- reads[rep_len(c(FALSE, TRUE, FALSE, FALSE), length(reads))]

result <- lapply(targets, function(u) {
  result <- grep(u, seqs, fixed = TRUE)
  cat(sprintf("target=%s  count=%d\n", u, length(result)))
  result
})
