# searchForTrisomy.R
library(GenomicRanges)
library(rtracklayer)
library(tidyverse)


total_counts <- function(u) {
  a <- import(u, genome = "sacCer3")
  sapply(split(a, seqnames(a)), function(u) sum(u$score * width(u)))[1:16]
}

bg_directory <- "/n/groups/churchman/rds19/data/S006/"
fn_extension <- "[.]pos[.]bedgraph[.]gz"
f <- tibble(bg_file = list.files(bg_directory, 
                                  str_glue("*", fn_extension, "$"), 
                                  full.names = TRUE),
            sample_id = map_chr(bg_file, function(u) str_replace(basename(u), 
                                    fn_extension, "")),
            strain = map_chr(sample_id, function(u) 
              str_split_fixed(u, "-", 2)[,1]),
            n_x_chr = map(bg_file, function(u) {
              result <- total_counts(u) + 
                total_counts(str_replace(u, fixed(".pos."), ".neg."))
              })
            # TODO group by strain
            )

m <- do.call(rbind, f$n_x_chr)
m_cor <- cor(t(m))
print(min(m_cor))
# Conclusion: no anuploidy
