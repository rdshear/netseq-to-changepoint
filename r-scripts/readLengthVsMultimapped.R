# read length vs mmap
library(tidyverse)
library(GenomicAlignments)

data_root <- "/n/groups/churchman/rds19/data/S006"

bam_fnames <- list.files(file.path(data_root),
                         pattern = "*[.]bam$", full.names = TRUE)

r <- tibble(fn = bam_fnames) %>%
  mutate(sampleId = map_chr(fn, function(u) 
    str_remove(basename(u), ".fastp.json"))) %>%
  mutate(profile = map(fn, function(u){
    GRanges(readGAlignments(u, 
         param = ScanBamParam(tag = c("NH", "HI"),
                  what = c("qname", "cigar", "qwidth"),
                  which = GRanges("chrI:1-230218")))) %>%
      as.tibble %>% 
      filter(HI == 1) %>% 
      mutate(multimapped = NH > 1) %>% 
      group_by(multimapped, width) %>% 
      summarise(n = n())
  }))
    


