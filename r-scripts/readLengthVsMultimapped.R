# read length vs mmap
library(tidyverse)
library(GenomicAlignments)
library(GenomicRanges)

data_root <- "/n/groups/churchman/rds19/data/S006"
rds_filename <- file.path(data_root, "GPosExp-100Rand.rds")
e <- readRDS(rds_filename)
which <- rowRanges(e)
bam_fnames <- list.files(file.path(data_root),
                         pattern = "*[.]bam$", full.names = TRUE)

r <- tibble(fn = bam_fnames) %>%
  mutate(sampleId = map_chr(fn, function(u)
    str_remove(basename(u), ".bam"))) %>%
  mutate(profile = map(fn, function(u) {
    GRanges(readGAlignments(u,
         param = ScanBamParam(tag = c("NH", "HI"),
                  what = c("qname", "cigar", "qwidth"),
                  which = which))) %>%
      as_tibble %>%
      group_by(NH, width) %>%
      summarise(n = n())
  })) %>% 
  mutate(totmask = apply(GPosExperiment::mask(e), 2, 
                         function(u) sum(sum(width(GRangesList(u))))),
         cover_mask = map_dbl(totmask, 
                        function(u) round(u / sum(width(which)) * 100, 2)))

# TODO
# Cumm Sum of coverage by mask width (decreasing)

# Cover_mask is ratio of multimapped multimapped to unqiue mapped reads
r %>%
  unnest(profile) %>%
  ggplot(aes(x = width, y = n)) +
    geom_col(position = "dodge", aes(group = multimapped, fill = multimapped)) +
    facet_wrap(facets = vars(sampleId, cover_mask), scales = "free_y",
                           labeller = "label_both")

r %>%
  unnest(profile) %>%
  filter(multimapped == TRUE) %>%
  ggplot(aes(x = width, y = n)) +
  geom_col() +
  facet_wrap(facets = vars(sampleId, cover_mask), scales = "free_y",
             labeller = "label_both")
