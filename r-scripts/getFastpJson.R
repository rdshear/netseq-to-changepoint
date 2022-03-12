# Get fastp run data
library(tidyverse)
library(jsonlite)
library(GPosExperiment)

set.seed(20220311)

matrix_apply <- function(x, f, ...) {
  result <- sapply(x, f, ...)
  dim(result) <- dim(x)
  dimnames(result) <- dimnames(x)
  result
}

data_root <- "/n/groups/churchman/rds19/data/S006"
rds_filename <- file.path(data_root, "GPosExp-100Rand.rds")
e <- readRDS(rds_filename)


fn <- list.files(path = data_root, pattern = "*[.]json$", full.names = TRUE)
x <- read_json(fn[1])

r <- tibble(fn) %>%
  mutate(sampleIds = map_chr(fn, function(u) 
      str_remove(basename(u), ".fastp.json")),
    j = map(fn, function(u) read_json(u)),
    mean_read_length = map_int(j, 
     function(u) u$summary$after_filtering$read1_mean_length),
    reads_mb = round(map_int(j, 
     function(u) u$read1_after_filtering$total_reads) / 1E6, 2),
    prop.unmasked =  colMeans(assay(e, "prop.unmasked")),
    totmask = apply(mask(e), 2, function(u) sum(sum(width(GRangesList(u))))))

print(sort(r$mean_read_length))
print(sort(r$prop.unmasked))
plot(prop.unmasked ~ mean_read_length, data = r)
lr <- glm(prop.unmasked ~ mean_read_length, data = r)
summary(lr)

plot(totmask ~ mean_read_length, data = r)
lr1 <- glm(totmask ~reads_mb + mean_read_length, data = r)
summary(lr1)

kableExtra::kable(select(r, c("sampleIds", "reads_mb", "mean_read_length", "totmask")), format = "pipe")
