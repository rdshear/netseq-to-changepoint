Sys.setenv(TAR = "/bin/tar")
options(repos = "https://cloud.r-project.org/")
# install.packages(c("breakpoint", "doMC", "foreach", "foreach", "kableExtra"),
#     quiet = TRUE)
BiocManager::install(c("plyranges", "Biostrings",
    "GenomicRanges", "rtracklayer",
    "breakpoint", "doMC", "foreach", "foreach", "kableExtra"),
     quietly = TRUE, update = TRUE, ask = FALSE)



#     BiocManager::install(c(
#     "AnVIL", "brio", "cli", "cpp11", "credentials", "devtools", "digest",
#     "fansi", "fs", "gert", "glue", "jsonlite", "knitr", "littler", "memoise",
#     "openssl", "pkgbuild", "pkgload", "remotes", "sessioninfo", "stringi",
#     "testthat", "usethis", "withr", "xml2"
#   ), update = TRUE, ask = FALSE)