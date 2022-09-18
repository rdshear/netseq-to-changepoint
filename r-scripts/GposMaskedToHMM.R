# GposMaskedToHMM
library(GenomicRanges)
library(rtracklayer)
library(GPosExperiment)
library(AnVIL)
library(fitdistrplus)
library(tidyverse)
library(kableExtra)
library(MCMCpack)
source("./r-scripts/changepoint-util.R")
library(rethinking) # from McElreath "Rethinking" book

cp2020 <- import.gff("./test/wt-1.cp.base.gff3")
names(cp2020) <- cp2020$tx_name
cpt <- tibble(id = cp2020$tx_name, index = cp2020$seq_index, length = width(cp2020)) %>%
  group_by(id) %>% nest()

# HACK One sample only!
samples = c("wt-1")

v <- readRDS("./test/GPosMasked.rds")
rownames(v) <- rowData(v)$ID
dim(v)
v <- v[cpt$id, samples]

x.bar <- assay(v, "x.bar")
pum <- assay(v, "prop.unmasked")
nomasks <- names(pum[rowMax(pum) == 0,])
x <- x.bar[nomasks,]
# HACK for one column only!
genes <- names(x[x > 0.5])
u <- v[genes, ]
cpt <- cpt %>% filter(id %in% genes)
# 
# target.genes <- names(head(sort(s.len[s.len > 500]), n = 7))

# The target.genes have no masked positions, average density > 0.5 and a lenght ~ 500nt
# df <- data.frame(score = s[target.genes[1], 1][[1]]$score)

target.genes <- "YCL059C"
cpt %>% filter(id == target.genes) %>% 
  unnest("data") %>% ungroup() %>%
  transmute(bp = cumsum(length)) %>% unlist() -> ref_bp

df <- data.frame(score = scores(v[target.genes[1], samples])[[1]]$score)
df$x <- seq(length(df$score))
df <- df[!is.na(df$score),]

timestamp()
posterior <- HDPHSMMnegbin(score ~ 1, data = df, K = 10, verbose = 500,
                           a.omega = 1, b.omega = 100, rho.step = 0.1)
timestamp()
plotHDPChangepoint(posterior, ylab="Density")

post.mean <- post.rho <- rep(NA, times = nrow(df))
post.rho <- matrix(NA, nrow = nrow(posterior), ncol = nrow(df))
post.beta <- matrix(NA, nrow = nrow(posterior), ncol = nrow(df))
for (j in 1:nrow(df)) {
  inds <- c(1:nrow(posterior)) + (nrow(posterior)) * (attr(posterior, "s.store")[,j]-1)
  post.beta[,j] <- posterior[inds]
  post.rho[,j] <- attr(posterior, "rho.store")[inds]
}
ch.pt <- rbind(0, diff(t(attr(posterior, "s.store"))) != 0)
pr.st <- rowMeans(ch.pt)
post.mean <- colMeans(post.beta)
nb <- names(which.max(table(attr(posterior, "num.regimes"))))
nb <- as.numeric(nb) - 1
ints <- find.intervals(pr.st, prob = 0.89)
par(mfrow=c(2,1))
plot(pr.st, type = "l", lwd = 3)
for (i in seq_along(ints)) {
print(c(range(ints[[i]]), diff(range(ints[[i]]))))
  abline(v = range(ints[[i]]), col = i + 1, lwd = 2)
}
plot(df$x,log(df$score), type = "h", col = "pink")
abline(v = ref_bp, col = "grey", lwd = 3)
ref_bp

# posterior mean beta
plot(post.mean, type = "l", col = "blue")
abline(v = ref_bp, col = "grey", lwd = 1)

