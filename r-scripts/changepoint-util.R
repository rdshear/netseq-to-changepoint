# From Blackwell2017
find.intervals <- function(x, prob = 0.95) {
  if (sum(x) < prob) return(NULL)
  max.ints <- floor(sum(x)/prob)
  one.points <- which(x > prob)
  if (length(one.points) == max.ints) return(one.points)
  xhold <- x
  xhold[one.points] <- NA
  dv <- div.vector(x, one.points)
  max.w <- max(sapply(dv, length))
  if (max.w == length(x)) max.w <- max.w - 1
  count <- length(one.points) + 1
  out <- as.list(one.points)
  for (w in 1:max.w) {
    for (i in 1:(length(x) - w)) {
      this.sum <- sum(xhold[i:(i+w)])
      if (is.na(this.sum)) this.sum <- 0
      if (this.sum > prob) {
        out[[count]] <- i:(i+w)
        xhold[i:(i+w)] <- NA
        count <- count + 1
      }
    }
  }
  out
}

div.vector <- function(x, cuts) {
  if (length(cuts) == 0) return(list(x))
  out <- list()
  out[[1]] <- x[1:(cuts[1] - 1)]
  for(i in 1:(length(cuts))) {
    lo <- cuts[[i]] + 1
    hi <- ifelse(i == length(cuts), length(x), cuts[[i+1]]-1)
    out[[i+1]] <- x[lo:hi]
  }
  out
}
