rhypoexp <- function(n = 1, rate = 1) {
    nlen <- length(n)
    if (nlen > 1) {
        n <- nlen
    }
    rowSums(cbind(sapply(rate, stats::rexp, n = n)))
}
