rhypoexp <- function(n = 1L, rate = 1.0) {
    if ((nlen <- length(n)) > 1L) {
        n <- nlen
    }
    rowSums(vapply(X         = rate,
                   FUN       = stats::rexp,
                   FUN.VALUE = numeric(n),
                   n         = n))
}
