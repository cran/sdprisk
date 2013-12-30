rhypoexp <- function(n = 1L, rate = 1.0) {
    if ((nlen <- length(n)) > 1L) {
        n <- nlen
    }
    rowSums(matrix(vapply(X         = rate,
                          FUN       = rexp,
                          FUN.VALUE = numeric(n),
                          n         = n),
                   nrow = n,
                   ncol = length(rate)))
}
