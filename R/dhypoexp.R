dhypoexp <- function(x, rate = 1, log = FALSE) {
    res <- drop(outer(x, rate, stats::dexp) %*% ratetoalpha(rate))
    if (log) {
        return(log(res))
    } else {
        return(res)
    }
}
