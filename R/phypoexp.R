phypoexp <- function(q, rate = 1, lower.tail = TRUE, log.p = FALSE, tailarea = FALSE) {
    mycoef <- ratetoalpha(rate)
    if (tailarea) {
        mycoef <- mycoef / (rate * sum(1 / rate))
    }
    res <- drop(
        outer(q, rate, stats::pexp, lower.tail = lower.tail, log.p = FALSE) %*% mycoef
    )
    if (log.p) {
        return(log(res))
    } else {
        return(res)
    }
}
