mean.claiminfo <- function(x, ...) {
    mymu <- NA_real_
    try(mymu <- get('mu', x), silent = TRUE)
    return(mymu)
}
