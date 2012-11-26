claiminfo <- function(...) {
    arglist <- list(...)

    if (length(arglist) == 1L && is.riskproc(arglist[[1L]])) {
        ## If we have a riskproc object, just return its claiminfo part
        return(arglist[[1L]][['claims']])
    } else {
        if ('hypoexp' %in% names(arglist) && is.numeric(arglist[['hypoexp']][['rates']])) {

            arglist <- within(arglist, {
                mu           <- sum(1.0 / arglist[['hypoexp']][['rates']])
                hypoexp$coef <- ratetoalpha(arglist[['hypoexp']][['rates']])

                mgf <- function(x) {
                    mgfhypoexp(x         = x,
                               rate      = arglist[['hypoexp']][['rates']],
                               difforder = 0L)
                }

                mgf.d1 <- function(x) {
                    mgfhypoexp(x         = x,
                               rate      = arglist[['hypoexp']][['rates']],
                               difforder = 1L)
                }

                mgf.d2 <- function(x) {
                    mgfhypoexp(x         = x,
                               rate      = arglist[['hypoexp']][['rates']],
                               difforder = 2L)
                }

                cdf <- function(x) {
                    phypoexp(q    = x,
                             rate = arglist[['hypoexp']][['rates']])
                }

                cdf.tailarea <- function(x) {
                    phypoexp(q        = x,
                             rate     = arglist[['hypoexp']][['rates']],
                             tailarea = TRUE)
                }

                pdf <- function(x) {
                    dhypoexp(x    = x,
                             rate = arglist[['hypoexp']][['rates']])
                }
            })
        }

        return(structure(.Data = arglist,
                         class = c('claiminfo', 'list')))
    }
}
