claiminfo <- function(...) {
    arglist <- list(...)

    if (length(arglist) == 1 && is.riskproc(arglist[[1]])) {
        ## If we have a riskproc object, just return its claiminfo part
        return(get('claims', arglist[[1]]))
    } else {
        if ("hypoexp" %in% names(arglist) && is.numeric(arglist$hypoexp$rates)) {

            arglist <- within(arglist, {
                mu <- sum(1 / arglist$hypoexp$rates)
                hypoexp$coef <- ratetoalpha(arglist$hypoexp$rates)

                mgf <- function(x) {
                    mgfhypoexp(x, rate = arglist$hypoexp$rates, difforder = 0)
                }

                mgf.d1 <- function(x) {
                    mgfhypoexp(x, rate = arglist$hypoexp$rates, difforder = 1)
                }

                mgf.d2 <- function(x) {
                    mgfhypoexp(x, rate = arglist$hypoexp$rates, difforder = 2)
                }

                cdf <- function(x) {
                    phypoexp(x, rate = arglist$hypoexp$rates)
                }

                cdf.tailarea <- function(x) {
                    phypoexp(x, rate = arglist$hypoexp$rates, tailarea = TRUE)
                }

                pdf <- function(x) {
                    dhypoexp(x, rate = arglist$hypoexp$rates)
                }
            })
        }

        return(structure(arglist, class = c('claiminfo', 'list')))
    }
}
