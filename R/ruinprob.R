ruinprob <- function(process, method = c('saddlepoint', 'fft', 'bounds', 'hypoexp', 'lundberg'), ...) {
    stopifnot(is.riskproc(process))
    method <- match.arg(method)

    switch(method,
        `saddlepoint` = saddlepointRuinprob(process, ...),
        `fft`         = fftRuinprob(process, ...),
        `bounds`      = boundsRuinprob(process, ...),
        `hypoexp`     = hypoexpRuinprob(process),
        `lundberg`    = lundbergRuinprob(process, ...)
    )
}
