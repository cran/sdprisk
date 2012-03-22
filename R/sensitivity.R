sensitivity <- function(process, method = c('saddlepoint', 'hypoexp'), ...) {
    stopifnot(is.riskproc(process))
    method <- match.arg(method)

    switch(method,
        `saddlepoint` = saddlepointSensitivity(process, ...),
        `hypoexp`     = hypoexpSensitivity(process, ...)
    )
}
