boundsRuinprob <- function(process, interval, maxreserve, richardson = TRUE, use.splines = FALSE) {
    stopifnot(is.riskproc(process), is.numeric(interval), is.numeric(maxreserve),
              is.logical(richardson), is.logical(use.splines))

    .p <- get('p', process)
    .q <- get('q', process)
    .zeta <- get('zeta', process)
    .claims <- get('claims', process)

    .n <- floor(maxreserve / interval) + 3

    try(.mytailarea <- get('cdf.tailarea', .claims), silent = TRUE)
    if (exists('.mytailarea')) {
        h2l <- diff(sapply(seq(0, by = interval, length.out = .n + 1), .mytailarea))
    } else {
        warning('Integrated tail area distribution function not supplied.\n',
                'Trying to use a rough approximation based on the CDF.\n',
                immediate. = TRUE, call. = FALSE)
        .mu <- mean(.claims)
        try(mycdf <- get('cdf', .claims), silent = TRUE)
        if (!is.na(.mu) & exists('mycdf')) {
            h2l <- sapply(seq(interval, by = interval, length.out = .n),
                          function(.y) { 1 - mycdf(.y) }) * interval / .mu
        } else {
            stop('Claim CDF or mean claim size not available. Please refer to the help.',
                 call. = FALSE)
        }
    }

    h1l <- diff(sapply(seq(0, by = interval, length.out = .n + 1), pexp, rate = .zeta))

    rp <- .C('rpbounds',
             h1l = as.double(h1l),
             h1u = as.double(c(0, h1l)),
             h2l = as.double(h2l),
             h2u = as.double(c(0, h2l)),
             q = as.double(.q),
             n = as.integer(.n),
             fl = double(.n),
             fu = double(.n),
             lowerbound = double(.n),
             upperbound = double(.n)
             )

    .x <- seq(from = 0, by = interval, along.with = rp$lowerbound)

    psi.lower <- stepfun(.x, c(rp$lowerbound, 0), right = TRUE)
    psi.upper <- stepfun(.x, c(1, rp$upperbound), right = FALSE)
    if (use.splines) {
        psi <- splinefun(x = .x, y = (rp$lowerbound + rp$upperbound) / 2)
    } else {
        psi <- approxfun(x = .x, y = (rp$lowerbound + rp$upperbound) / 2,
                         method = 'linear', rule = 2)
    }

    psi.x <- (rp$lowerbound + rp$upperbound) / 2
    diff1 <- -diff(psi.x) / (interval * .zeta * .p)

    if (richardson) {
        diff3 <- -diff(psi.x, lag = 3) / (3 * interval * .zeta * .p)
        diff1 <- rowSums(cbind(
            c(1, rep.int( 9 / 8, .n - 3)) * diff1[1 - .n],
            c(0, rep.int(-1 / 8, .n - 3)) * c(0, diff3)
        ))
    }

    .x.1 <- c(0, .x[1:(.n - 3)] + interval / 2)
    if (use.splines) {
        psi.1 <- splinefun(x = .x.1, y = c(1, diff1[1:(.n - 3)]))
    } else {
        psi.1 <- approxfun(x = .x.1, y = c(1, diff1[1:(.n - 3)]), method = 'linear', rule = 2)
    }

    return(
        structure(
            list(psi = psi, psi.1 = psi.1, psi.2 = function(x) { psi(x) - psi.1(x) },
                 psi.lower = psi.lower, psi.upper = psi.upper),
            compmethod  = 'bounds',
            riskproc    = process,
            parameters  = list(interval = interval, maxreserve = maxreserve),
            diagnostics = list(rp = rp)
        )
    )
}
