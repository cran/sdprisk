fftRuinprob <- function(process, interval, maxreserve, n, use.splines = FALSE) {
    stopifnot(is.logical(use.splines))

    .misspar <- c(missing(interval), missing(maxreserve), missing(n))

    if (sum(.misspar) > 1) {
        stop('At least two parameters of ', sQuote('interval'), ', ',
             sQuote('maxreserve'), ' or ', sQuote('n'), ' are required.')
    } else {
        if (missing(n)) {
            n <- nextn(maxreserve / interval)
        }
        if (missing(maxreserve)) {
            maxreserve <- n * interval
        }
        if (missing(interval)) {
            interval <- maxreserve / n
        }
        stopifnot(is.numeric(interval), is.numeric(maxreserve), is.numeric(n))
    }

    .p <- get('p', process)
    .q <- get('q', process)
    .zeta <- get('zeta', process)
    .claims <- get('claims', process)

    .x <- seq(0, by = interval, length.out = n + 1)

    a <- diff(pexp(q = .x, rate = .zeta))
    b <- diff(sapply(.x, .claims$cdf.tailarea))

    a.fft <- fft(a)
    w <- Re(fft(.p * a.fft / (1 - .q * a.fft * fft(b)), inverse = TRUE) / n)

    .x <- .x[-(n + 1)]

    if (use.splines) {
        psi   <- splinefun(.x, rev(cumsum(rev(w))))
        psi.1 <- splinefun(.x, w / (.zeta * .p * interval))
    } else {
        psi   <- approxfun(.x, rev(cumsum(rev(w))), method = 'linear', rule = 2)
        psi.1 <- approxfun(.x, w / (.zeta * .p * interval), method = 'linear', rule = 2)
    }

    return(
        structure(
            list(psi = psi, psi.1 = psi.1, psi.2 = function(x) { psi(x) - psi.1(x) }),
            compmethod  = 'fft',
            riskproc    = process,
            parameters  = list(interval    = interval,
                               maxreserve  = maxreserve,
                               n           = n,
                               use.splines = use.splines),
            diagnostics = list(w = w)
        )
    )
}
