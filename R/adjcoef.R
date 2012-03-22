adjcoef <- function(process) {
    stopifnot(is.riskproc(process))
    if ('adjcoef' %in% names(process)) {
        .adjcoef <- get('adjcoef', process)
    } else {
        .adjcoef <- NULL
    }

    if (!is.null(.adjcoef)) {
        return(.adjcoef)
    } else {
        if (is.hypoexp(process)) {
            .r <- attr(ruinprob(process, method = 'hypoexp'), 'diagnostics')$r
            .adjcoef <- min(Re(.r[which(Im(.r) == 0 & Re(.r) > 0)]))
        } else {
            .freq <- get('freq', process)
            .premium <- get('premium', process)
            .variance <- get('variance', process)
            .claims <- get('claims', process)
            .mgfMaxAggregateLoss <-  function(.x) {
                 .freq * .claims$mgf(.x) + 0.5 * .variance * .x^2 - .premium * .x - .freq
            }
            try({
                # Find the largest root to get an upper bound of the interval
                .adjcoef <- max(rootSolve::uniroot.all(
                    f     = .mgfMaxAggregateLoss,
                    lower = .Machine$double.eps^(1/4),
                    upper = 1e+8,
                ))
                # Find the smallest non-zero root with high accuracy
                .adjcoef <- min(rootSolve::uniroot.all(
                    f     = .mgfMaxAggregateLoss,
                    lower = .Machine$double.eps^(1/4),
                    upper = .adjcoef + 1,
                    n = 1e+5
                ))
            }, silent = TRUE)

            if (is.null(.adjcoef)) {
                stop('The adjustment coefficient could not be computed.')
            }
        }

        ## Now that we computed the thing, let us add it to the riskproc object supplied by the user.
        ## Attention: This modifies an object in the global environment!
        .process <- process
        .process$adjcoef <- .adjcoef
        assign(deparse(substitute(process)), .process, pos = globalenv())

        return(.adjcoef)
    }
}
