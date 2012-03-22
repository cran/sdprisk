lundbergRuinprob <- function(process, use.factor = FALSE) {
    stopifnot(is.logical(use.factor))

    .adjcoef <- adjcoef(process)
    if (is.null(.adjcoef)) {
        stop('Unable to compute the adjustment coefficient.')
    }

    .const <- 1
    if (use.factor) {
        try({
            .pdf <- process$claims$pdf            
        }, silent = TRUE)
        if (!is.null(.pdf)) {
            .const <- get('p', process) / (.adjcoef * (get('q', process) / mean(process$claim) *
                      stats::integrate(f = function(.x) {
                                            .res <- .pdf(.x) * .x * exp(.x * .adjcoef)
                                            .res[which(is.nan(.res))] <- 0
                                            return(.res)
                                       },
                                       lower = 0,
                                       upper = Inf)$value
                      + 1 / get('zeta', process)
                      ))
        } else {
            warning('Cannot compute the coefficient without the claim density.\n',
                    'Using 1 as a replacement as with use.factor = FALSE.')
        }
    }

    .rp <- function(x) {
        return(.const * exp(-x * .adjcoef))
    }

    return(structure(.rp, C = .rp, r = .adjcoef))
}
