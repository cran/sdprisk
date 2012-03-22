hypoexpVaru <- function(process) {
    .ruinprob <- hypoexpRuinprob(process)
    .psi.r <- attr(.ruinprob, 'diagnostics')$r
    .psi.C <- attr(.ruinprob, 'diagnostics')$C
    .psi <- .ruinprob$psi

    if (length(.psi.r) == 1) {
        .log.C <- log(.psi.C)
        .varu <- function(prob) {
            return(
                (.log.C - log(prob)) / .psi.r
            )
        }
    } else {
        .adjcoef <- min(Re(.psi.r[which(Im(.psi.r) == 0 & Re(.psi.r) > 0)]))
        .varu <- function(prob) {
            # Compute the Lundberg approximation (without multiplicative constant), multiplied by 4.0 (more
            # generally: a constant bigger than 1) to err on the safe side, and use this a upper bound for the
            # interval to be considered for the VaRu.
            .x.initial <- 4.0 * (-log(prob) / .adjcoef)
            .difffun <- function(.reserve, .prob) {
                .psi(.reserve) - .prob
            }
            .res <- mapply(FUN = uniroot,
                           upper = .x.initial,
                           .prob = prob,
                           MoreArgs = list(f = .difffun, lower = 0))
            return(unname(unlist(.res['root',])))
        }
    }

    return(.varu)
}
