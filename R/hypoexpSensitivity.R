hypoexpSensitivity <- function(process) {
    .p <- get('p', process)
    .zeta <- get('zeta', process)
    .psi.1 <- get('psi.1', hypoexpRuinprob(process))
    .varu <- hypoexpVaru(process)
    .tvaru <- hypoexpTvaru(process)

    .varu.sens <- function(prob) {
        return(-1 / (.p * .zeta * .psi.1(.varu(prob))))
    }

    .tvaru.sens <- function(prob) {
        return((.varu(prob) - .tvaru(prob)) / prob)
    }

    return(list(varu.sens  = .varu.sens,
                tvaru.sens = .tvaru.sens))
}
