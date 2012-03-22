saddlepointSensitivity <- function(process, ...) {

    .p <- get('p', process)
    .zeta <- get('zeta', process)
    .KL.d2 <- get('KL.d2', process)
    .rv <- get('rv', process)
    .varu <- saddlepointVaru(process)
    .tvaru <- saddlepointTvaru(process, ...)

    .varu.sens <- function(prob) {
        .v <- attr(.varu(prob), 'saddlepoint')
        return(-sqrt(.KL.d2(.v)) / stats::dnorm(.rv(.v)))
    }
    .tvaru.sens <- function(prob) {
        return((.varu(prob) - .tvaru(prob)) / prob)
    }

    return(list(varu.sens  = .varu.sens,
                tvaru.sens = .tvaru.sens))
}
