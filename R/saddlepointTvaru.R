saddlepointTvaru <- function(process, type = c('tail', 'density'), ...) {
    type <- match.arg(type)

    .q     <- get('q', process)
    .p     <- get('p', process)
    .zeta  <- get('zeta', process)
    .KL.d1 <- get('KL.d1', process)
    .KL.d2 <- get('KL.d2', process)

    .psi     <- attr(saddlepointRuinprob(process, ...), 'diagnostics')
    .adjcoef <- adjcoef(process)
    .varu    <- saddlepointVaru(process, type = 2)

    .tvaru <- switch(type,
        tail = function(prob, n = 4) {
            .quant <- .varu(prob, n)
            .intmin <- attr(.quant, 'saddlepoint')
            .quant <- as.vector(.quant)
            .aux <- lower.tri(x = matrix(1, .temp <- length(prob), .temp), diag = TRUE)

            drop(
                int.multi(
                    f = function(.v) {
                        .psi$psi.v(.v) * .KL.d2(.v)
                    },
                    nodes = c(.intmin, .adjcoef)
                ) %*% .aux
            ) / prob + .quant
        },
        density = function(prob, n = 4) {
            .quant <- .varu(prob, n)
            .intmin <- attr(.quant, 'saddlepoint')
            .quant <- as.vector(.quant)
            .aux <- lower.tri(x = matrix(1, .temp <- length(prob), .temp), diag = TRUE)

            drop(
                int.multi(
                    f = function(.v) {
                        .psi$psi.1.v(.v) * .KL.d1(.v) * .KL.d2(.v)
                    },
                    nodes = c(.intmin, .adjcoef)
                ) %*% .aux * .p * .zeta
            ) / prob
        }
    )

    return(.tvaru)
}
