saddlepointRuinprob <- function(process, jensen = FALSE, normalize = FALSE) {
    stopifnot(is.logical(jensen), is.logical(normalize))

    .p <- get('p', process)
    .q <- get('q', process)
    .zeta <- get('zeta', process)
    .claims <- get('claims', process)
    .mu <- mean(.claims)
    .adjcoef <- adjcoef(process)

    .KL <- get('KL', process)
    .KL.d1 <- get('KL.d1', process)
    .KL.d2 <- get('KL.d2', process)

    .vx <- get('vx', process)
    .rv <- get('rv', process)
    .sv <- get('sv', process)
    .zv <- get('zv', process)

    .epsilon <- .Machine$double.eps^(2/3)

    if (!jensen && normalize) {
        .corrconst <- integrate(
            f = function(.v) {
                # exp(.KL(.v) - .v * .KL.d1(.v)) * .KL.d2(.v) / sqrt(2 * pi * .KL.d2(.v))
                dnorm(.rv(.v)) * sqrt(.KL.d2(.v))
            },
            lower = -Inf,
            upper = .adjcoef - .epsilon
#             upper = optimize(
#                 function(.x){
#                     abs(.KL.d1(.x) - 1e+10)  
#                 },
#                 interval = c(-100000, 100)
#             )$minimum,
#             subdivisions = length(reserve)
        )$value
    } else {
        .corrconst <- 1
    }

    if (jensen) { ## Jensen (1992)
        .psi.v <- function(v) {
            return(pnorm(.zv(v), lower.tail = FALSE))
        }

        .psi.1.v <- function(v) {
            .rv.v <- .rv(v)
            #.sv.v <- .sv(v)
            .zv.v <- .zv(v)
            # return(pnorm(.zv.v, lower.tail = FALSE) - pnorm(.rv.v + log(.sv.v / (.rv.v * (1 - v / (.p * .zeta)))) / .rv.v, lower.tail = FALSE))
            return(pnorm(.zv.v - log(1 - v / .zeta / .p) / .rv.v) - pnorm(.zv.v))
        }

        .psi.2.v <- function(v) {
            .rv.v <- .rv(v)
            .sv.v <- .sv(v)
            .zv.v <- .zv(v)
            # return(pnorm(.rv.v + log(.sv.v / (.rv.v * (1 - .v / (.p * .zeta)))) / .rv.v,
            #              lower.tail = FALSE))
            return(pnorm(.zv.v - log(1 - v / .zeta / .p) / .rv.v, lower.tail = FALSE))
        }
    } else { # Lugannani and Rice (1980), Daniels (1954)
        .psi.v <- function(v) {
            .rv.v <- .rv(v)
            return(pnorm(.rv.v, lower.tail = FALSE) - dnorm(.rv.v) * (1 / .rv.v - 1 / .sv(v)))
        }

        .psi.1.v <- function(v) {
            return(dnorm(.rv(v)) * v / (.sv(v) * .p * .zeta * .corrconst))
        }

        .psi.2.v <- function(v) {
            .rv.v <- .rv(v)
            return(
                pnorm(.rv.v, lower.tail = FALSE)
                - dnorm(.rv.v) * (1 / .rv.v - 1 / .sv(v) * (1 - v / (.p * .zeta * .corrconst)))
            )
        }
    }

    .psi <- function(x) {
        ok <- which(x > .epsilon)
        .v <- .vx(x[ok], vmax = .adjcoef - .epsilon)
        res <- rep.int(1, length(x))
        res[ok] <- .psi.v(.v)
        return(res)
    }

    .psi.1 <- function(x) {
        ok <- which(x > .epsilon)
        .v <- .vx(x[ok], vmax = .adjcoef - .epsilon)
        if (jensen) {
            res <- rep.int(0, length(x))
        } else {
            res <- rep.int(1, length(x))
        }
        res[ok] <- .psi.1.v(.v)
        return(res)
    }

    .psi.2 <- function(x) {
        ok <- which(x > .epsilon)
        .v <- .vx(x[ok], vmax = .adjcoef - .epsilon)
        if (jensen) {
            res <- rep.int(1, length(x))
        } else {
            res <- rep.int(0, length(x))
        }
        res[ok] <- .psi.2.v(.v)
        return(res)
    }

    return(
        structure(
            list(psi = .psi, psi.1 = .psi.1, psi.2 = .psi.2),
            compmethod  = 'saddlepoint',
            riskproc    = process,
            parameters  = list(jensen = jensen,
                               normalize = normalize),
            diagnostics = list(psi.v = .psi.v,
                               psi.1.v = .psi.1.v,
                               psi.2.v = .psi.2.v,
                               corrconst = .corrconst)
        )
    )
}
