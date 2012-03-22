hypoexpRuinprob <- function(process) {
    stopifnot(is.hypoexp(process))

    .p <- get('p', process)
    .q <- get('q', process)
    .zeta <- get('zeta', process)
    .claims <- get('claims', process)

    .coef <- get('hypoexp', .claims)$coef
    .rates <- get('hypoexp', .claims)$rates
    .mu <- mean(.claims)

    mypoly.factors <- polynom::as.polylist(
        lapply(.rates, function(arg) { c(arg, -1) })
    )

    mypoly.rhs <- .mu * polynom::polynomial(c(.zeta, -1)) * prod(mypoly.factors)
    mypoly.lhs <- .zeta * .q * sum(polynom::as.polylist(lapply(
        seq.int(along.with = mypoly.factors),
        function(.index) {
            .coef[.index] * prod(mypoly.factors[-.index])
        }
    )))

    .r <- polynom:::solve.polynomial(mypoly.lhs - mypoly.rhs)

    .const <- solve(
        a = rbind(outer(.rates, .r, function(..rates, ..r) { ..rates / (..rates - ..r) }),
                  rep.int(1, length(.r))),
        b = rep.int(1, length(.r))
    )

    .const1 <- .r * .const / (.p * .zeta)
    .const2 <- .const - .const1

    psi <- function(x) {
        pmin.int(1, Re(drop(exp(outer(x, -.r)) %*% .const)))
    }

    psi.1 <- function(x) {
        pmin.int(1, Re(drop(exp(outer(x, -.r)) %*% .const1)))
    }

    psi.2 <- function(x) {
        pmin.int(1, Re(drop(exp(outer(x, -.r)) %*% .const2)))
    }

    return(
        structure(
            list(psi = psi, psi.1 = psi.1, psi.2 = psi.2),
            compmethod  = 'hypoexp',
            riskproc    = process,
            parameters  = list(NULL),
            diagnostics = list(C = .const, C1 = .const1, C2 = .const2, r = .r)
        )
    )
}
