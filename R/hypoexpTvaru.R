hypoexpTvaru <- function(process) {
    .psi <- get('psi', hypoexpRuinprob(process))

    .tvaru <- function(prob) {
        .quant <- hypoexpVaru(process)(prob)
        return(
            rev(cumsum(rev(int.multi(f = .psi, nodes = c(.quant, Inf))))) / prob + .quant
        )
    }

    return(.tvaru)
}
