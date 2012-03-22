riskproc <- function(claims, premium, freq, variance) {
    stopifnot(is.numeric(premium), premium > 0, is.numeric(freq), freq > 0,
              is.numeric(variance), variance >= 0, "mu" %in% names(claims),
              is.claiminfo(claims))

    .mu <- mean(claims)
    .p <- 1 - freq * .mu / premium
    .q <- freq * .mu / premium
    .beta <- premium / (freq * .mu) - 1
    .zeta <- 2 * premium / variance

    try({
        mgf    <- claims$mgf
        mgf.d1 <- claims$mgf.d1
        mgf.d2 <- claims$mgf.d2
    }, silent = TRUE)

    if (all(exists('mgf')    && is.function(mgf),
            exists('mgf.d1') && is.function(mgf.d1),
            exists('mgf.d2') && is.function(mgf.d2))) {

        KL <- function(x) {
#             if (!is.na(.adjcoef)) {
#                 x[which(x > .adjcoef)] <- NA   
#             }
            val <- .p * x / (x - x^2 / .zeta + .q / .mu * (1 - mgf(x)))
            val[which(x == 0)] <- 1
            return(log(val))
        }

        KL.d1 <- function(x) {
#             if (!is.na(.adjcoef)) {
#                 x[which(x > .adjcoef)] <- NA   
#             }

            mgf.x <- mgf(x)
            mgf.d1.x <- mgf.d1(x)

            res <- (x^2 * .mu + .q * .zeta - .q * .zeta * mgf.x + x * mgf.d1.x * .q
                    * .zeta) / (x * .zeta * .mu - x^2 * .mu + .q * .zeta - .q * .zeta *
                    mgf.x) / x

            mynan <- which(is.infinite(res))
            res[mynan] <- sapply(x[mynan], numDeriv::grad, func = KL)
            return(unlist(res))
        }

        KL.d2 <- function(x) {
#             if (!is.na(.adjcoef)) {
#                 x[which(x > .adjcoef)] <- NA   
#             }

            mgf.x <- mgf(x)
            mgf.d1.x <- mgf.d1(x)
            mgf.d2.x <- mgf.d2(x)

            return(-(-2 * x * .zeta^2 * .mu * .q * mgf.x + 2 * x^2 * .zeta^2 * .mu *
                     mgf.d1.x * .q + 4 * x^2 * .mu * .q * .zeta * mgf.x - 4 * x^3 *
                     .mu * mgf.d1.x * .q * .zeta - 2 * .q^2 * .zeta^2 * mgf.x + .q^2 *
                     .zeta^2 * mgf.x^2 - mgf.d2.x * .q * .zeta^2 * x^3 * .mu +
                     mgf.d2.x * .q * .zeta * x^4 * .mu + x^2 * mgf.d2.x * .q^2 *
                     .zeta^2 * mgf.x + 2 * x * .zeta^2 * .mu * .q - 4 * x^2 * .mu * .q
                     * .zeta - x^2 * mgf.d1.x^2 * .q^2 * .zeta^2 - x^2 * mgf.d2.x *
                     .q^2 * .zeta^2 - x^4 * .mu^2 + .q^2 * .zeta^2) / ( - x * .zeta *
                     .mu + x^2 * .mu - .q * .zeta + .q * .zeta * mgf.x)^2 / x^2)
        }

    } else {
        KL    <- NULL
        KL.d1 <- NULL
        KL.d2 <- NULL
    }

    if (!(is.null(KL) | is.null(KL.d1) | is.null(KL.d2))) {

        vx <- function(x, vmin = -1e+64, vmax) {
            sapply(x,
                   function(.x) {
                       stats::uniroot(function(.y){KL.d1(.y) - .x},
                                      interval = c(vmin, vmax))$root
                   }) 
        }

        rv <- function(v) {
            return(sign(v) * sqrt(2 * (v * KL.d1(v) - KL(v))))
        }

        sv <- function(v) {
            return(v * sqrt(KL.d2(v)))
        }

        zv <- function(v) {
            rv.v <- rv(v)
            return(rv.v + log(sv(v) / rv.v) / rv.v)
        }

    }

    structure(
        list(premium   = premium,
             freq      = freq,
             variance  = variance,
             diffusion = (variance != 0),
             p         = .p,
             q         = .q,
             beta      = .beta,
             zeta      = .zeta,
             KL        = KL,
             KL.d1     = KL.d1,
             KL.d2     = KL.d2,
             vx        = vx,
             rv        = rv,
             sv        = sv,
             zv        = zv,
             claims    = claims),
        class = c('riskproc', 'list')
    )
}
