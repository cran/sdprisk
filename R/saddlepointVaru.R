saddlepointVaru <- function(process, type = 2) {
    .adjcoef <- adjcoef(process)
    .KL.d1 <- get('KL.d1', process)
    .KL.d2 <- get('KL.d2', process)
    .zv <- get('zv', process)
    .vx <- get('vx', process)

    switch(type,
           { # TYPE 1 (switching back and forth between the monetary domain and the frequency domain)
            .phi <- function(x, prob, i) {
                if (i == 0) {
                    return(stats::qexp(p = prob, rate = .adjcoef, lower.tail = FALSE))
                } else {
                    v <- .vx(x, vmax = .adjcoef - .Machine$double.eps^(2/3))
                    return(x + (stats::qnorm(prob)^2 - .zv(v)^2) / (2 * v))
                }
            }
            .varu <- function(prob, n = 4) {
                stopifnot(n >= 0)
                for(i in 0:n) {
                    x <- .phi(x, prob, i)
                }
                return(x)
            }
        },
        { # TYPE 2 (iteration in the frequency domain)
            .phi <- function(x, prob, i) {
                if (i == 0) {
                    return(.adjcoef * (1 + 1 / log(prob)))
                } else {
                    if (i == 1) {
                        .div <- 2 * stats::qexp(p = prob, rate = .adjcoef, lower.tail = FALSE)
                    } else {
                        .div <- 2 * x * .KL.d2(x)
                    }
                    return(x + (stats::qnorm(prob)^2 - .zv(x)^2) / .div)
                }
            }
            .varu <- function(prob, n = 4) {
                stopifnot(n >= 0)
                for(i in 0:n) {
                    x <- .phi(x, prob, i)
                }
                return(structure(.KL.d1(x), saddlepoint = x))
            }
        }
    )

    return(.varu)
}
