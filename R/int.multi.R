int.multi <- function(f, nodes, ...) {
    vapply(X         = mapply(FUN      = stats::integrate,
                              lower    = utils::head(nodes, -1L),
                              upper    = utils::tail(nodes, -1L),
                              MoreArgs = list(f = match.fun(f), ...),
                              SIMPLIFY = FALSE),
           FUN       = `[[`,
           FUN.VALUE = numeric(1L),
           index     = 'value')
}
