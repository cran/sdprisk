int.multi <- function(f, nodes, ...) {
    unname(unlist(mapply(
        FUN      = stats::integrate,
        lower    = utils::head(nodes, -1),
        upper    = utils::tail(nodes, -1),
        MoreArgs = list(f = f, ...)
    )[1, ]))
}
