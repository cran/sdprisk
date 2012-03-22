mgfhypoexp <- function(x, rate = 1, difforder = 0) {
    stopifnot(difforder %in% c(0, 1, 2))
    inv.diff <- 1 / outer(rate, x, '-')
    prod(rate) * apply(inv.diff, 2, prod) * switch(difforder + 1,
        1,
        colSums(inv.diff),
        colSums(inv.diff)^2 + colSums(inv.diff^2)
    )
}
