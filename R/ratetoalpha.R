ratetoalpha <- function(rate) {
    stopifnot(is.numeric(rate), all(!is.na(rate)))

    sapply(seq.int(along.with = rate),
           function(i) {
               prod(rate[-i] / (rate[-i] - rate[i]))
           }
    )
}
