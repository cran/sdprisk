# is.hypoexp <- function(x) {
#     if (is.riskproc(x)) {
#         Recall(get('claims', x))
#     } else {
#         stopifnot(is.claiminfo(x))
#         return('hypoexp' %in% names(x))
#     }
# }

is.hypoexp <- function(x) {
    UseMethod('is.hypoexp')
}

is.hypoexp.default <- function(x) {
    'hypoexp' %in% names(as.list(x))
}

is.hypoexp.claiminfo <- function(x) {
    'hypoexp' %in% names(x)
}

is.hypoexp.riskproc <- function(x) {
    is.hypoexp(get('claims', x))
}
