#' p-quantile from the hypo-exponential distribution
#' 
#' The p-quantile from the hypo-exponential distribution is computed implicitly by computing the root of H(x) = F(x) - p,
#' where F is the distribution function of the hypo-exponential distribution
#' 
#' @param p vector of probabilities
#' @param rate vector of (unique) rates 
#' @param interval search interval in which \code{uniroot} searches for the quantile
#' @return p-quantile from the hypo-exponential distribution
#' @author Sebastian Szugat \email{szugat@@statistik.tu-dortmund.de}
qhypoexp <- Vectorize(function(p, rate, interval = c(0, 10^10)){
  require(sdprisk)
  stopifnot(p >= 0 && p <= 1,
            all(is.finite(rate)))
  H <- function(b){
    phypoexp(q = b, rate = rate) - p
  }
  uniroot(H, interval = interval)$root
}, vectorize.args = "p")