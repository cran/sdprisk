\name{hypoexp}
\alias{phypoexp}
\alias{dhypoexp}
\alias{rhypoexp}
\alias{mgfhypoexp}

\title{Hypo-Exponential Distribution}

\description{
    Density, distribution function, random generation and moment-generating
    function (and its first two derivatives) for the hypo-exponential
    distribution with rates \code{rate}.
}

\usage{
dhypoexp(x, rate = 1, log = FALSE)
phypoexp(q, rate = 1, lower.tail = TRUE, log.p = FALSE, tailarea = FALSE)
rhypoexp(n = 1, rate = 1)
mgfhypoexp(x, rate = 1, difforder = 0)
}

\arguments{
    \item{x, q}{vector of quantiles.}
    \item{n}{number of observations. If \code{length(n) > 1}, the length is
        taken to be the number required.}
    \item{difforder}{the order of derivative for the moment-generating function;
        currently only implemented for 0, 1, 2.}
    \item{rate}{vector of (unique) rates.}
    \item{lower.tail}{logical; if \code{TRUE}, probabilities are
        \eqn{\mathbf{P}(X \le x)}{P(X \le x)}, otherwise \eqn{\mathbf{P}(X >
            x)}{P(X > x)}.}
    \item{log, log.p}{logical; if \code{TRUE}, probabilities \eqn{p} are given
        as \eqn{\log(p)}{log(p)}.}
    \item{tailarea}{logical; if \code{TRUE}, probabilities are given for the
        integrated tail area distribution.}
}

\details{
    The sum of \eqn{n} independent exponentially distributed random variables
    \eqn{X_{i}} with rate parameters \eqn{\lambda_{i}} has a hypo-exponential
    distribution with rate vector \eqn{(\lambda_{1}, \dots, \lambda_{n})}.

    The hypo-exponential distribution is a generalization of the Erlang
    distribution (a Gamma distribution with an integer-valued shape parameter)
    and a special case of the phase-type distribution (see References section).
}

\value{
    \code{dhypoexp} gives the density, \code{phypoexp} gives the distribution
    function (or the integrated tail area distribution function),
    \code{rhypoexp} generates random deviates and \code{mgfhypoexp} gives the
    moment-generating function (or its derivative up to the second order)
}

\references{
    Neuts, M. F. (1981)  \emph{Matrix-Geometric Solutions in Stochastic Models:
    An Algorithmic Approach}, reprinted and corrected.
}

\note{
    If \code{length(rate) == 1}, \code{dhypoexp}, \code{phypoexp} and
    \code{rhypoexp} are equivalent to \code{\link[stats]{dexp}},
    \code{\link[stats]{pexp}} and \code{\link[stats]{rexp}} with rate parameter
    \code{rate} and should, in fact, be replaced by the latter ones for
    computation speed.
}

\seealso{
    \code{\link[stats]{dexp}}, \code{\link[stats]{dgamma}}
}

\examples{
## Random generation
rhypoexp(10, c(3, 5))

## Mean
mu <- mgfhypoexp(0, c(3, 5), difforder = 1)

## Variance
mgfhypoexp(0, c(3, 5), difforder = 2) - mu^2
}

\keyword{datagen}
\keyword{distribution}
\keyword{univar}

% vim: ft=rhelp tw=80
