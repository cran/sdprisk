\name{tvaru}
\alias{tvaru}
\alias{hypoexpTvaru}
\alias{saddlepointTvaru}

\title{Tail Value at Ruin}

\description{
    The tail value at ruin for a given probability level
    \eqn{\varepsilon}{\epsilon} is defined as the conditional expectation of the
    maximal aggregate loss given that it is above the value at ruin of level
    \eqn{\varepsilon}{\epsilon}.
}

\usage{
tvaru(process, method = c("saddlepoint", "hypoexp"), ...)
hypoexpTvaru(process)
saddlepointTvaru(process, type = c("tail", "density"), ...)
}

\arguments{
  \item{process}{a \code{"riskproc"} object.}
  \item{method}{character string indicating the calculation or approximation
      method.}
  \item{type}{character string indicating which function is to be used for the
      approximation.}
  \item{\dots}{further arguments that are passed on to
      \code{\link{saddlepointRuinprob}}, depending on the value of
      \code{method}.}
}

\details{
    \code{tvaru} is a wrapper function for \code{hypoexpTvaru} and
    \code{saddlepointTvaru}.

    \code{hypoexpTvaru} calculates the tail value at ruin in the case of
    hypo-exponentially distributed claim amounts by numerical integration of the
    probability of ruin, which can be computed exactly.

    \code{saddlepointTvaru} uses saddlepoint techniques for the approximation of
    the tail value at ruin.  More precisely, the saddlepoint approximation to
    the probability of is numerically integrated in the frequency domain, and
    implicitly also the saddlepoint approximation to the value at ruin (see
    \code{\link{varu}}) is used.  If \code{type = "tail"} the integrand is the
    probability of ruin (as function in the frequency domain), otherwise
    (\code{type = "density"}) it is essentially a re-scaled version of the
    probability of ruin due to claims.  The former requires fewer calculations
    and seems to produce slightly more accurate results.
}

\value{
    A function returning the tail value at ruin of a given probability level is
    returned.

    If \code{method = "saddlepoint"} or if \code{saddlepointTvaru} is used, the
    returned function has an additional second argument giving the number of
    iterations used for the approximation of the value at ruin (i. e., the lower
    integration limit), see \code{\link{varu}}.
}

\seealso{
    \code{\link{varu}}
}

\keyword{models}
\keyword{univar}

% vim: ft=rhelp tw=80
