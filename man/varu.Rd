\name{varu}
\alias{varu}
\alias{hypoexpVaru}
\alias{saddlepointVaru}

\title{Value at Ruin}

\description{
    The value at ruin at a given probability level \eqn{\varepsilon}{\epsilon}
    is defined as the minimal capital that is required in order to have a ruin
    probability of at most \eqn{\varepsilon}{\epsilon}.  This is equivalent to
    the \eqn{(1-\varepsilon)}{(1-\epsilon)}-quantile of the maximal aggregate
    loss.
}

\usage{
varu(process, method = c("saddlepoint", "hypoexp"), ...)
hypoexpVaru(process)
saddlepointVaru(process, type = 2)
}

\arguments{
  \item{process}{a \code{"riskproc"} object.}
  \item{method}{character string indicating the calculation or approximation
      method.}
  \item{type}{number indicating the type of approximation; possible choices are
      1 and 2.}
  \item{\dots}{further arguments that are passed on to \code{saddlepointVaru},
      depending on the value of \code{method}.}
}

\details{
    \code{varu} is a wrapper function for \code{hypoexpVaru} and
    \code{saddlepointVaru}.

    \code{hypoexpVaru} calculates the value at ruin in the case of
    hypo-exponentially distributed claim amounts by numerical inversion of the
    probability of ruin, which can be computed exactly.

    \code{saddlepointVaru} uses saddlepoint techniques for the approximation of
    the value at ruin, more specifically, the inversion algorithms provided by
    \cite{Wang (1995)}.  The first one (\code{type = 1}) is only given for
    completeness (or comparison purposes), because, due to repeatedly switching
    back and forth between the monetary domain the frequency (saddlepoint)
    domain, it is much slower than the second one (\code{type = 2}), which is
    performed entirely in the frequency domain.  Refer to the references given
    below for more details.
}

\value{
    A function returning the value at ruin of a given probability level is
    returned.

    If \code{method = "saddlepoint"} or if \code{saddlepointVaru} is used, the
    returned function has an additional second argument giving the number of
    iterations.
}

\references{
    Wang, Suojin (1995) One-Step Saddlepoint Approximations for Quantiles.
    \emph{Computational Statistics and Data Analysis} \bold{20}(1), pp. 65--74.
}

\keyword{models}
\keyword{univar}

% vim: ft=rhelp tw=80
