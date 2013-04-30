\name{adjcoef}
\alias{adjcoef}

\title{Adjustment Coefficient}

\description{
    Returns the adjustment coefficient of a risk process with Gaussian
    diffusion.
}

\usage{
adjcoef(process)
}

\arguments{
  \item{process}{An object of the class \code{"riskproc"}.}
}

\details{
    The moment-generating function of the individual claim amounts is required
    to be contained within the \code{"claiminfo"} element of \code{process} (see
    \code{\link{claiminfo}}).  The only exception to this is the case of
    hypo-exponentially distributed claim amounts.
}

\value{
    The adjustment coefficient of the risk process represented by
    \code{process}.  If the adjustment coefficient does not exist or cannot be
    computed for any reason, \code{NULL} is returned instead.
}

\seealso{
    \code{\link{riskproc}} and \code{\link{claiminfo}} for more details on how
    to provide the information necessary to compute the adjustment coefficient.
}

\examples{
## Setting up a risk process with hypo-exponentially distributed claims
myprocess <- riskproc(
    claims   = claiminfo(hypoexp = list(rates = c(1, 10))),
    premium  = 2,
    freq     = 1,
    variance = 0.4
)

## Return the adjustment coefficient
adjcoef(myprocess)
}

\keyword{utilities}

% vim: ft=rhelp tw=80
