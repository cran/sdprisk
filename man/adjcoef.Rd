\name{adjcoef}
\alias{adjcoef}

\title{Adjustment Coefficient}

\description{
    Computes the adjustment coefficient of a risk process with Gaussian
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
    \code{process}.

    If \code{process} already contains the adjustment coefficient -- regardless
    of whether it has originally been supplied at the time of construction or of
    whether it stems from a previous run --, then this value is returned without
    repeating the computation.  Re-computation can be forced by manually
    removing the list element \code{"adjcoef"} of \code{process} or setting it
    to \code{NULL}.
}

\note{
    As a side effect, the object that has been passed as \code{process} argument
    is modified in the global environment such that it contains the adjustment
    coefficient (as a list element named \code{"adjcoef"}).  The rationale
    behind this is avoiding repeated calculations, because other functions, such
    as \code{\link{tvaru}}, may require the adjustment coefficient.
}

\seealso{
    \code{\link{riskproc}} and \code{\link{claiminfo}} for more details on how
    to provide the information necessary to compute the adjustment coefficient.
}

\examples{
## Setting up a risk process with hypo-exponentially distrubuted claims
myprocess <- riskproc(
    claims   = claiminfo(hypoexp = list(rates = c(1, 10))),
    premium  = 2,
    freq     = 1,
    variance = 0.4
)

## Compute the adjustment coefficient
adjcoef(myprocess)

## Note that the computed information is now stored inside the process object
myprocess$adjcoef
}

\keyword{utilities}

% vim: ft=rhelp tw=80
