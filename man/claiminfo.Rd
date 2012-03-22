\name{claiminfo}
\alias{claiminfo}
\alias{is.claiminfo}
\alias{is.hypoexp}

\title{Distribution Information about Individual Claim Amounts}

\description{
    Creates or tests for claim information objects.
}

\usage{
claiminfo(\dots)
is.claiminfo(x)
is.hypoexp(x)
}

\arguments{
    \item{x}{An \R object}
    \item{\dots}{various objects determining the individual claim amount
        distribution. Refer to the details section.}
}

\details{
    Typicall usages are:
\preformatted{
claiminfo(mgf, mgf.d1, mgf.d2, pdf, cdf, mean)
}

To be completed.
}

\value{
    \code{claiminfo} returns an object of the class \code{"claiminfo"} (see
    details section).

    \code{is.claiminfo} returns \code{TRUE} if \code{x} is a \code{"claiminfo"}
    object, and \code{FALSE} otherwise.

    \code{is.hypoexp} returns \code{TRUE} if \code{x} is a \code{"claiminfo"}
    object describing hypo-exponentially distributed individual claim amounts
    (see \link{dhypoexp}), and \code{FALSE} otherwise.  If \code{x} is a
    \code{\link{riskproc}} object, the function is applied to the
    \code{"claiminfo"} object contained within it.
}

\note{
    To be completed.
}

\seealso{
    \code{\link{riskproc}}
}

\examples{
## For hypo-exponentially distribution claim amounts
claiminfo(hypoexp = list(rates = c(1, 10)))

## A more complicated example
\dontrun{claiminfo()}
}

\keyword{distribution}
\keyword{models}
\keyword{utilities}

% vim: ft=rhelp tw=80
