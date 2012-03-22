\name{riskproc}
\alias{riskproc}
\alias{is.riskproc}

\title{Compound Poisson Risk Process with Diffusion}

\description{
    Creates an \R object representing a compound Poisson risk process with
    Gaussian diffusion, which contains some or all information necessary for
    further processing.
}

\usage{
riskproc(claims, premium, freq, variance)
is.riskproc(x)
}

\arguments{
  \item{claims}{a \link{claiminfo} object.}
  \item{premium}{premium rate.}
  \item{freq}{claim frequency.}
  \item{variance}{squared volatility of the Wiener component; currently only
      implemented for \code{variance > 0}.}
  \item{x}{an \R object.}
}

\details{
    Given the arguments, most prominently \code{claims}, various auxialiary parameters and
    functions associated with the risk process to be represented are calculated.
}

\value{
    \code{riskproc} returns an object of the class \code{"riskproc"}.
    Internally, this is a list containing various elements (depending on the
    information provided in the arguments).

    \code{is.riskproc} returns \code{TRUE} if \code{x} is a \code{"riskproc"}
    object, and \code{FALSE} otherwise.
}

\seealso{
    \code{\link{claiminfo}} for more details about passing on information about
    the distribution of the individual claim amounts.
}

\examples{
## A risk process with hypo-exponentially distributed individual claim amounts
riskproc(
    claims   = claiminfo(hypoexp = list(rates = c(1, 10))),
    premium  = 2,
    freq     = 1,
    variance = 0.4
)
}

\keyword{models}
\keyword{utilities}

% vim: ft=rhelp tw=80
