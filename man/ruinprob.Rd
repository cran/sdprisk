\name{ruinprob}
\alias{ruinprob}
\alias{boundsRuinprob}
\alias{fftRuinprob}
\alias{hypoexpRuinprob}
\alias{saddlepointRuinprob}

\title{Calculation or Approximation of the Probability of Ruin}

\description{
    This functions provide various approximation methods for the (total)
    probability of ruin, the probability of ruin due to oscillation and the
    probability of ruin due to a claim.  Exact calculations are possible in the
    case of hypo-exponentially distrubuted claim amounts.
}

\usage{
ruinprob(process, method = c("saddlepoint", "fft", "bounds", "hypoexp", "lundberg"), \dots)
boundsRuinprob(process, interval, maxreserve, richardson = TRUE, use.splines = FALSE)
fftRuinprob(process, interval, maxreserve, n, use.splines = FALSE)
hypoexpRuinprob(process)
saddlepointRuinprob(process, jensen = FALSE, normalize = FALSE)
}

\arguments{
  \item{process}{a \code{"riskproc"} object.}
  \item{method}{character string indicating the method used for approximation or
      calculation.}
  \item{interval}{interval width for the discretization of the claim
      distribution.}
  \item{maxreserve}{maximal value of the initial reserve for which the
      approximation can be calculated.}
  \item{n}{Length of the probability vectors resulting from the discretization.}
  \item{richardson}{logical; if \code{TRUE}, Richardson extrapolation is used
      for the approximation of the probability of ruin due to oscillation.}
  \item{use.splines}{logical; if \code{TRUE}, a cubic spline interpolation is
      used instead of step functions.}
  \item{jensen}{logical; if \code{TRUE}, the formulae of \cite{Jensen (1992)}
      are used instead of the ones by \cite{Lugannani and Rice (1980)} and
      \cite{Daniels (1954)} (see references).}
  \item{normalize}{logical; if \code{TRUE}, the saddlepoint approximation to the
      density is re-normalized such that it integrates to 1; only used if
      \code{jensen} is \code{FALSE}.}
  \item{\dots}{further arguments that are passed on to \code{boundsRuinprob},
      \code{fftRuinprob}, \code{hypoexpRuinprob} or \code{saddlepointRuinprob},
      depending on the value of \code{method}.}
}

\details{
    \code{ruinprob} is a wrapper function for the other ones given here.
}

\value{
    \item{psi}{the total probability of ruin (as a function of the initial
        reserve).}
    \item{psi.1}{the probability of ruin due to oscillation (as a function of
        the initial reserve).}
    \item{psi.2}{the probability of ruin due to a claim (as a function of the
        initial reserve).}
    \item{\dots}{}
}

\references{
    Daniels, H. E. (1954)  Saddlepoint Approximations in Statistics.
    \emph{Annals of Mathematical Statistics} \bold{25}(4), pp. 631--650.

    Gatto, R. and Mosimann, M. (2012)  Four Approaches to Compute the
    Probability of Ruin in the Compound Poisson Risk Process with Diffusion.
    \emph{Mathematical and Computer Modelling} \bold{55}(3--4), pp. 1169--1185

    Jensen, J. L. (1992)  The Modified Signed Likelihood Statistic and
    Saddlepoint Approximations.  \emph{Biometrika} \bold{79}(4), pp. 693--703.

    Lugannani, R. and Rice, S. (1980)  Saddle Point Approximation for the
    Distribution of the Sum of Independent Random Variables.  \emph{Advances in
    Applied Probability} \bold{12}(2), pp. 475--490.
}

\seealso{
    \code{\link{riskproc}}, \code{\link{claiminfo}}
}

\keyword{models}
\keyword{univar}

% vim: ft=rhelp tw=80
