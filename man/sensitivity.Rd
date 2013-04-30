\name{sensitivity}
\alias{sensitivity}
\alias{hypoexpSensitivity}
\alias{saddlepointSensitivity}

\title{Sensitivity of the Value and Tail Value at Ruin}

\description{
    The sensitivities of both the value and the tail value at ruin are defined
    as their respective derivatives with respect to the probability level.
}

\usage{
sensitivity(process, method = c("saddlepoint", "hypoexp"), ...)
hypoexpSensitivity(process)
saddlepointSensitivity(process, ...)
}

\arguments{
    \item{process}{a \link{riskproc} object.}
    \item{method}{character string indicating the calculation or approximation
        method.}
    \item{\dots}{further arguments that are passed on to
        \code{\link{saddlepointTvaru}}.}
}

\value{
    \item{varu}{a function returning the sensitivity of the value at ruin.}
    \item{tvaru}{a function returning the sensitivity of the tail value at
        ruin.}
}

\keyword{models}
\keyword{univar}

% vim: ft=rhelp tw=80
