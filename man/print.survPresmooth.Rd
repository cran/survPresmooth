\name{print.survPresmooth}
\Rdversion{1.1}
\alias{print.survPresmooth}
\title{
Print a survPresmooth Object
}
\description{
Print method for \code{survPresmooth} objects.
}
\usage{
\method{print}{survPresmooth}(x, more = NULL, ...)
}
\arguments{
  \item{x}{
An object of class 'survPresmooth'.
}
  \item{more}{
A character vector specifying the names of additional components of the \code{x} object that will be printed. The default is  NULL.
}
   \item{\dots}{
Optional arguments for the default method (i.e. \code{print.default}) of the \code{print} generic function.
}
}
\author{
Lopez-de-Ullibarri, I. and Jacome, M. A.
}
\seealso{
\code{\link{presmooth}}
}

\examples{
## Load example dataset

data(pscheck)

## Printing an object of class \code{survPremooth} with 4 digits, including
# in the output the value of the pilot bandwidth

print(presmooth(t, delta, pscheck), more = c("pilot.bw", "p.hat"), digits = 4)
}
