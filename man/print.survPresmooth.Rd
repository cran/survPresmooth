\name{print.survPresmooth}
\alias{print.survPresmooth}
\title{Print a survPresmooth Object}
\description{Print method for `survPresmooth' objects.}
\usage{\method{print}{survPresmooth}(x, long = FALSE, more = NULL, \dots)}
\arguments{
	\item{x}{
An object of class `survPresmooth'.
	}
	\item{long}{
A logical value. If TRUE the entire estimated curve is printed; if FALSE, the default, only the head of the estimated curve is printed if it consists of more than 100 values.
	}
	\item{more}{
A character vector specifying the names of additional components of the \code{x} object that will be printed. The valid values for the components of the character vector are "data", "q.weight", "mise", "grid.pil", "pilot.bw", "grid.bw" and "p.hat". The default is NULL.
	}
	\item{\dots}{
Optional arguments for the default method (i.e. \code{print.default}) of the \code{print} generic function.
	}
}
\author{Lopez de Ullibarri, I. and Jacome, M. A.}
\seealso{\code{\link{presmooth}}}
\examples{
	\dontrun{
## Printing an object of class survPresmooth with 4 digits, including in the output the
# value of the pilot bandwidth and the estimated p function

print(presmooth(t, delta, pscheck), more = c("pilot.bw", "p.hat"), digits = 4)
	}
}

