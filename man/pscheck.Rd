\name{pscheck}
\docType{data}
\alias{pscheck}
\title{Example Dataset}
\description{An artificial data frame used to illustrate the techniques implemented in the package. Pseudorandom numbers were drawn from Weibull distributions of scale parameter 1 and shape parameter 3 (for failure times) and 5 (for censoring times). Each observed time is the minimum of the corresponding failure and censoring times, and its censoring status indicates a failure time smaller than a censoring time (see \code{format} below).}
\usage{pscheck}
\format{
The data frame contains 2 variables:
	\describe{
		\item{\code{t}}{
		Observed time	
		}
		\item{\code{delta}}{
		Censoring status (0 = censored, 1 = uncensored)
		}
	}
}
\author{Lopez de Ullibarri, I. and Jacome, M. A.}
