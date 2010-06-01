\name{control.presmooth}
\Rdversion{1.1}
\alias{control.presmooth}
\title{
Control Values for presmooth() Function
}
\description{
This function returns a list of values for control parameters of the \code{presmooth} function.
}
\usage{
control.presmooth(n.boot = 500, q.weight = c(0.2, 0.8), 
	length.grid.bw.plugin = 100, length.grid.ise = 100,
	pilot.par.ini= NULL, save.data = FALSE, save.mise = FALSE,
	na.action = na.omit)
}
\arguments{
  \item{n.boot}{
A numeric value specifying the number of bootstrap resamples used in bootstrap bandwidth selection. The default is 500.
}
  \item{q.weight}{
A numeric vector of length 2 specifying the order of quantiles of the observed times that determine the left- and right-ends of the support of the weight function. The default value is c(0.2, 0.8).
}
  \item{length.grid.bw.plugin}{
An integer specifying the number of points of the grid used for numerical integration (Simpson's rule) of integrals involved in plug-in bandwidth selection. The default is 100.
}
  \item{length.grid.ise}{
An integer specifying the number of points of the grid used for computing integrated squared error by numerical integration (trapezium rule) in bootstrap bandwidth selection. The default is 100.
}
  \item{pilot.par.ini}{
A numeric vector of length 6 specifying the initial values of optimization routines used for pilot bandwidth computation in plug-in and bootstrap bandwidth selection (although, in the last case, it is not needed for survival and cumulative hazard function estimation). If NULL, the default, values are internally computed.
}
  \item{save.data}{
A logical value specifying if the data must be saved as a component of the value returned by the \code{presmooth} function. The default is FALSE.
}
  \item{save.mise}{
A logical value specifying if the MISE computed for bootstrap bandwidth selection must be saved as a component of the value returned by the \code{presmooth} function. The default is FALSE.
}
  \item{na.action}{
A function specifying how to handle missing values. The default value (na.omit), causes the deletion of the rows containing missing values in the (internal) data frame formed by the \code{times} and \code{status} arguments of the call to the \code{presmo}-\code{oth} function.
}
}
\value{
A list whose components are the arguments of the function, its defaults being replaced with the values the function was called with.
}
\author{
Lopez-de-Ullibarri, I. and Jacome, M. A.
}
\seealso{
\code{\link{presmooth}}
}

