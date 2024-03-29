\name{control.presmooth}

\alias{control.presmooth}

\title{Control Values for the presmooth() Function}

\description{This function returns a list of values for control
  parameters of the \code{presmooth} function.
}

\usage{control.presmooth(n.boot = c(5000, 1000), q.weight = c(0.2, 0.8),
   k = 1, length.grid.bw.plugin = 100, length.grid.ise = 100,
   pilot.par.ini = NULL, save.data = FALSE, save.mise = FALSE,
   na.action = na.omit)
}

\arguments{
  \item{n.boot}{A numeric vector of length 2 or, alternatively, a
  numeric value specifying the number of bootstrap resamples used in
  bootstrap bandwidth selection. If a vector, the first element is used
  for S and H estimation, and the second for f and h. The default is
  c(5000, 1000).
  }
  \item{q.weight}{A numeric vector of length 2 specifying the order of
  quantiles of the observed times that determine the left- and
  right-ends of the support of the weight function. The default is
  c(0.2, 0.8).
}
  \item{k}{A numeric value specifying the multiple of the data range
  used as the maximum possible value for the selected plug-in or
  bootstrap bandwidth. The default is 1.
  }
  \item{length.grid.bw.plugin}{An integer specifying the number of
  points of the grid used for numerical integration (Simpson's rule) of
  integrals involved in plug-in bandwidth selection. The default is 100.
  }
  \item{length.grid.ise}{An integer specifying the number of points of
  the grid used for computing the integrated squared error by numerical
  integration (Simpson's rule) in bootstrap bandwidth selection. The
  default is 100.
  }
  \item{pilot.par.ini}{A numeric vector of length 8 specifying the
  initial values of optimization routines used for pilot bandwidth
  computation in plug-in and bootstrap bandwidth selection (although, in
  the last case, it is not needed for survival and cumulative hazard
  function estimation). If \code{NULL}, the default, values are
  internally computed.
  }
  \item{save.data}{A logical value specifying if the data must be saved
  as a component of the value returned by the \code{presmooth}
  function. The default is \code{FALSE}.
  }
  \item{save.mise}{A logical value specifying if the MISE computed for
  bootstrap bandwidth selection must be saved as a component of the
  value returned by the \code{presmooth} function. The default is
  \code{FALSE}.
}
  \item{na.action}{A function specifying how to handle missing
  values. The default value, na.omit, deletes the rows containing
  missing values in the (internal) data frame formed by the \code{times}
  and \code{status} arguments of the call to the \code{presmooth}
  function.
}
}

\value{A list whose components are the arguments of the function, its
  defaults being replaced with the values the function was called with.
}

\details{The default values of \code{n.boot} represent a trade-off
  between computational speed and precision in bootstrap bandwidth
  selection. If enough computing power or time are available, it is
  recommended to increase the number of bootstrap resamples up to 10000
  (which could be very slow for hazard and density estimation).

  The six first values of \code{pilot.par.ini} are three pairs of
  parameters specifying three Weibull distributions (in each pair, the
  first element is the shape parameter; the second one, the scale
  parameter). These three Weibull distributions form a mixture, the
  weights of its first two components being given by the last two values
  of \code{pilot.par.ini}.
}

\author{Ignacio Lopez-de-Ullibarri [aut, cre],
  Maria Amalia Jacome [aut]
}

\seealso{\code{\link{presmooth}}}
