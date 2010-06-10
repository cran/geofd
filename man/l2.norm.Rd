\name{l2.norm}
\alias{l2.norm}
\alias{geofd.l2.norm}
\docType{data}
\title{
  Calculates L2 norm among functions
}
\description{
  Calculates the integral of the squared differences between functions
}
\usage{
  l2.norm(s, datafd, M)
}
\arguments{
  \item{s}{
    number of sites where the original dataset was measured
  }
  \item{datafd}{
    a functional data object representing a smoothed dataset. See DETAILS below.
  }
  \item{M}{
    symmetric matrix defining the roughness penalty for functions expressed in terms of a B-spline or Fourier basis. See DETAILS below.
  }
}
\details{
  In the details section of the \code{\link{okfd}} function can be reviewed the details of the paremeters used when calling this function.

  \bold{Functional data object}

  This is an object of the class \code{\link{fd}} it can be created using some functions like \code{\link{Data2fd}} or \code{\link{smooth.basis}}, take in count if a penalization parameter is going to be used.

  \bold{Roughness penalty matrix}

  This matrix is the output of one of the following functions: \code{\link{fourierpen}} y \code{\link{bsplinepen}}. The used function depends upon the smoothing type which is going to be applied.

}
\value{
  The calculated matrix of squared differences between each observation for each measured site. This matrix has two properties:
  \itemize{
  \item
  Is symmetric.
  \item
  It's diagonal is filled with zeros.
  }
}
\seealso{
\code{\link{okfd}} for doing ordinary kriging for function-value data, \code{\link{variog.okfd}} for functional empirical trace variogram calculation, \code{\link{variofit.okfd}} for fitting a variogram model in the funcional scenario.
}
\keyword{math}
