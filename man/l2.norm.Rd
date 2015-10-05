\encoding{utf8}
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

  \bold{Roughness penalty matrix}

  This matrix is the output of one of the following functions: \code{\link{fourierpen}} y \code{\link{bsplinepen}}. The used function depends upon the smoothing type which is going to be applied.

  When the roughness penalty matrix is being calculated, the following considerations are taked in count:
  \itemize{
  \item The differential operator passed as parameter for both \code{fourierpen} and \code{bsplinepen} is always zero.
  \item When the selected smooth method is bsplines, the basis object passed to \code{\link{bsplinepen}} is the output of the function \code{\link{create.bspline.basis}} using 
\code{argvals} as the \code{rangeval} parameter, \code{nbasis} as the number of basis functions parameter and the default order of b-splines, which is four, a cubic spline, as the \code{norder} parameter.
  \item When the selected smooth method is fourier, the basis object is the output of the function \code{fourierpen}. The parameters \code{rangeval} and \code{nbasis} are the same as for \code{\link{create.bspline.basis}}, and the \code{period} parameter as the number of observations on each curve.
  }

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
\code{\link{okfd}} for doing Ordinary Kriging for function-value data, \code{\link{trace.variog}} for functional empirical trace variogram calculation, \code{\link{fit.tracevariog}} for fitting a variogram model in the functional scenario.
}
\keyword{math}
