\encoding{utf8}
\name{variofit.okfd}
\alias{variofit.okfd}
\alias{geofd.variofit.okfd}
\title{
  Variogram Estimation
}
\description{
  Fits a parametric model to a empirical variogram and estimates covariance parameters. Aditionally all fitted variogram models are plotted for verification purpose.
}
\usage{
variofit.okfd(emp.trace.vari, models, sigma2.0, phi.0, fix.nugget, nugget,
            fix.kappa, kappa, max.dist.variogram)
}
\arguments{
  \item{emp.trace.vari}{
    empirical trace-variogram. An object returned from the \code{\link{variog}} function.
}
  \item{models}{
    vector of correlation functions against which empirical trace variogram will be fitted.
}
  \item{sigma2.0}{
    initial value for the covariance parameter \eqn{\sigma^2}{sigma^2} (partial sill). For further details see documentation for the parameter \code{ini.cov.pars} from the \code{\link{variofit}} function.
}
  \item{phi.0}{
    initial value for the covariance parameter \eqn{\phi}{phi} (range).  For further details see documentation for the parameter \code{ini.cov.pars} from the \code{\link{variofit}} function.
}
  \item{fix.nugget}{
    logical, indicating whether the \code{nugget} parameter should be estimated or not. For further details see documentation for the parameter \code{nugget} from the \code{\link{variofit}} function.
}
  \item{nugget}{
    value for the nugget parameter. For further details see documentation for the parameter \code{nugget} from the \code{\link{variofit}} function.
}
  \item{fix.kappa}{
    logical, indicating whether the \code{kappa} parameter should be estimated or not. For further details see documentation for the parameter \code{fix.kappa} from the \code{\link{variofit}} function.
}
  \item{kappa}{
    value of the smoothness parameter. For further details see documentation for the parameter \code{kappa} from the \code{\link{variofit}} function.
}
  \item{max.dist.variogram}{
    a numerical value defining the maximum distance considered when fitting the variogram. For further details see documentation for the parameter \code{max.dist} from the \code{\link{variofit}} function.
}
}
\details{
  This function makes exhaustive use of the \code{\link{variofit}} function.
}
\value{
  A list with the following components:
  \item{trace.vari}{choosed theoretical variogram model}
  \item{trace.vari.array}{vector of all fitted theoretical variogram models}
}
\references{
  Giraldo, R. (2009) \emph{Geostatistical Analysis of Functional Data.} Ph.D. thesis. Universitat Politècnica de Catalunya.
}
\author{
  Ramon Giraldo \email{rgiraldoh@unal.edu.co}, \cr
  Pedro Delicado \email{pedro.delicado@upc.edu}, \cr
  Jorge Mateu \email{mateu@mat.uji.es}.
}
\keyword{spatial}
