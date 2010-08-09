\encoding{utf8}
\name{fit.tracevariog}
\alias{fit.tracevariog}
\alias{geofd.fit.tracevariog}
\title{
  Variogram Estimation
}
\description{
  Fits a parametric model to a empirical variogram and estimates covariance parameters. Aditionally all fitted variogram models are plotted for verification purpose.
}
\usage{
fit.tracevariog(emp.trace.vari, models, sigma2.0, phi.0,
                fix.nugget=FALSE, nugget=0,
                fix.kappa=TRUE, kappa=0.5,
                max.dist.variogram=NULL)
}
\arguments{
  \item{emp.trace.vari}{
    empirical trace-variogram. An object returned from the \code{\link{variog}} function.
}
  \item{models}{
    a character vector of correlation function names used in geoR against which empirical trace variogram will be fitted.
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
  This function makes use of the \code{\link{variofit}} function.

  \bold{Variogram models an parameters}

  When the \code{cov.model} parameter is \code{NULL} a function determines the optimal model between spherical, exponential gaussian and matern using the received parameters. The arguments \code{sigma2.0} and \code{phi.0} are used as initial values for fitting each variogram model.

  The parameters \code{fix.nugget}, \code{nugget}, \code{fix.kappa}, \code{kappa} and \code{max.dist.variogram} are the same for each variogram model specified in \code{models}.

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
