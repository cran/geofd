\encoding{utf8}
\name{okfd}
\alias{okfd}
\alias{geofd.okfd}
\title{
  Function for doing Ordinary Kriging for function-value Data
}
\description{
  This function allows to carry out prediction by ordinary kriging for function-value data by considering a Fourier or B-splines basis for smoothing the observed data set
}
\usage{
okfd(new.coords, coords, data,
    smooth.type=NULL, nbasis=max(50,dim(data)[1]), argvals=seq(0, 1, len = dim(data)[1]), lambda=0,
    cov.model=NULL, fix.nugget=FALSE, nugget=0, fix.kappa=TRUE, kappa=0.5, max.dist.variogram=NULL)
% argnames = c("argument", "sites", "values"),
}
\arguments{
  \item{new.coords}{
    an n x 2 matrix containing the coordinates of the new n sites where functional kriging has to be done % dim: new.s by 2
  }
  \item{coords}{
    an s x 2 matrix containing the coordinates of the n sites where functional data are observed % dim: s by 2
  }
  \item{data}{
    an m x s matrix with values for the observed functions % dim: m by s
  }
  \item{smooth.type}{
    a string with the name of smoothing method to be applied to \code{data}. Available choices are: "bsplines" and "fourier". 
  }
  \item{nbasis}{
    a numeric value defining the number of basis functions used to smooth the discrete data set recorded at each site
  }
  \item{argvals}{
    a vector of argument values corresponding to the observations in matrix data % length: m
  }
  \item{lambda}{
    optional. Penalization parameter for smoothing the observed functions.
  }
  \item{cov.model}{
    a string with the name of the correlation function. Default is NULL, see \code{DETAILS} below. For further details see documentation for the parameter \code{cov.model} from the \code{\link{variofit}} function.
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
%  \item{argnames}{
%    optional. A character vector of length 3 containing: the name of the parameter (argvals), a description of the sites (coord) and the name of the observed function values
%  }
}
\details{
  \bold{Fitting covariance model}

  When the \code{cov.model} parameter is \code{NULL} a function determines the optimal model between spherical, exponential gaussian and matern using the received parameters. The arguments \code{fix.nugget}, \code{nugget}, \code{fix.kappa}, \code{kappa} and \code{max.dist.variogram} are used as initial values for fitting each variogram model.

  \bold{Penalization parameter}

  The penalization parameter \code{lambda} is used in both smoothing methods. When the selected smooth method is:
  \enumerate{
  \item
  bsplines, the function which uses it is \code{\link{fdPar}}
  \item
  fourier, the function which uses it is \code{\link{Data2fd}}
  }
  
  \bold{Functional data object}

  The function which creates the functional data object is determined based on the selected smooth method:
  \itemize{
  \item When it is bsplines, the functional data object must be created using two different functions, \code{\link{fdPar}} and \code{\link{smooth.basis}} in order to include the penalization parameter lambda.
  \item When it is fourier, the functional data object is directly returned by \code{\link{Data2fd}} because it includes the penalization parameter, the basis object, the argument values and the data, all at the same time.
  }

  \bold{Roughness penalty matrix and Basis Object}

  When the roughness penalty matrix is being calculated using the \code{\link{l2.norm}} function, the following considerations are taked in count:
  \itemize{
  \item The differential operator passed as parameter for both \code{fourierpen} and \code{bsplinepen} is always zero.
  \item When the selected smooth method is bsplines, the basis object passed to \code{\link{bsplinepen}} is the output of the function \code{\link{create.bspline.basis}} using 
\code{argvals} as the \code{rangeval} parameter, \code{nbasis} as the number of basis functions parameter and the default order of b-splines, which is four, a cubic spline, as the \code{norder} parameter.
  \item When the selected smooth method is fourier, the basis object is the output of the function \code{fourierpen}. The parmeters \code{rangeval} and \code{nbasis} are the same as for \code{\link{create.bspline.basis}} and the number of observations on each curve as the \code{period} parameter.
  }

}
\value{
  A list with the following components:
  \item{coords}{a matrix containing the coordinates of the sites where functional data are observed.}
  \item{data}{a matrix with values for the observed functions.}
  \item{argvals}{a vector of argument values corresponding to the observations in matrix data}
  \item{nbasis}{a numeric value defining the number of basis functions used to smooth the discrete data set recorded at each site.}
  \item{lambda}{penalization parameter for smoothing the observed functions.}
  \item{new.coords}{matrix containing the coordinates of the new sites where functional kriging has to be done.}
  \item{emp.trace.vari}{empirical trace-variogram.}
  \item{trace.vari}{choosed theoretical variogram model}
%  \item{Eu.d}{euclidian distance among observed sites.}
  \item{new.Eu.d}{distance matrix among sampling and new sites.}
  \item{functional.kriging.weights}{functional kriging weights.}
  \item{krig.new.data}{predicted values for the new sites.}
  \item{pred.var}{prediction variance.}
  \item{trace.vari.array}{vector of all fitted variogram models.}
  \item{datafd}{a functional data object containing a smooth of the data.}

%  \item{smooth.type}{
%    a string with the name of smoothing method to be applied to \code{data}. Available choices are: "bsplines" and "fourier". 
%  }
%  \item{cov.model}{
%    a string with the name of the correlation function. Default is NULL, see \code{DETAILS} below. For further details see documentation for the parameter \code{cov.model} from the \code{\link{variofit}} function.
%  }
%  \item{fix.nugget}{
%    logical, indicating whether the \code{nugget} parameter should be estimated or not. For further details see documentation for the parameter \code{nugget} from the \code{\link{variofit}} function.
%  }
%  \item{nugget}{
%    value for the nugget parameter.
%  }
%  \item{fix.kappa}{
%    logical, indicating whether the \code{kappa} parameter should be estimated or not.
%  }
%  \item{kappa}{
%    value of the smoothness parameter. For further details see documentation for the parameter \code{kappa} from the \code{\link{variofit}} function.
%  }
%  \item{max.dist.variogram}{
%    a numerical value defining the maximum distance considered when fitting the variogram. For further details see documentation for the parameter \code{max.dist} from the \code{\link{variofit}} function.
%  }
}
\references{
  Giraldo, R. (2009) \emph{Geostatistical Analysis of Functional Data.} Ph.D. thesis. Universitat Politècnica de Catalunya.
}
\author{
  Ramon Giraldo \email{rgiraldoh@unal.edu.co}, \cr
  Pedro Delicado \email{pedro.delicado@upc.edu}, \cr
  Jorge Mateu \email{mateu@mat.uji.es}.
}
\seealso{
  \code{\link{l2.norm}} for calculating L2 norm among functions, \code{\link{variog.okfd}} for functional empirical trace variogram calculation, \code{\link{variofit.okfd}} for fitting a variogram model in the funcional scenario.
}
\examples{

# First example: prediction in only one site and smoothing using B-splines
# Data
data(maritimes)
coord.cero <- matrix(c(-64.06, 45.79),nrow=1,ncol=2)
n<-dim(maritimes.data)[1]
argvals<-seq(1,n, by=1)
# Prediction by okfd
okfd.res<-okfd(new.coords=coord.cero, coords=maritimes.coords, data=maritimes.data,
              nbasis=65, argvals=argvals, fix.nugget=TRUE, kappa=0.7)
# The smoothed and predicted curves, and predicted site average values are plotted
plot(okfd.res$datafd, lty=1,col=8, xlab="Day", ylab="Temperature (Degrees C)", main="Smoothed")
lines(okfd.res$argvals, okfd.res$krig.new.data, col=1, lwd=2, type="l", lty=1, main="Predictions", xlab="Day", ylab="Temperature (Degrees C)")
lines(maritimes.avg,  type="p", pch=20,cex=0.5, col=2, lwd=1)

# Second example: prediction in multiple sites and smoothing using Fourier basis functions
# Preparing example parameters
data(maritimes)
n<-dim(maritimes.data)[1]
argvals<-seq(1,n, by=1)
col1<-sample( (min(maritimes.coords[,1])*100):(max(maritimes.coords[,1])*100), 10, replace=TRUE)/100
col2<-sample( (min(maritimes.coords[,2])*100):(max(maritimes.coords[,2])*100), 10, replace=TRUE)/100
new.coords <- cbind(col1,col2)
# Prediction by okfd
okfd.res<-okfd(new.coords=new.coords, coords=maritimes.coords, data=maritimes.data,
              smooth.type="fourier", nbasis=65, argvals=argvals)
# The smoothed and predicted curves are plotted
par(mfrow=c(1,2))
plot(okfd.res$datafd, lty=1,col=8, xlab="Day", ylab="Temperature (Degrees C)", main="Smoothed")
matplot(okfd.res$argvals, okfd.res$krig.new.data, col=1, lwd=1, type="l", lty=1, main="Predictions", xlab="Day", ylab="Temperature (Degrees C)")

}
\keyword{spatial}
\keyword{datagen}
\keyword{smooth}
