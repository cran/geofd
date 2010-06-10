\encoding{utf8}
\name{okfd.cv}
\alias{okfd.cv}
\alias{geofd.okfd.cv}
\title{
  Function for doing Cross-Validation analysis for Ordinary Kriging for function-value data
}
\description{
  Unreviewed
}
\usage{
  okfd.cv(coords, data, argnames=c("argument", "sites", "values"),
          smooth.type=NULL, array.nbasis=max(50,dim(data)[1]), argvals=seq(0,1,len=dim(data)[1]), array.lambda=0,
          cov.model=NULL, fix.nugget=FALSE, nugget=0, fix.kappa=TRUE, kappa=0.5, max.dist.variogram=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coords}{
    coordinates of the sites where functional data are observed (dim: s by 2)
  }
  \item{data}{
    matrix with values for the observed functions (dim: m by s)
  }
  \item{argnames}{
    a character vector of length three containing: the name of the argument (argvals), a description of the sites (coord), the name of the observed function values.
  }
  \item{smooth.type}{
    a string with the name of smoothing method to be applied to \code{data}. Available choices are: "bsplines" and "fourier". 
  }
  \item{array.nbasis}{
    array with values for the number of elements in the cubic B-spline basis.
  }
  \item{argvals}{
    a set of argument values. (length: m)
  }
  \item{array.lambda}{
    array of penalization parameters for smoothing the observed functions.
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
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A list with the following components:
  \item{k.opt}{unreviewed}
  \item{l.opt}{unreviewed}
  \item{krig.cv}{unreviewed}
  \item{mse.cv}{unreviewed}
  \item{mse.cv.opt}{unreviewed}
}
\references{
  Giraldo, R. (2009) \emph{Geostatistical Analysis of Functional Data.} Ph.D. thesis. Universitat Politècnica de Catalunya.
}
\author{
  Ramon Giraldo \email{rgiraldoh@unal.edu.co}, \cr
  Pedro Delicado \email{pedro.delicado@upc.edu}, \cr
  Jorge Mateu \email{mateu@mat.uji.es}.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.

%## The function is currently defined as

%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}
\keyword{datagen}
\keyword{smooth}
