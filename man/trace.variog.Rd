\encoding{utf8}
\name{trace.variog}
\alias{trace.variog}
\alias{geofd.trace.variog}
\title{
  Empirical Variograms for function-value data
}
\description{
  Computes empirical trace-variograms using the L2 norm matrix for the semivariance values. Output a variogram cloud.
}
\usage{
  trace.variog(coords, L2norm, bin=FALSE, max.dist, uvec="default",
              breaks="default", nugget.tolerance)
}
\arguments{
  \item{coords}{an s x 2 matrix containing the coordinates of the n sites where functional data are observed.}
  \item{L2norm}{L2 norm among functions.}
  \item{bin}{logical, indicating whether the output is the binned variogram.}
  \item{max.dist}{a numerical value defining the maximum distance for the variogram.}
  \item{uvec}{a vector with values defining the centers of the bins or the number of bins as in \code{\link{variog}}. Only used when 'bin = TRUE'.}
  \item{breaks}{a vector with values defining the variogram binning as in \code{\link{variog}}. Only used when 'bin = TRUE'.}
  \item{nugget.tolerance}{a numeric value defining the shortest lag distance. Only used when 'bin = TRUE'.}
}
\details{
  \bold{Binned variogram}
  This is just a visual feature adapted from the cloud variogram and it doesn't have any relation against the fitting of the variogram model or the calculation of the predictions.

  The binning is made in the same manner as in \code{\link{variog}} because it uses the \code{.define.bins} function.
}
\value{
  An object of the class variogram as the the \code{\link{variog}} function which is a list with the following components:
  \item{u}{a vector with distances.}
  \item{v}{a vector with estimated variogram values at distances given in u.}
  \item{max.dist}{maximum distance of the variogram.}
  \item{output.type}{variogram type.}
  \item{Eu.d}{euclidian distance array among sites.}
  \item{L2norm}{echoes the 'L2norm' argument.}
  \item{bins.lim}{limits defining the interval spanned by each bin as in \code{\link{variog}}. Only returned when 'bin = TRUE'.}
  \item{nugget.tolerance}{echoes the 'nugget.tolerance' argument.}
%  \item{uvec}{variogram bin centers.}
}
\references{
Giraldo, R. (2009) \emph{Geostatistical Analysis of Functional Data.} Ph.D. thesis. Universitat Politecnica de Catalunya.

Giraldo, R., Delicado, P. and Mateu, J. (2012) \emph{\pkg{geofd}: An R package for function-valued geostatistical prediction}. 
Revista Colombiana de Estadistica. 35, 385-407.
}
\author{
  Ramon Giraldo \email{rgiraldoh@unal.edu.co}, \cr
  Pedro Delicado \email{pedro.delicado@upc.edu}, \cr
  Jorge Mateu \email{mateu@mat.uji.es}.
}
\examples{

# First example: creating a binned variogram
# okfd first example

data(maritimes.avg)
data(maritimes.coords)
data(maritimes.data)
coord.cero <- matrix(c(-64.06, 45.79),nrow=1,ncol=2)
n<-dim(maritimes.data)[1]
argvals<-seq(1,n, by=1)
okfd.res<-okfd(new.coords=coord.cero, coords=maritimes.coords,
              data=maritimes.data, nbasis=65, argvals=argvals, fix.nugget=TRUE,
              kappa=0.7)

# Calculating the empiricial trace bin variogram
new.emp.trace.vari <- trace.variog(coords=okfd.res$coords,
                      L2norm=okfd.res$emp.trace.vari$L2norm, bin=TRUE)

# The empirical trace cloud variogram is replaced with the trace bin variogram
okfd.res$emp.trace.vari <- new.emp.trace.vari

# The modified okfd result is plotted
plot(okfd.res)

}
\keyword{spatial}
