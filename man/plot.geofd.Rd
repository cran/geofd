\name{plot.geofd}
\alias{plot.geofd}
\title{Plot Trace Variogram and adjusted models}
\description{
  This funcion produces a plot of an object of class geofd which contains...
}
\usage{
\method{plot}{geofd}(x, emp.trace.vari=x$emp.trace.vari,
                    trace.vari.array=x$trace.vari.array,
                    colors=bpy.colors(length(trace.vari.array)), \dots)
}
\arguments{
  \item{x}{a list containing elements \code{emp.trace.vari} and \code{trace.vari.array} described below. Typically an object of the class \code{"geofd"}. If not provided the arguments \code{emp.trace.vari} and \code{trace.vari.array} must be provided instead.}
  \item{emp.trace.vari}{empirical trace-variogram.}
  \item{trace.vari.array}{vector of variogram models.}
  \item{colors}{a character vector of color names used to plot each variogram model. Dimensions must be the same of \code{trace.vari.array}.}
  \item{\dots}{graphical arguments to be passed to \code{\link{plot}}.}
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
%\examples{
%data(plot.geofd)
%## maybe str(plot.geofd) ; plot(plot.geofd) ...
%}
\keyword{spatial}
\keyword{dplot}
