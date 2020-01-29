\encoding{utf8}
\name{geofd-package}
\alias{geofd-package}
\alias{geofd}
\docType{package}
\title{
Spatial Prediction for Function Value Data
}
\description{
Kriging based methods are used for predicting functional data (curves) with spatial dependence.
Initially the curves are pre-processed by fitting a Fourier or B-splines basis functions.
Then the spatial dependence among curves is estimated 
by means of the trace-variogram function. 
Finally the parameters for performing prediction by Ordinary Kriging 
at unsampled locations are estimated by solving a linear system based on the estimated trace-variogram. 
}

\details{
Index:
\preformatted{
geofd-package           Spatial prediction for function value data
}
}
\author{
Ramon Giraldo <rgiraldoh@unal.edu.co>, 
Pedro Delicado <pedro.delicado@upc.edu>, 
Jorge Mateu <mateu@mat.uji.es>

Maintainer: Pedro Delicado <pedro.delicado@upc.edu>
}
\keyword{spatial}
\keyword{datagen}

