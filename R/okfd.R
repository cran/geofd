# general parameters
# smoothing parameters
# variofit parameters
okfd <- function(new.coords, coords, data,
    smooth.type=NULL, nbasis=max(50,dim(data)[1]), argvals=seq(0, 1, len = dim(data)[1]), lambda=0,
    cov.model=NULL, fix.nugget=FALSE, nugget=0, fix.kappa=TRUE, kappa=0.5, max.dist.variogram=NULL)
# argnames=c("argument", "sites", "values"),
{

  # Loading required libraries
  require(fda)
  require(geoR)

  # Argument validation
  smooth.type <- match.arg(smooth.type, c("bsplines","fourier"))
  # The argument cov.model is validated if it is NOT null
  if(!is.null(cov.model)){
    cov.model <- match.arg(cov.model, c("spherical","exponential","gaussian","matern"))
  }
  if(is.null(new.coords)) stop("new.coords is not an optional parameter")
  if(ncol(new.coords)!=2) stop("new.coords must be an n x 2 matrix")
  # nbasis, argvals and lambda are validated in runtime on their corresponding using function
  # max.dist, fix.nugget and nugget does not seem to be validated

  # Argument type conversion
  new.coords <- as.matrix(new.coords)
  coords <- as.matrix(coords)

  # Number of sites
  s <- dim(data)[2]

  fdmodel <- .simple.fdmodel(new.coords, coords, data, smooth.type, nbasis, argvals, lambda, cov.model, fix.nugget, nugget, fix.kappa, kappa, max.dist.variogram)

  ##################################################################
  # Doing prediction
  ##################################################################

  prediction <- .okfd.predict(argvals, fdmodel$fdobjects$datafd, coords, new.coords, fdmodel$trace.vari.objects$best, fdmodel$emp.trace.vari$Eu.d)

  ##################################################################
  # Return
  ##################################################################

  return.list <- list(
    coords=coords,
    data=data,
    argvals=argvals,
    nbasis=nbasis,
    lambda=lambda,
    new.coords=new.coords,
    emp.trace.vari=fdmodel$emp.trace.vari,
    trace.vari=fdmodel$trace.vari.objects$best,
#    Lo que viene en el parámetro $u son las distancias 
#    Eu.d=Eu.d, 
    new.Eu.d=prediction$new.Eu.d,
    functional.kriging.weights=prediction$functional.kriging.weights,
    krig.new.data=prediction$pred,
    pred.var=prediction$var,
    trace.vari.array=fdmodel$trace.vari.objects$fitted,
    datafd=fdmodel$fdobjects$datafd
  )
  #return.list$argnames <- argnames
  class(return.list) <- "geofd"

  return(return.list)

}

".okfd.predict" <- function(argvals, datafd, coords, new.coords, trace.vari, Eu.d){

  s <- dim(coords)[1]

  # Distances among sites and distances to the prediction site
  # Distance matrix among sampling sites and distance to NEW SITES
  new.s <- dim(new.coords)[1]
  new.Eu.d <- as.matrix(dist(rbind(coords,new.coords), method="euclidean"))
  new.Eu.d <- matrix(new.Eu.d[1:s,(s+1):(s+new.s)], nrow=s, ncol=new.s)

  # Solving the system
  sigma2 <- trace.vari$cov.pars[1]
  leftmatrix <- sigma2 - cov.spatial(Eu.d, cov.model=trace.vari$cov.model, cov.pars=trace.vari$cov.pars, kappa=trace.vari$kappa)
  unosfila <- rep(1,s)
  leftmatrix <- rbind(leftmatrix, unosfila)
  unosycerocolumna <- c(rep(1,s),0)
  leftmatrix <- cbind(leftmatrix, unosycerocolumna)

  rightmatrix <- sigma2 - cov.spatial(new.Eu.d, cov.model=trace.vari$cov.model, cov.pars=trace.vari$cov.pars, kappa=trace.vari$kappa)
  unosfila <- rep(1, new.s)
  rightmatrix <- rbind(rightmatrix, unosfila)

  functional.kriging.weights <- solve(leftmatrix, rightmatrix)
  functional.kriging.weights.sinlagrange <- matrix(functional.kriging.weights[-(s+1),], nrow=s, ncol=new.s)
  sum(functional.kriging.weights.sinlagrange)

  # Solution
  eval.data <- eval.fd(argvals, datafd)
  krig.new.data <- eval.data%*%functional.kriging.weights.sinlagrange

  # Prediction variance
  vect.semiv <-rightmatrix[-(s+1),]
  varianza <-functional.kriging.weights.sinlagrange*vect.semiv 
  suma.varianza<-sum(varianza)
  pred.var <- suma.varianza + functional.kriging.weights[s+1,]

  return( list(pred=krig.new.data, var=pred.var, new.Eu.d=new.Eu.d, functional.kriging.weights=functional.kriging.weights) )

}

# This function is a simple use of the geofd modelling functions
# it should not be used from other functions than those of this package
".simple.fdmodel"<-function(new.coords, coords, data,
    smooth.type=NULL, nbasis=max(50,dim(data)[1]), argvals=seq(0, 1, len = dim(data)[1]), lambda=0,
    cov.model=NULL, fix.nugget=FALSE, nugget=0, fix.kappa=TRUE, kappa=0.5, max.dist.variogram=NULL){

  # Number of sites
  s <- dim(data)[2]

  # Smoothing original data
  fdobjects <- .create.fd.object(data, smooth.type, argvals, nbasis, lambda)

  # Calculating L2 norm among functions
  L2norm <- l2.norm(s, fdobjects$datafd, fdobjects$M)

  # Empirical trace-variogram
  emp.trace.vari <- trace.variog(coords, L2norm)

  # Fitting a theoretical variogram model to empirical trace-variogram
  if(smooth.type == "fourier"){
    # partial sill is set to the quantile 0.75
    sigma2.0 <- quantile(emp.trace.vari$v, 0.75)
  }else{
    # partial sill is set to the variance of the data
    sigma2.0 <- var(as.vector(data))
  }
  trace.vari.objects <- fit.tracevariog(emp.trace.vari=emp.trace.vari, models=cov.model, sigma2.0=sigma2.0, phi.0=quantile(emp.trace.vari$Eu.d, 0.75), fix.nugget, nugget, fix.kappa, kappa, max.dist.variogram)

  return( list(fdobjects=fdobjects, emp.trace.vari=emp.trace.vari, trace.vari.objects=trace.vari.objects) )

}

