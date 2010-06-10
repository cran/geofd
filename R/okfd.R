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
  # If cov.model is null then the best cov.model will be evaluated between spherical, exponential and gaussian
  if(!is.null(cov.model)){
    cov.model <- match.arg(cov.model, c("spherical","exponential","gaussian","matern"))
  }
  # nbasis, argvals and lambda are validated in runtime on their corresponding using function
  # max.dist, fix.nugget and nugget does not seem to be validated
#  if(is.null(new.coords)) stop("new.coords is not an optional parameter")
  if(ncol(new.coords)!=2) stop("new.coords must be an n x 2 matrix")

  # Argument type conversion
  new.coords <- as.matrix(new.coords)
  coords <- as.matrix(coords)

  if(smooth.type == "fourier"){
    n <- dim(data)[1]
  }

  s <- dim(data)[2] # number of sites

##################################################################
#  Creating a fd object
##################################################################

  rangeval <- range(argvals)

  if(smooth.type == "fourier"){

    period <- n
    basis <- create.fourier.basis(rangeval, nbasis, period)
    datafd <- Data2fd(argvals=argvals, y=data, basisobj=basis, lambda=lambda)

  }else{
    # 4 define un spline cubico
    norder <- 4
    basis <- create.bspline.basis(rangeval, nbasis, norder)

    datafdPar <- fdPar(basis, Lfdobj=2, lambda)
    smooth.datafd <- smooth.basis(argvals, data, datafdPar)
    datafd <- smooth.datafd$fd   
  }

##################################################################
# Calculating L2 norm among functions
# Integral del cuadrado de las diferencias entre las funciones
##################################################################

  if(smooth.type == "fourier"){
    M <- fourierpen(basis,Lfdobj=0)
  }else{
    M <- bsplinepen(basis,Lfdobj=0)
  }

  L2norm <- l2.norm(s, datafd, M)

##################################################################
# Empirical trace-variogram
##################################################################

  emp.trace.vari <- variog.okfd(coords, L2norm)
  Eu.d <- emp.trace.vari$Eu.d

##################################################################
# Fitting a theoretical variogram model to empirical trace-variogram
##################################################################

  if(smooth.type == "fourier"){
    # partial sill is set to the quantile 0.75
    sigma2.0 <- quantile(emp.trace.vari$v, 0.75)
  }else{
    # partial sill is set to the variance of the data
    sigma2.0 <- var(as.vector(data))
  }

  # range is set to the quantile 0.75
  phi.0 <- quantile(Eu.d, 0.75)

  if(is.null(max.dist.variogram)){
    max.dist.variogram <- max(emp.trace.vari$u)
  }

  # If cov.model is null then the best cov.model will be evaluated between spherical, exponential and gaussian  
  # If cov.model is NOT null then the argument received will be evaluated
  if(is.null(cov.model)){
    # The proposed variogram models are loaded
    models <- c("spherical","exponential","gaussian","matern")
  }else{
    # The proposed variogram model is loaded
    models <- cov.model
  }

  trace.vari.objects <- variofit.okfd(emp.trace.vari, models, sigma2.0, phi.0, fix.nugget, nugget, fix.kappa, kappa, max.dist.variogram)
  trace.vari <- trace.vari.objects$best
  trace.vari.array <- trace.vari.objects$fitted

##################################################################
# Distances among sites and distances to the prediction site
# Distance matrix among sampling sites and distance to NEW SITES
##################################################################

  new.s <- dim(new.coords)[1]
  new.Eu.d <- as.matrix(dist(rbind(coords,new.coords), method="euclidean"))
  new.Eu.d <- matrix(new.Eu.d[1:s,(s+1):(s+new.s)],nrow=s,ncol=new.s)

###################################################################
# Solving the system
##################################################################

  sigma2 <- trace.vari$cov.pars[1]
  leftmatrix <- sigma2 - cov.spatial(Eu.d,cov.model=trace.vari$cov.model, cov.pars=trace.vari$cov.pars, kappa=trace.vari$kappa)
  unosfila <- rep(1,s)
  leftmatrix <- rbind(leftmatrix,unosfila)
  unosycerocolumna <- c(rep(1,s),0)
  leftmatrix <- cbind(leftmatrix,unosycerocolumna)

  rightmatrix <- sigma2 - cov.spatial(new.Eu.d, cov.model=trace.vari$cov.model, cov.pars=trace.vari$cov.pars, kappa=trace.vari$kappa)
  unosfila <- rep(1,new.s)
  rightmatrix <- rbind(rightmatrix,unosfila)

  functional.kriging.weights <- solve(leftmatrix,rightmatrix)
  functional.kriging.weights.sinlagrange <- matrix(functional.kriging.weights[-(s+1),],nrow=s,ncol=new.s)
  sum(functional.kriging.weights.sinlagrange)

##################################################################
# Solution and prediction plot
##################################################################

  eval.data <- eval.fd(argvals,datafd)
  krig.new.data <- eval.data%*%functional.kriging.weights.sinlagrange

##################################################################
# Prediction variances
# Prediction variance
##################################################################

  vect.semiv <-rightmatrix[-(s+1),]
  varianza <-functional.kriging.weights.sinlagrange*vect.semiv 
  suma.varianza<-sum(varianza)
  pred.var <- suma.varianza + functional.kriging.weights[s+1,]

##################################################################
# Return:
##################################################################

  return.list <- list(
    coords=coords, 
    data=data, 
    argvals=argvals, 
    nbasis=nbasis, 
    lambda=lambda, 
    new.coords=new.coords,
    emp.trace.vari=emp.trace.vari,
    trace.vari=trace.vari,
#    Lo que viene en el parámetro $u son las distancias 
#    Eu.d=Eu.d, 
    new.Eu.d=new.Eu.d,  
    functional.kriging.weights=functional.kriging.weights, 
    krig.new.data=krig.new.data,
    pred.var=pred.var,
    trace.vari.array=trace.vari.array,
    datafd=datafd
  )
  #return.list$argnames <- argnames
  class(return.list) <- "geofd"

  return(return.list)

}
