# Cross-Validation analysis for ordinary kriging for function value-data.

# coords, data, argvals=seq(0,1,len=dim(data)[1]),
# argnames=c("argument", "sites", "values"),
# array.nbasis=max(50,dim(data)[1]), array.lambda = 0,
# max.dist.variogram=NULL, nugget.fix=NULL

okfd.cv <- function(coords, data, argnames=c("argument", "sites", "values"),
    smooth.type=NULL, array.nbasis=max(50,dim(data)[1]), argvals=seq(0,1,len=dim(data)[1]), array.lambda=0,
    cov.model=NULL, fix.nugget=FALSE, nugget=0, fix.kappa=TRUE, kappa=0.5, max.dist.variogram=NULL)
{
  # Argument validation
  if(!is.vector(array.nbasis)) stop("the argument \"array.nbasis\" must be a vector")
  if(!is.vector(array.lambda)) stop("the argument \"array.lambda\" must be a vector")
  if(!is.vector(argnames) || length(argnames)!=3) stop("argnames must be a character vector of length 3")

  # Init local variables 
  m <- dim(data)[1]
  diff.argvals <- diff(argvals)
  s <- dim(coords)[1]
  mse.cv <- matrix(0,nrow=length(array.nbasis),ncol=length(array.lambda))
  krig.cv <- array(0,dim=c(length(array.nbasis),length(array.lambda),s,m),
  dimnames <- c("nbasis","lambda",argnames[2],argnames[1]))

  k <- 0
  k.opt <- 1
  l.opt <- 1

  # Loop over all number of basis functions parameters
  for (nbasis.k in array.nbasis){
    k <- k+1
    l <- 0
    # Loop over all smoothing penalization parameters
    for (lambda.l in array.lambda){
      l <- l+1
      # Loop over all sites number
      for (i in 1:s){
        # Prediction for the first 'i' sites using the last 's-i'
        res.okfd <- okfd( new.coords=as.data.frame(coords)[i,], coords=as.data.frame(coords)[-i,], data=data[,-i],
                          smooth.type=smooth.type, nbasis=nbasis.k, argvals=argvals, lambda=lambda.l,
                          cov.model=cov.model, fix.nugget=fix.nugget, nugget=nugget, fix.kappa=fix.kappa, kappa=kappa, max.dist.variogram=max.dist.variogram)
        # Predicted data is saved 
        krig.cv[k,l,i,] <- res.okfd$krig.new.data
        # An error measure is calculated
        aux <- (data[,i]-res.okfd$krig.new.data)^2
        aux <- diff.argvals * (aux[1:(m-1)]+aux[2:m])/2
        mse.cv[k,l] <- mse.cv[k,l] + sum( aux ) 
      }
      if (mse.cv[k,l] <= mse.cv[k.opt,l.opt]){
        k.opt <- k
        l.opt <- l
      }
    }
  }

  mse.cv.opt <- mse.cv[k.opt,l.opt]

  return(list(
    k.opt= k.opt,
    l.opt= l.opt,
    krig.cv= krig.cv,
    mse.cv= mse.cv,
    mse.cv.opt= mse.cv.opt
  ))

}

