#the binning is
#     defined as follows:
#
#       1. read the argument ‘max.dist’. If not provided it is set to
#          the maximum distance between the pairs of points.
#
#       2. the center of the bins are initially defined by the sequence
#          ‘u = seq(0, max.dist, l = 13)’.
#
#       3. the interval spanned by each bin is given by the mid-points
#          between the centers of the bins.
#
#     If an vector is passed to the argument ‘breaks’ its elements are
#     taken as the limits of the bins (classes of distance) and the
#     argument ‘uvec’ is ignored.

# nugget.tolerance, uvec y breaks son solo para bin variogram
"trace.variog"<-function(coords, L2norm, bin=FALSE, max.dist, uvec="default", breaks="default", nugget.tolerance){

  # Argument validation
  if(is.null(coords)) stop("coords is not an optional parameter")
  if(ncol(coords)!=2) stop("coords must be an n x 2 matrix")
  if(!isSymmetric(L2norm)) stop("L2norm must be a symmetric matrix")
  if(sum(diag(L2norm))!=0) stop("each element of the diagonal of L2norm must zero")

  # Euclidian distance among sites
  Eu.d <- as.matrix(dist(coords, method="euclidean"))

  # Los valores del variograma son las diferencias punto a punto entre todas
  # las curvas de L2norm. La extracción de los valores del variograma
  # se hace quitando la parte superior de la matriz L2norm excluyendo la
  # diagonal y pasando secuencialmente por filas el resultado a un arreglo
  vtemp <- array(as.dist(L2norm))

  # Los valores de distancia se extraen igual que los valores del variograma
  utemp <- array(as.dist(Eu.d))

  # Se define temporalmente el valor de máxima distancia para utilizarse
  # en .define.bins y como valor definitivo de max.dist en cloud
  if(missing(max.dist)){
    umax <- max(utemp)
  }else{
    umax <- max(utemp[utemp < max.dist])
  }

  # |bin|<---lag--->|bin|<---lag--->|bin|<---lag--->|bin|<---lag--->|bin|

  if(bin){
    # Si no hay valor de nugget.tolerance o es mas pequeño que cierto valor
    # se define el nugget.tolerance y se indica que se debe remover el valor del primer bin
    if(missing(nugget.tolerance) || nugget.tolerance < 1e-11){
      nugget.tolerance <- 1e-12
      nt.ind <- FALSE
    }else{
      if(mode(nugget.tolerance) != "numeric") stop("nugget.tolerance must be numeric")
      nt.ind <- TRUE
    }
    # If the minimum distance is minor than the nugget.tolerance then it is
    # indicated that the first bin should NOT be removed
    min.dist <- min(utemp)
    if(min.dist < nugget.tolerance){
      nt.ind <- TRUE
    }
    # The bins are defined
    dbins <- .define.bins(max.dist=umax, nugget.tolerance=nugget.tolerance, uvec=uvec, breaks=breaks)
    # The max.dist value is defined
    if(missing(max.dist)){
      max.dist <- max(dbins$bins.lim)
    }
    # The first bin is adjusted
    if(dbins$bins.lim[1] < 1e-16){
      dbins$bins.lim[1] <- 0
    }
    if(!nt.ind){
      dbins$uvec <- dbins$uvec[-1]
      dbins$bins.lim <- dbins$bins.lim[-1]
    }

    u <- rep(NA, length(dbins$bins.lim)-1)
    v <- rep(NA, length(dbins$bins.lim)-1)
    for(i in seq(1,length(dbins$bins.lim)-1)){
      # vector defining points belong which lag
      # ToDo: specify what happens if a point is exactly on the edge of the lag
      lagpoints <- dbins$bins.lim[i]<utemp & utemp<dbins$bins.lim[i+1]
      # If there is any point inside the lag then it is used in the variogram
      if(any(lagpoints)){
        u[i] <- mean(utemp[lagpoints])
        v[i] <- mean(vtemp[lagpoints])
      }
    }
    output.type <- "bin"
    emp.trace.vari <- list(bins.lim=dbins$bins.lim, nugget.tolerance=nugget.tolerance)
#    emp.trace.vari <- list(bins.lim=dbins$bins.lim, nugget.tolerance=nugget.tolerance, uvec=dbins$uvec)
  }else{
    u <- utemp
    v <- vtemp
    max.dist <- umax
    output.type <- "cloud"
    emp.trace.vari <- list()
  }

#  plot(u, v)

  # ToDo: que hacer si los bins son muy pequeños, manejar con nugget.tolerance

  # Se crea el objeto trace-variogram a retornar
  # ToDo: Se está retornando los objetos 'Eu.d' y 'L2norm' que son formas aumentadas de los objetos 'u' y 'v' cuando el tipo de variograma es cloud
  # El objeto Eu.d luego es utilizado para calcular el rango del semivariograma
  # Eu.d y L2norm pueden ser obtenidos de u y v haciendo un ciclo que itere sobre todos los valores:
  # [a, b, c]
  # y cree una matrix de la forma:
  # 0 a b
  # a 0 c
  # b c 0
  emp.trace.vari <- c(emp.trace.vari, list(u=u, v=v, output.type=output.type, max.dist=max.dist, Eu.d=Eu.d, L2norm=L2norm))

  # Se asigna la clase variogram al objeto de retorno para que luego
  # pueda ser utilizado como cualquier variomodel de geoR
  class(emp.trace.vari) <- "variogram"

  return(emp.trace.vari)
}

".create.fd.object"<-function(data, smooth.type, argvals, nbasis, lambda){

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

  return( list(M=M, datafd=datafd) )

}

