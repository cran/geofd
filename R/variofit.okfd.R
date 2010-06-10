"variofit.okfd"<-function(emp.trace.vari, models, sigma2.0, phi.0, fix.nugget,
                        nugget, fix.kappa, kappa, max.dist.variogram){

  # Argument validation
  if(!identical(class(emp.trace.vari),"variogram")) stop("the parameter emp.trace.vari must be of class variogram")
  if(!identical(dim(colors),dim(models))) stop("dimensions of parameters colors and models must be identical")
  if(!is.numeric(sigma2.0)) stop("the partial sill paramter must be a number")
  if(!is.numeric(phi.0)) stop("the range parameter must be a number")
  if(!is.logical(fix.nugget)) stop("the parameter fix.nugget must be logical")
  if(!is.numeric(nugget)) stop("the nugget parameter must be a number")
  if(!is.logical(fix.kappa)) stop("the parameter fix.kappa must be logical")
  if(!is.numeric(kappa)) stop("the kappa parameter must be a number")
  if(!is.numeric(max.dist.variogram)) stop("the max.dist.variogram parameter must be a number")

  # A comparison variable for choosing the best variogram model is loaded
  trace.vari <- c(NULL)
  trace.vari["value"] <- Inf
  # Is created an array where the variofit results will be pushed
  trace.vari.array = c()
  # A loop is made in order to prove each model and choose the best using as criterion the variofit minimised sum of squares
  for(cont in c(1:length(models))){
    # The variofit function is called
    trace.vari.tmp <- variofit(emp.trace.vari, ini.cov.pars=c(sigma2.0,phi.0), max.dist=max.dist.variogram,
    fix.nugget=fix.nugget, nugget=nugget, fix.kappa=fix.kappa, kappa=kappa,
    cov.model=models[cont], messages=FALSE)
    # Each calculated trace variogram is pushed inside trace.vari.array
    trace.vari.array[[length(trace.vari.array)+1]] <- trace.vari.tmp
    # The last caculated variofit is compared against the last optimal variofit
    # If it is more optimal, then it is loaded in the variable wich will be returned
    if(as.numeric(trace.vari.tmp["value"])<as.numeric(trace.vari["value"])){
      trace.vari <- trace.vari.tmp
    }
  }

##################################################################
# Return:
##################################################################

  return(list(best=trace.vari, fitted=trace.vari.array))

}

"plot.geofd"<-function(x, emp.trace.vari=x$emp.trace.vari,
                      trace.vari.array=x$trace.vari.array,
                      colors=bpy.colors(length(trace.vari.array)), ... )
{

  if(missing(x)) x <- list(emp.trace.vari=emp.trace.vari, trace.vari.array=trace.vari.array)
  if(!is.list(trace.vari.array)) stop("the trace.vari.array parameter must be a list containing Variogram Model elements")
  if(!is.list(emp.trace.vari) || class(emp.trace.vari)!="variogram") stop("the emp.trace.vari parameter must be a list object of the class \"variogram\"")
  if(length(colors)!=length(trace.vari.array)) stop("the arguments \"colors\" and \"trace.vari.array\" must have the same length")

  # The empirical trace variogram is plotted
  plot(emp.trace.vari, xlab="Distance", ylab="Trace-Variogram", main="Trace-Variogram Cloud")

  # The argument legend parameter for the legend is loaded
  legend <- c("empirical trace variogram")

  cont <- 1
  # Each calculated trace variogram is plotted
  for(trace.vari in trace.vari.array){
    lines(trace.vari,col=colors[cont],lwd=2)
    legend <- c(legend, trace.vari$cov.model)
    cont <- cont+1
  }

  # Other arguments for the legend function are loaded
  colors <- c(1, colors)
  pch <- c(21, array(-1,length(trace.vari.array)) )
  lty <- c(0, array(1,length(trace.vari.array)) )
  pt.cex <- c(1, array(1,length(trace.vari.array)) )

  # The legend is plotted
  legend("topleft", "(x,y)", legend, col=colors, pch=pch, lty=lty, pt.cex=pt.cex)
}

