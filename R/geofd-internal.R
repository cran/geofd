.create.fd.object <-
function(data, smooth.type, argvals, nbasis, lambda){

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

.okfd.predict <-
function(argvals, datafd, coords, new.coords, trace.vari, Eu.d){

  s <- dim(coords)[1]

  # Distances among sites and distances to the prediction site
  # Distance matrix among sampling sites and distance to NEW SITES
  new.s <- dim(new.coords)[1]
  new.Eu.d <- as.matrix(dist(rbind(coords,new.coords), method="euclidean"))
  new.Eu.d <- matrix(new.Eu.d[1:s,(s+1):(s+new.s)], nrow=s, ncol=new.s)

  # Solving the system
  sigma2 <- trace.vari$cov.pars[1]
  leftmatrix <- sigma2 - .cov.spatial(Eu.d, cov.model=trace.vari$cov.model, cov.pars=trace.vari$cov.pars, kappa=trace.vari$kappa)
  unosfila <- rep(1,s)
  leftmatrix <- rbind(leftmatrix, unosfila)
  unosycerocolumna <- c(rep(1,s),0)
  leftmatrix <- cbind(leftmatrix, unosycerocolumna)

  rightmatrix <- sigma2 - .cov.spatial(new.Eu.d, cov.model=trace.vari$cov.model, cov.pars=trace.vari$cov.pars, kappa=trace.vari$kappa)
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

.simple.fdmodel <-
function(new.coords, coords, data,
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

# Functions copied from geoR in January 2020, because geoR was close to be uncontinued
.cov.spatial <-  # copied from geoR in January 2020, because geoR was close to be uncontinued
  function(obj, cov.model = "matern",
           cov.pars = stop("no cov.pars argument provided"),
           kappa = 0.5)
  {
    fn.env <- sys.frame(sys.nframe())
    .check.cov.model(cov.model=cov.model, cov.pars=cov.pars, kappa=kappa,
                     env=fn.env, output=FALSE)
    phi <- get("phi", envir=fn.env)
    sigmasq <- get("sigmasq", envir=fn.env)
    ##
    ## computing correlations/covariances
    ##
    #  covs <- array(0, dim = dim(obj))
    covs <- obj; covs[] <- 0 
    for(i in 1:get("ns", envir=fn.env)) {
      if(phi[i] < 1e-16)
        cov.model[i] <- "pure.nugget"
      obj.sc <- obj/phi[i]
      cov.values <- switch(cov.model[i],
                           pure.nugget = rep(0, length(obj)),
                           wave = (1/obj) * (phi[i] * sin(obj.sc)),
                           exponential = exp( - (obj.sc)),
                           matern = {
                             if(kappa[i] == 0.5) exp( - (obj.sc))
                             else matern(u = obj, phi = phi[i], kappa = kappa[i])},
                           gaussian = exp( - ((obj.sc)^2)),
                           spherical = ifelse(obj < phi[i], (1 - 1.5 * (obj.sc) +
                                                               0.5 * (obj.sc)^3), 0),
                           circular = {
                             obj.sc[obj.sc > 1] <- 1;
                             ifelse(obj < phi[i], (1 - (2 * ((obj.sc) *
                                                               sqrt(1 - ((obj.sc)^2)) +
                                                               asin(obj.sc)))/pi), 0)
                           },
                           cubic = {
                             ifelse(obj < phi[i], (1 - (7 * (obj.sc^2) -
                                                          8.75 * (obj.sc^3) +
                                                          3.5 * (obj.sc^5) -
                                                          0.75 * (obj.sc^7))), 0)
                           },
                           power = (obj)^phi,
                           powered.exponential = exp( - ((obj.sc)^kappa[i])),
                           cauchy = (1 + (obj.sc)^2)^(-kappa[i]),
                           gneiting = {
                             obj.sc <- 0.301187465825 * obj.sc;   
                             t2 <- (1 - obj.sc);
                             t2 <- ifelse(t2 > 0, (t2^8), 0);
                             (1 + 8 * obj.sc + 25 * (obj.sc^2) + 32 * (obj.sc^3)) * t2
                           },
                           gencauchy = (1 + (obj.sc)^kappa[2])^(-kappa[1]/kappa[2]),
                           gneiting.matern = {
                             obj.sc <- 0.301187465825 * obj.sc/kappa[2] ;
                             t2 <- (1 - obj.sc);
                             t2 <- ifelse(t2 > 0, (t2^8), 0);
                             cov.values <- (1 + 8 * obj.sc + 25 * (obj.sc^2) + 32 * (obj.sc^3)) * t2;
                             cov.values * matern(u = obj, phi = phi[i], kappa = kappa[1])
                             
                           },
                           stop("wrong or no specification of cov.model")
      )
      if(cov.model[i] == "power"){
        A <- (max(cov.values)/sqrt(pi))*gamma((1+phi[i])/2)*gamma(1-(phi[i]/2))
        ## 1:2 below ensures valid results for 1 and 2D
        A <- A * max(gamma(phi[i]+(1+(1:2))/2)/(gamma(1+phi[i])*gamma((1+(1:2))/2)))
        cov.values <- A - cov.values
        cov.values <- cov.values/max(cov.values)
      }
      cov.values <- ifelse(obj < 1e-16, sigmasq[i], sigmasq[i] * cov.values)
      covs <- covs + cov.values
    }
    #  if(all(cov.model == "power"))
    #    covs <- max(covs) - covs
    #  else covs[obj < 1e-16] <- sum(sigmasq)
    if(sum(sigmasq) < 1e-16) covs[obj < 1e-16] <- 1
    if(any(!is.finite(covs))) warning("Infinity value in .cov.spatial")
    if(any(is.na(covs))) warning("NA value in .cov.spatial")
    if(any(is.nan(covs))) warning("NaN value in .cov.spatial")
    return(covs)
}


.variofit <- # copied from geoR in January 2020, because geoR was close to be uncontinued
  function (vario, ini.cov.pars, cov.model,
            fix.nugget = FALSE, nugget = 0, 
            fix.kappa = TRUE, kappa = 0.5,
            simul.number = NULL,  max.dist = vario$max.dist,
            weights, minimisation.function,
            limits = pars.limits(), messages, ...) 
  {
    call.fc <- match.call()
    if(missing(messages))
      messages.screen <- as.logical(ifelse(is.null(getOption("geoR.messages")), TRUE, getOption("geoR.messages")))
    else messages.screen <- messages
    if(length(class(vario)) == 0 || all(class(vario) != "variogram"))
      warning("object vario should preferably be of the geoR's class \"variogram\"")
    if(!missing(ini.cov.pars)){
      if(any(class(ini.cov.pars) == "eyefit"))
        cov.model <- ini.cov.pars[[1]]$cov.model
      if(any(class(ini.cov.pars) == "variomodel"))
        cov.model <- ini.cov.pars$cov.model
    }
    if(missing(cov.model)) cov.model <- "matern"
    cov.model <- match.arg(cov.model, choices =  c("matern", "exponential", "gaussian", "spherical",
      "circular", "cubic", "wave", "linear", "power",
      "powered.exponential", "stable", "cauchy", "gencauchy",
      "gneiting", "gneiting.matern", "pure.nugget"))
    if(cov.model == "stable") cov.model <- "powered.exponential"
    if(cov.model == "powered.exponential")
      if(limits$kappa["upper"] > 2) limits$kappa["upper"] <- 2
    ##  if(cov.model == "matern" | cov.model == "    powered.exponential" | 
    ##     cov.model == "cauchy" | cov.model == "gneiting.matern")
    ##    fix.kappa <- TRUE
    if(missing(weights)){
      if(vario$output.type == "cloud") weights <- "equal"
      else weights <- "npairs"
    }
    else
      weights <- match.arg(weights, choices = c("npairs", "equal", "cressie"))
    if(messages.screen){
      cat(paste(".variofit: covariance model used is", cov.model, "\n"))
      cat(paste(".variofit: weights used:", weights, "\n"))
    }
    #  if(missing(minimisation.function)){
    #    if(weights == "equal") minimisation.function <- "nls"
    #    else minimisation.function <- "optim"
    #  }
    if(missing(minimisation.function))
      minimisation.function <- "optim"
    if(any(cov.model == c("linear", "power")) & minimisation.function == "nls"){
      cat("warning: minimisation function nls can not be used with given cov.model.\n          changing for \"optim\".\n")
      minimisation.function <- "optim"
    }
    if(minimisation.function == "nls" & weights != "equal"){
      warning(".variofit: minimisation function nls can only be used with weights=\"equal\".\n          changing for \"optim\".\n")
      minimisation.function <- "optim"
    }
    if (is.matrix(vario$v) & is.null(simul.number)) 
      stop("object in vario$v is a matrix. This function works for only 1 empirical variogram at once\n")
    if (!is.null(simul.number)) 
      vario$v <- vario$v[, simul.number]
    ##
    ## Setting maximum distance
    ##
    if(mode(max.dist) != "numeric" || length(max.dist) > 1)
      stop("a single numerical value must be provided in the argument max.dist") 
    if (max.dist == vario$max.dist) 
      XY <- list(u = vario$u, v = vario$v, n=vario$n)
    else
      XY <- list(u = vario$u[vario$u <= max.dist],
                 v = vario$v[vario$u <= max.dist],
                 n = vario$n[vario$u <= max.dist])
    if(cov.model == "pure.nugget"){
      ##
      ## parameter estimation for model which does not require numerical minimisation
      ##
      minimisation.function <- "not used"
      message <- "correlation function does not require numerical minimisation"
      if(weights == "equal") lm.wei <- rep(1, length(XY$u))
      else lm.wei <- XY$n
      if(cov.model == "pure.nugget"){
        if(fix.nugget){
          temp <- lm((XY$v-nugget) ~ 1, weights = lm.wei)
          cov.pars <- c(temp$coef, 0)
        }
        else{
          temp <- lm(XY$v ~ 1, weights = lm.wei)
          nugget <- temp$coef
          cov.pars <- c(0,0)
        }
      }
      value <- sum((temp$residuals)^2)
    }
    else{
      if(messages.screen)
        cat(paste(".variofit: minimisation function used:", minimisation.function, "\n"))
      ##
      ## setting things for numerical minimisation
      ##
      ##  Checking initial values
      ##
      umax <- max(vario$u)
      vmax <- max(vario$v)
      if(missing(ini.cov.pars)){
        ini.cov.pars <- as.matrix(expand.grid(c(vmax/2, 3*vmax/4, vmax),
                                              seq(0, 0.8*umax, len=6)))
        if(!fix.nugget)
          nugget <- unique(c(nugget, vmax/10, vmax/4, vmax/2))
        if(!fix.kappa)
          kappa <- unique(c(kappa, 0.25, 0.5, 1, 1.5, 2))
        if(messages.screen)
          warning("initial values not provided - running the default search")
      }
      else{
        if(any(class(ini.cov.pars) == "eyefit")){
          init <- nugget <- kappa <- NULL
          for(i in 1:length(ini.cov.pars)){
            init <- drop(rbind(init, ini.cov.pars[[i]]$cov.pars))
            nugget <- c(nugget, ini.cov.pars[[i]]$nugget)
            if(cov.model == "gneiting.matern")
              kappa <- drop(rbind(kappa, ini.cov.pars[[i]]$kappa))
            else
              kappa <- c(kappa, ini.cov.pars[[i]]$kappa)
          }
          ini.cov.pars <- init
        }
        if(any(class(ini.cov.pars) == "variomodel")){
          nugget <- ini.cov.pars$nugget
          kappa <- ini.cov.pars$kappa
          ini.cov.pars <- ini.cov.pars$cov.pars
        }
      }
      if(is.matrix(ini.cov.pars) | is.data.frame(ini.cov.pars)){
        ini.cov.pars <- as.matrix(ini.cov.pars)
        if(nrow(ini.cov.pars) == 1)
          ini.cov.pars <- as.vector(ini.cov.pars)
        else{
          if(ncol(ini.cov.pars) != 2)
            stop("\nini.cov.pars must be a matrix or data.frame with 2 components: \ninitial values for sigmasq (partial sill) and phi (range parameter)\n")
        }
      }
      else
        if(length(ini.cov.pars) > 2)
          stop("\nini.cov.pars must provide initial values for sigmasq and phi\n")
      ##
      ## Preparing grid of initial values and choosing the best
      ##
      if(is.matrix(ini.cov.pars) | (length(nugget) > 1) | (length(kappa) > 1)) {
        if(messages.screen)
          cat(".variofit: searching for best initial value ...")
        ini.temp <- matrix(ini.cov.pars, ncol=2)
        grid.ini <- as.matrix(expand.grid(sigmasq=unique(ini.temp[,1]),
                                          phi=unique(ini.temp[,2]),
                                          tausq=unique(nugget), kappa=unique(kappa)))
        ##  loss function:
        v.loss <- function(parms, u, v, n, cov.model, weights){
          sigmasq <- parms[1]
          phi <- parms[2]
          if(cov.model == "power") phi <- 2 * exp(phi)/(1+exp(phi))
          tausq <- parms[3]
          kappa <- parms[4]
          if(cov.model == "power")
            v.mod <- tausq +
            .cov.spatial(u, cov.pars=c(sigmasq, phi), cov.model="power", kappa=kappa)
          else
            v.mod <- (sigmasq + tausq) -
            .cov.spatial(u, cov.pars=c(sigmasq, phi), cov.model = cov.model,
                        kappa = kappa)
          if(weights == "equal")
            loss <- sum((v - v.mod)^2)
          if (weights == "npairs") 
            loss <- sum(n * (v - v.mod)^2)
          if (weights == "cressie") 
            loss <- sum((n/(v.mod^2)) * (v - v.mod)^2)
          return(loss)
        }
        grid.loss <- apply(grid.ini, 1, v.loss, u=XY$u, v=XY$v, n=XY$n, cov.model = cov.model, weights = weights)
        ini.temp <- grid.ini[which(grid.loss == min(grid.loss))[1],, drop=FALSE]
        if(is.R()) rownames(ini.temp) <- "initial.value"
        if(messages.screen){
          cat(" selected values:\n")
          print(rbind(round(ini.temp, digits=2), status=ifelse(c(FALSE, FALSE, fix.nugget, fix.kappa), "fix", "est")))
          cat(paste("loss value:", min(grid.loss), "\n"))
        }
        names(ini.temp) <- NULL
        ini.cov.pars <- ini.temp[1:2]
        nugget <- ini.temp[3]
        kappa <- ini.temp[4]
        grid.ini <- NULL
      }
      ##
      ## checking for unreasonable initial values
      ##
      if(ini.cov.pars[1] > 2*vmax)
        warning("unreasonable initial value for sigmasq (too high)")
      if(ini.cov.pars[1] + nugget > 3*vmax)
        warning("unreasonable initial value for sigmasq + nugget (too high)")
      if(vario$output.type != "cloud"){
        if(ini.cov.pars[1] + nugget < 0.3*vmax)
          warning("unreasonable initial value for sigmasq + nugget (too low)")
      }
      if(nugget > 2*vmax)
        warning("unreasonable initial value for nugget (too high)")
      if(ini.cov.pars[2] > 1.5*umax)
        warning("unreasonable initial value for phi (too high)")
      ##
      ## transforming kappa for constraint minimisation
      ##
      if(!fix.kappa){
        if(cov.model == "powered.exponential")
          Tkappa.ini <- log(kappa/(2-kappa))
        else
          Tkappa.ini <- log(kappa)
      }
      ##
      ## minimisation using "nls"
      ##
      if (minimisation.function == "nls") {
        if(ini.cov.pars[2] == 0) ini.cov.pars <- max(XY$u)/10
        if(kappa == 0) kappa <- 0.5
        if(cov.model == "power")
          Tphi.ini <- log(ini.cov.pars[2]/(2-ini.cov.pars[2])) 
        else Tphi.ini <- log(ini.cov.pars[2])
        XY$cov.model <- cov.model
        ##
        if (fix.nugget) {
          XY$nugget <- as.vector(nugget)
          if(fix.kappa){
            XY$kappa <- as.vector(kappa)
            res <- nls((v-nugget) ~ matrix((1-.cov.spatial(u,cov.pars=c(1,exp(Tphi)),
                                                          cov.model=cov.model, kappa=kappa)),
                                           ncol=1),
                       start=list(Tphi=Tphi.ini), data=XY, algorithm="plinear", ...)
          }
          else{
            if(cov.model == "powered.exponential")
              res <- nls((v-nugget) ~ matrix((1-.cov.spatial(u,cov.pars=c(1,exp(Tphi)),
                                                            cov.model=cov.model,
                                                            kappa=(2*exp(Tkappa)/(1+exp(Tkappa))))),
                                             ncol=1),
                         start=list(Tphi=Tphi.ini, Tkappa = Tkappa.ini),
                         data=XY, algorithm="plinear", ...)
            else
              res <- nls((v-nugget) ~ matrix((1-.cov.spatial(u,cov.pars=c(1,exp(Tphi)),
                                                            cov.model=cov.model,
                                                            kappa=exp(Tkappa))), ncol=1),
                         start=list(Tphi=Tphi.ini, Tkappa = Tkappa.ini),
                         data=XY, algorithm="plinear", ...)       
            kappa <- exp(coef(res)["Tkappa"])
            names(kappa) <- NULL
          }
          cov.pars <- coef(res)[c(".lin", "Tphi")]
          names(cov.pars) <- NULL
        }
        else{
          if(fix.kappa){
            XY$kappa <- kappa
            res <- nls(v ~ cbind(1,(1- .cov.spatial(u, cov.pars=c(1,exp(Tphi)),
                                                   cov.model = cov.model, kappa=kappa))),
                       start=list(Tphi=Tphi.ini), algorithm="plinear", data=XY, ...)
          }
          else{
            if(cov.model == "powered.exponential")
              res <- nls(v ~ cbind(1, (1-.cov.spatial(u, cov.pars=c(1, exp(Tphi)),
                                                     cov.model = cov.model,
                                                     kappa=(2*exp(Tkappa)/(1+exp(Tkappa)))))),
                         start=list(Tphi=Tphi.ini, Tkappa = Tkappa.ini),
                         algorithm="plinear", data=XY, ...)
            else
              res <- nls(v ~ cbind(1, (1-.cov.spatial(u, cov.pars=c(1, exp(Tphi)),
                                                     cov.model = cov.model,
                                                     kappa=exp(Tkappa)))),
                         start=list(Tphi=Tphi.ini, Tkappa = Tkappa.ini),
                         algorithm="plinear", data=XY, ...)
            kappa <- exp(coef(res)["Tkappa"]);names(kappa) <- NULL
          }
          nugget <- coef(res)[".lin1"];names(nugget) <- NULL
          cov.pars <- coef(res)[c(".lin2", "Tphi")]
          names(cov.pars) <- NULL
        }
        if(cov.model == "power")
          cov.pars[2] <- 2 * exp(cov.pars[2])/(1+exp(cov.pars[2]))  
        else cov.pars[2] <- exp(cov.pars[2])
        if(nugget < 0 | cov.pars[1] < 0){
          warning("\n.variofit: negative variance parameter found using the default option \"nls\".\n        Try another minimisation function and/or fix some of the parameters.\n")
          temp <- c(sigmasq=cov.pars[1], phi=cov.pars[2], tausq=nugget, kappa=kappa)
          print(rbind(round(temp, digits=4),
                      status=ifelse(c(FALSE, FALSE, fix.nugget, fix.kappa), "fix", "est")))
          return(invisible())
        }
        value <- sum(resid(res)^2)
        message <- "nls does not provides convergence message"
      }
      ##
      ## minimisation using "optim" or "nlm"
      ##
      if (minimisation.function == "nlm" | minimisation.function == "optim") {
        ##
        ## Preparing lists for the minimiser
        ##
        .global.list <- list(u = XY$u, v = XY$v, n=XY$n, fix.nugget = fix.nugget,
                             nugget = nugget, fix.kappa = fix.kappa, kappa = kappa,
                             cov.model = cov.model, m.f = minimisation.function,
                             weights = weights)
        ##
        ## Preparing initial value
        ##
        ini <- ini.cov.pars
        if(cov.model == "power") ini[2] <- log(ini[2]/(2-ini[2])) 
        if(cov.model == "linear") ini <- ini[1] 
        if(fix.nugget == FALSE) ini <- c(ini, nugget)
        ## setting kappa > 0 for both methods
        if(!fix.kappa) ini <- c(ini, Tkappa.ini)
        names(ini) <- NULL
        # Pedro Delicado, January 29th, 2020:
        # disabling the "nlm" option because the sentences
        #    assign(".temp.theta",  NULL, pos=1)
        #    assign(".temp.theta", theta, pos=1)
        # do not pass the CRAN checking
        # if(minimisation.function == "nlm"){
        #   result <- nlm(.loss.vario, ini, g.l = .global.list, ...)
        #   result$par <- result$estimate
        #   result$value <- result$minimum
        #   result$convergence <- result$code
        #   if(!is.null(get(".temp.theta", pos =1)))
        #     result$par <- get(".temp.theta", pos=1)
        # }
        # else{
        {
          #        if(fix.kappa == FALSE) ini <- c(ini, kappa)
          #        names(ini) <- NULL
          lower.l <- sapply(limits, function(x) x[1])
          upper.l <- sapply(limits, function(x) x[2])
          if(fix.kappa == FALSE){
            if(fix.nugget){
              lower <- lower.l[c("sigmasq.lower", "phi.lower","kappa.lower")]
              upper <- upper.l[c("sigmasq.upper", "phi.upper","kappa.upper")]
            }
            else{
              lower <- lower.l[c("sigmasq.lower", "phi.lower",
                                 "tausq.rel.lower", "kappa.lower")]
              upper <- upper.l[c("sigmasq.upper", "phi.upper",
                                 "tausq.rel.upper", "kappa.upper")]
            }
          }
          else{
            if(cov.model == "power"){
              if(fix.nugget){
                lower <- lower.l[c("sigmasq.lower", "phi.lower")]
                upper <- upper.l[c("sigmasq.upper", "phi.upper")]
              }
              else{
                lower <- lower.l[c("sigmasq.lower", "phi.lower", "tausq.rel.lower")]
                upper <- upper.l[c("sigmasq.upper", "phi.upper", "tausq.rel.upper")]
              }
            }
            else{
              lower <- lower.l["phi.lower"]
              upper <- upper.l["phi.upper"]
            }
          }
          result <- optim(ini, .loss.vario, method = "L-BFGS-B",
                          hessian = TRUE, lower = lower,
                          upper = upper, g.l = .global.list, ...)
          #        require(methods)
          #        if(exists("trySilent"))
          #          hess <- trySilent(solve(as.matrix(result$hessian)))
          #        else{
          #          op.sem <- options()$show.error.messages
          #          options(show.error.messages = FALSE)
          #          hess <- try(solve(as.matrix(result$hessian)))
          #          options(show.error.messages = op.sem)
          #        }
          #        if(!inherits(hess, "try-error")) hess <- sqrt(diag(hess))
          #        else print("WARNING: unable to compute the hessian")
        }
        value <- result$value
        message <- paste(minimisation.function, "convergence code:", result$convergence)
        if(cov.model == "linear")
          result$par <- c(result$par[1],1,result$par[-1])
        cov.pars <- as.vector(result$par[1:2])
        if(cov.model == "power")
          cov.pars[2] <- 2 * exp(cov.pars[2])/(1+exp(cov.pars[2]))  
        if(!fix.kappa){
          if (fix.nugget)
            kappa <- result$par[3]
          else{
            nugget <- result$par[3]
            kappa <- result$par[4]
          }
          ## kappa now > 0 for both nlm() and optim()
          ##        if(minimisation.function == "nlm"){
          if(.global.list$cov.model == "powered.exponential")
            kappa <- 2*(exp(kappa))/(1+exp(kappa))
          else kappa <- exp(kappa)
          ##        }
        }
        else
          if(!fix.nugget)
            nugget <- result$par[3]        
      }
    }
    ##
    ## Estimating implicity beta
    ##
    
    ##
    ## Preparing output
    ##
    estimation <- list(nugget = nugget, cov.pars = cov.pars, 
                       cov.model = cov.model, kappa = kappa, value = value, 
                       trend = vario$trend, beta.ols = vario$beta.ols,
                       practicalRange = practicalRange(cov.model=cov.model,
                                                       phi = cov.pars[2], kappa = kappa),
                       max.dist = max.dist, 
                       minimisation.function = minimisation.function)
    #  if(exists("hess")) estimation$hessian <- hess
    estimation$weights <- weights
    if(weights == "equal") estimation$method <- "OLS"
    else estimation$method <- "WLS"
    estimation$fix.nugget <- fix.nugget
    estimation$fix.kappa <- fix.kappa
    estimation$lambda <- vario$lambda
    estimation$message <- message
    estimation$call <- call.fc
    oldClass(estimation) <- c("variomodel", "variofit")
    return(estimation)
  }

.loss.vario <- # copied from geoR in January 2020, because geoR was close to be uncontinued
  function (theta, g.l) 
  {
    if(g.l$cov.model == "linear")
      theta <- c(theta[1], 1, theta[-1])
    # Pedro Delicado, January 29th, 2020:
    # disabling the "nlm" option because the sentences
    #    assign(".temp.theta",  NULL, pos=1)
    #    assign(".temp.theta", theta, pos=1)
    # do not pass the CRAN checking
    # ##
    # ## Imposing constraints for nlm
    # ##
    # if(g.l$m.f == "nlm"){  
    #   assign(".temp.theta",  NULL, pos=1)
    #   if(!g.l$fix.kappa){
    #     if(g.l$fix.nugget){
    #       if(g.l$cov.model == "power")
    #         theta.minimiser <- theta[1]
    #       else          
    #         theta.minimiser <- theta[1:2]
    #       Tkappa <- theta[3]
    #     }
    #     else{
    #       if(g.l$cov.model == "power")
    #         theta.minimiser <- theta[c(1:3)]
    #       else          
    #         theta.minimiser <- theta[1:3]
    #       Tkappa <- theta[4]
    #     }
    #   }
    #   else theta.minimiser <- theta
    #   penalty <- 10000 * sum(0 - pmin(theta.minimiser, 0))
    #   theta <- pmax(theta.minimiser, 0)
    #   if(!g.l$fix.kappa) theta <- c(theta.minimiser, Tkappa)
    #   if (any(theta.minimiser < 0)) assign(".temp.theta", theta, pos=1)
    #   else penalty <- 0
    # }
    # else penalty <- 0
    penalty <- 0 
    ##
    ## reading parameters
    ##
    if(!g.l$fix.kappa){
      if (g.l$fix.nugget){
        tausq <- g.l$nugget
        Tkappa <- theta[3]
      }
      else{
        tausq <- theta[3]
        Tkappa <- theta[4]
      }
      ## kappa now > 0 for both nlm() and optim()
      ##if(g.l$m.f == "nlm"){
      if(g.l$cov.model == "powered.exponential")
        kappa <-  2*(exp(Tkappa))/(1+exp(Tkappa))
      else kappa <- exp(Tkappa)
      ##}
      ##else kappa <- Tkappa
    }
    else{
      kappa <- g.l$kappa
      if (g.l$fix.nugget) tausq <- g.l$nugget
      else tausq <- theta[3]
    }
    ##
    sigmasq <- theta[1]
    phi <- theta[2]
    if(g.l$cov.model == "power") phi <- 2 * exp(phi)/(1+exp(phi))
    sill.total <- sigmasq + tausq
    ##
    ## Computing values for the theoretical variogram 
    ##
    if(any(g.l$cov.model == c("linear", "power")))
      gammaU <- tausq + sigmasq * (g.l$u^phi)
    else
      gammaU <- sill.total - .cov.spatial(g.l$u, cov.model = g.l$cov.model, 
                                         kappa = kappa, cov.pars = c(sigmasq, phi))
    ##
    ## Computing loss function
    ##
    if(g.l$weight == "equal")
      loss <- sum((g.l$v - gammaU)^2)
    if (g.l$weights == "npairs") 
      loss <- sum(g.l$n * (g.l$v - gammaU)^2)
    if (g.l$weights == "cressie") 
      loss <- sum((g.l$n/(gammaU^2)) * (g.l$v - gammaU)^2)
    if(loss > (.Machine$double.xmax^0.5) | loss == Inf | loss == -Inf | is.nan(loss))
      loss <- .Machine$double.xmax^0.5
    return(loss + penalty)
  }

.geoR.cov.models <-
   c("matern", "exponential", "gaussian", "spherical",
     "circular", "cubic", "wave", "linear", "power",
     "powered.exponential", "stable", "cauchy", "gencauchy",
     "gneiting", "gneiting.matern", "pure.nugget")

"geoRCovModels" <- .geoR.cov.models

"practicalRange" <-
  function (cov.model, phi, kappa=0.5, correlation = 0.05, ...) 
  {
    cov.model <- match.arg(cov.model, choices = .geoR.cov.models)
    .check.cov.model(cov.model = cov.model, cov.pars=c(1,phi), kappa=kappa, output=FALSE)
    if(cov.model %in% c("circular","cubic","spherical"))
      return(phi)
    if(any(cov.model %in% c("pure.nugget")))
      return(0)  
    if(any(cov.model %in% c("linear")))
      return(Inf)  
    if(any(cov.model %in% c("power")))
      return(Inf)  
    findRange <- function(range, cm, p, k, cor)
      .cov.spatial(range, cov.model = cm, kappa = k, cov.pars = c(1, p))-cor
    pr <- uniroot(findRange, interval = c(0, 50 * phi + 1), 
                  cm = cov.model, p = phi, k = kappa, cor = correlation, 
                  ...)$root
    return(pr)
  }

".check.cov.model" <-
  function(cov.model, cov.pars, kappa, env=NULL, output=TRUE)
  {
    ## extracting covariance parameters
    if(is.vector(cov.pars)) sigmasq <- cov.pars[1]
    else sigmasq <- cov.pars[, 1]
    if(is.vector(cov.pars)) phi <-  cov.pars[2]
    else phi <- cov.pars[, 2]
    if(missing(kappa) || is.null(kappa)) kappa <- NA
    ## checking for nested models
    cov.pars <- drop(cov.pars)
    if(is.vector(cov.pars)) ns <- 1
    else{
      ns <- nrow(cov.pars)
      if(length(cov.model) == 1) cov.model <- rep(cov.model, ns)
      if(length(kappa) == 1) kappa <- rep(kappa, ns)
    }
    if(length(cov.model) != ns) stop("wrong length for cov.model")
    ##
    cov.model <- sapply(cov.model, match.arg, .geoR.cov.models)
    cov.model[cov.model == "stable"] <- "powered.exponential"
    if(any(cov.model == c("gneiting.matern", "gencauchy"))){
      if(length(kappa) != 2*ns)
        stop(paste("wrong length for kappa, ", cov.model, "model requires two values for the argument kappa")) 
    }
    else{
      if(length(kappa) != ns) stop('wrong length for kappa')
    }
    ## settings for power model (do not reverse order of the next two lines!)
    phi[cov.model == "linear"] <- 1
    cov.model[cov.model == "linear"] <- "power"
    ## checking input for cov. models with extra parameter(s)
    if(any(cov.model == 'gneiting.matern') && ns > 1)
      stop('nested models including the gneiting.matern are not implemented') 
    for(i in 1:ns){
      if(any(cov.model[i] == c("matern","powered.exponential", "cauchy",
                               "gneiting.matern", "gencauchy"))){
        if(any(cov.model[i] == c("gneiting.matern", "gencauchy"))){
          if(any(is.na(kappa)) || length(kappa) != 2*ns)
            stop(paste(cov.model[i],"correlation function model requires a vector with 2 parameters in the argument kappa"))
        }
        else{
          if(is.na(kappa[i]) | is.null(kappa[i]))
            stop("for matern, powered.exponential and cauchy covariance functions the parameter kappa must be provided")
        }
        if((cov.model[i] == "matern" | cov.model[i] == "powered.exponential" | 
            cov.model[i] == "cauchy") & length(kappa) != 1*ns)
          stop("kappa must have 1 parameter for this correlation function")
        if(cov.model[i] == "matern" & kappa[i] == 0.5) cov.model[i] == "exponential"
      }      
      if(cov.model[i] == "power")
        if(any(phi[i] >= 2) | any(phi[i] <= 0))
          stop("for power model the phi parameters must be in the interval ]0,2[")
    }
    if(!is.null(env)){
      assign("sigmasq", sigmasq, envir=env)
      assign("phi", phi, envir=env)
      assign("kappa", kappa, envir=env)
      assign("ns", ns, envir=env)
      assign("cov.model", cov.model, envir=env)
    }
    if(output)
      return(list(cov.model=cov.model, sigmasq=sigmasq, phi=phi, kappa=kappa, ns=ns))
    else return(invisible())
  }

"matern" <-
  function (u, phi, kappa) 
  {
    if(is.vector(u)) names(u) <- NULL
    if(is.matrix(u)) dimnames(u) <- list(NULL, NULL)
    uphi <- u/phi
    uphi <- ifelse(u > 0,
                   (((2^(-(kappa-1)))/ifelse(0, Inf,gamma(kappa))) *
                      (uphi^kappa) *
                      besselK(x=uphi, nu=kappa)), 1)    
    uphi[u > 600*phi] <- 0 
    return(uphi)
}

pars.limits <-
 function (phi = c(lower = 0, upper = +Inf), 
          sigmasq = c(lower = 0,upper = +Inf), 
          nugget.rel = c(lower = 0, upper = +Inf), 
          kappa = c(lower = 0, upper = +Inf), 
          kappa2 = c(lower = 0, upper = +Inf), 
          lambda = c(lower = -3, upper = 3), 
          psiR = c(lower = 1, upper = +Inf), 
          psiA = c(lower = 0, upper = 2 * pi), 
          tausq.rel = nugget.rel) 
{
  if (length(phi) != 2) 
    stop("phi must be a 2 components vector with lower and upper limits for the parameter phi")
  if (length(sigmasq) != 2) 
    stop("sigmasq must be a 2 components vector with lower and upper limits for the parameter sigmasq")
  if (length(tausq.rel) != 2) 
    stop("tausq.rel must be a 2 components vector with lower and upper limits for the parameter tausq.rel")
  if (length(kappa) != 2) 
    stop("kappa must be a 2 components vector with lower and upper limits for the parameter kappa")
  if (length(kappa2) != 2) 
    stop("kappa must be a 2 components vector with lower and upper limits for the parameter kappa")
  if (length(lambda) != 2) 
    stop("lambda must be a 2 components vector with lower and upper limits for the parameter lambda")
  if (length(psiR) != 2) 
    stop("psiR must be a 2 components vector with lower and upper limits for the parameter psiR")
  if (length(psiA) != 2) 
    stop("psiA must be a 2 components vector with lower and upper limits for the parameter psiA")
  if (phi[1] >= phi[2]) 
    stop("parameter phi: lower limit greater or equal upper limit")
  if (sigmasq[1] >= sigmasq[2]) 
    stop("parameter sigmasq: lower limit greater or equal upper limit")
  if (tausq.rel[1] >= tausq.rel[2]) 
    stop("parameter tausq.rel: lower limit greater or equal upper limit")
  if (kappa[1] >= kappa[2]) 
    stop("parameter kappa: lower limit greater or equal upper limit")
  if (kappa2[1] >= kappa2[2]) 
    stop("parameter kappa: lower limit greater or equal upper limit")
  if (lambda[1] >= lambda[2]) 
    stop("parameter lambda: lower limit greater or equal upper limit")
  if (psiR[1] >= psiR[2]) 
    stop("parameter psiR: lower limit greater or equal upper limit")
  if (psiA[1] >= psiA[2]) 
    stop("parameter psiA: lower limit greater or equal upper limit")
  names(phi) <- names(sigmasq) <- names(tausq.rel) <- names(kappa) <- names(kappa2) <- names(lambda) <- names(psiR) <- names(psiA) <- c("lower", 
                                                                                                                                        "upper")
  return(list(phi = phi, sigmasq = sigmasq, tausq.rel = tausq.rel, 
              kappa = kappa, kappa2 = kappa2, lambda = lambda, psiR = psiR, 
              psiA = psiA))
 }

"plot.variogram" <-
  function (x, max.dist, vario.col = "all", scaled = FALSE,  
            var.lines = FALSE,  envelope.obj = NULL,
            pts.range.cex, bin.cloud = FALSE,  ...) 
  {
    fc.call <- match.call()
    fc.call$max.dist <- fc.call$vario.col <- fc.call$scaled <- 
      fc.call$pts.range.cex <-  fc.call$bin.cloud <-
      fc.call$envelope.obj <- NULL
    Ldots <- list(...)
    argsnames <- c(names(formals(plot.default)), names(par()))
    names(Ldots) <- argsnames[charmatch(names(Ldots), argsnames)]
    if(missing(max.dist)) max.dist <- max(x$u)
    if(is.null(Ldots$xlim)) fc.call$xlim <- c(0, max.dist) 
    if (bin.cloud == TRUE && all(is.na(x$bin.cloud))) 
      stop("plot.variogram: object must be a binned variogram with option bin.cloud=TRUE")
    if (bin.cloud == TRUE && any(!is.na(x$bin.cloud))) 
      boxplot(x$bin.cloud, varwidth = TRUE, 
              xlab = "distance", ylab = paste(x$estimator.type, "variogram"))
    else {
      if(!missing(pts.range.cex)){
        cex.min <- min(pts.range.cex)
        cex.max <- max(pts.range.cex)
        if(cex.min != cex.max){
          pts.prop <- TRUE
          sqn <- sqrt(x$n[x$u <= max.dist])
          pts.cex <- cex.min + ((sqn - min(sqn)) * (cex.max - cex.min) / (max(sqn) - min(sqn)))
        }
        else pts.prop <- FALSE
      }
      else pts.prop <- FALSE 
      u <- x$u[x$u <= max.dist]
      v <- x$v
      if(is.vector(v) | length(v) == length(x$u))
        v <- matrix(v, ncol=1)
      v <- v[x$u <= max.dist,, drop=FALSE]
      if(vario.col == "all")
        vario.col <- 1:dim(v)[2]
      else
        if(mode(vario.col) != "numeric" | any(vario.col > ncol(v)))
          stop("argument vario.col must be equals to \"all\" or a vector indicating the column numbers to be plotted")
      v <- v[, vario.col, drop=FALSE]
      if (scaled) v <- t(t(v)/x$var.mark[vario.col])
      if(is.null(list(...)$ylim)){ 
        ymax <- max(v)
        if (!is.null(envelope.obj)) 
          ymax <- max(c(envelope.obj$v.upper, ymax))
        fc.call$ylim <- c(0, ymax)
      }
      if(ncol(v) == 1){
        fc.call[[1]] <- as.name("plot")
        fc.call$x <- data.frame(distance=u, semivariance = as.vector(v))
        if(pts.prop) fc.call$cex <- pts.cex
        eval(fc.call, sys.frame(sys.parent()))
      }
      else{
        fc.call[[1]] <- as.name("matplot")
        fc.call$x <- u
        fc.call$y <- v
        if(is.null(Ldots$xlab)) fc.call$xlab <- "distance"
        if(is.null(Ldots$ylab)) fc.call$ylab <- "semivariance"
        if(is.null(Ldots$ty)){
          if (x$output.type == "bin") fc.call$type <- "p"
          if (x$output.type == "smooth") fc.call$type <- "l"
          if (x$output.type == "cloud") fc.call$type <- "p"
        }
        #      if(is.null(Ldots$col)) fc.call$col <- 1:6
        #      if(is.null(Ldots$lty)) fc.call$lty <- 1:5
        #      if(is.null(Ldots$lwd)) $lwd <- 1
        #      if(is.null(Ldots$pch)) Ldots$pch <- NULL
        #      if(is.null(Ldots$cex)) Ldots$cex <- NULL
        #      if(is.null(Ldots$add)) Ldots$add <- FALSE
        #      
        #      matplot(x=u, y= v, xlim = c(0, max.dist), ylim = Ldots$ylim, 
        #              xlab = Ldots$xlab, ylab = Ldots$ylab, type = Ldots$type,
        #              add = Ldots$add, pch = Ldots$pch,
        #              lty = Ldots$lty, lwd = Ldots$lwd, col = Ldots$col)
        eval(fc.call, sys.frame(sys.parent()))
      }
      if (var.lines) {
        if (scaled) abline(h = 1, lty = 3)
        else abline(h = x$var.mark, lty = 3)
      }
      if (!is.null(envelope.obj)) {
        lines(u, envelope.obj$v.lower, lty = 4)
        lines(u, envelope.obj$v.upper, lty = 4)
      }
    }
    return(invisible())
  }

"lines.variomodel.variofit" <-
#  "lines.variomodel.likGRF" <-
  function (x, max.dist, scaled = FALSE, ...)
  {
    my.l <- list()
    if(missing(max.dist)){
      my.l$max.dist <- x$max.dist
      if (is.null(my.l$max.dist)) 
        stop("argument max.dist needed for this object")
    }
    else
      my.l$max.dist <- max.dist
    if (any(x$cov.model == c("matern","powered.exponential",
                             "cauchy", "gencauchy", "gneiting.matern"))) 
      my.l$kappa <- x$kappa
    else kappa <- NULL
    if (is.vector(x$cov.pars)) 
      my.l$sill.total <- x$nugget + x$cov.pars[1]
    else my.l$sill.total <- x$nugget + sum(x$cov.pars[, 1])
    my.l$nugget <- x$nugget
    my.l$cov.pars <- x$cov.pars
    my.l$cov.model <- x$cov.model
    if (scaled){
      if(is.vector(x$cov.model))
        my.l$cov.pars[1] <-  my.l$cov.pars[1]/my.l$sill.total
      else my.l$cov.pars[,1] <-  my.l$cov.cov.pars[,1]/my.l$sill.total
      my.l$sill.total <- 1
    }
    gamma.f <- function(x, my.l){
      if(any(my.l$cov.model == c("linear", "power")))
        return(my.l$nugget + my.l$cov.pars[1] * (x^my.l$cov.pars[2]))
      else
        return(my.l$sill.total -
                 .cov.spatial(x, cov.model = my.l$cov.model,
                             kappa = my.l$kappa,
                             cov.pars = my.l$cov.pars))
    }
    curve(gamma.f(x,my.l=my.l), from=0, to=my.l$max.dist, add=TRUE, ...)
    return(invisible())
  }

