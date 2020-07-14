setMethodS3("getSigmamin", "NPVI", function(this, ...) {
  this$.sigmamin;
})

setMethodS3("getSigmamax", "NPVI", function(this, ...) {
  this$.sigmamax;
})

setMethodS3("getSigma", "NPVI", function(this, tabulate, ...) {
  if (missing(tabulate)) {
    tabulate <- getTabulate(this);
  }
  if (!tabulate) {
    this$.sigma;
  } else {
    this$.sigmatab;
  }
})

setMethodS3("getSigmaAux", "NPVI", function(this, tabulate, ...) {
  if (missing(tabulate)) {
    tabulate <- getTabulate(this);
  }
  if (!tabulate) {
    this$.SigmaAux;
  } else {
    this$.SigmaAuxtab;
  } 
})


setMethodS3("setSigma", "NPVI", function(this, sigma, ...) {
  ## Argument 'sigma':
  if ((!is.null(sigma))  && (mode(sigma)!="function")) {
    throw("Argument \var{sigma} should be of mode 'function', not ", mode(sigma));
  }
  
  sigmamin <- getSigmamin(this)
  sigmamax <- getSigmamax(this)
  thresholdedSigma <- function(W) {
    threshold(sigma(W), min=sigmamin, max=sigmamax)
  }
  this$.sigma <- thresholdedSigma ;
})

setMethodS3("setSigmaAux", "NPVI", function(this, sigmaAux, ...) {
  ## Argument 'sigmaAux':
  if ((!is.null(sigmaAux))  && (mode(sigmaAux)!="function")) {
    throw("Argument \var{sigmaAux} should be of mode 'function', not ", mode(sigmaAux));
  }
  
  sigmamin <- getSigmamin(this)
  sigmamax <- getSigmamax(this)
  thresholdedSigmaAux <- function(W) {
    threshold(sigmaAux(W), min=sigmamin, max=sigmamax)
  }
  this$.sigmaAux <- thresholdedSigmaAux ;
})


setMethodS3("setSigmaTab", "NPVI", function(this, sigma, ...) {
  ## Argument 'sigma':
  if ((!is.null(sigma))  && (mode(sigma)!="function")) {
    throw("Argument \var{sigma} should be of mode 'function', not ", mode(sigma));
  }

  sigmamin <- getSigmamin(this)
  sigmamax <- getSigmamax(this)
  thresholdedSigma <- function(W) {
    threshold(sigma(W), min=sigmamin, max=sigmamax)
  }

  this$.sigmatab <- thresholdedSigma
})

setMethodS3("setSigmaAuxTab", "NPVI", function(this, sigmaAux, ...) {
  ## Argument 'sigmaAux':
  if ((!is.null(sigmaAux))  && (mode(sigmaAux)!="function")) {
    throw("Argument \var{sigmaAux} should be of mode 'function', not ", mode(sigmaAux));
  }

  sigmamin <- getSigmamin(this)
  sigmamax <- getSigmamax(this)
  thresholdedSigmaAux <- function(W) {
    threshold(sigmaAux(W), min=sigmamin, max=sigmamax)
  }

  this$.SigmaAuxtab <- thresholdedSigmaAux
})


setMethodS3("initializeSigma", "NPVI", function(this, sigmaAux, g, ...) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'sigmaAux':
  if (mode(sigmaAux) != "function") {
    throw("Argument 'sigmaAux' should be a function, not a ", mode(sigmaAux));
  }
  ## Argument 'g':
  if (mode(g) != "function") {
    throw("Argument 'g' should be a function, not a ", mode(g));
  }

  setSigmaAux(this, sigmaAux)

  ##sigma
  sigma <- function(W) {
    sigmaAux(W)* as.vector((1-g(W)))
  }
  setSigma(this, sigma)

  ## tabulated version of 'sigmaAux'
  fW <- getFW(this);
  obs <- getObs(this);  
  SIGMAAUXTAB <- sigmaAux(fW(obs)); ## a *vector*, not a function
  
  
  sigmaAuxtab <- function(ii) {
    if(is.matrix(SIGMAAUXTAB)){
      SIGMAAUXTAB[ii,];
    }
    else{
      SIGMAAUXTAB[ii];
    }
  }
  #<\Cabral>
  setSigmaAuxTab(this, sigmaAuxtab)
  
  ## tabulated version of 'sigma'
  SIGMATAB <- sigma(fW(obs)); ## a *vector*, not a function
  
  sigmatab <- function(ii) {    
    if(is.matrix(SIGMATAB)){
      SIGMATAB[ii,]
    }
    else{
      SIGMATAB[ii];
    }
  }
  
  setSigmaTab(this, sigmatab)
  
  ## calculating the value of sigma2
  derivH_beta <- getDerivH_beta(this)
  weights <-getWeightsW(this)
  
  ## To do : to remove at a later stage.... just ensuring that the weights are all set to 1
  #weights <- rep(1, length(weights)) 
  
  
  W <- fW(obs)
  derivH_betaW <- derivH_beta(W) 
  
  H.list <- lapply(1:nrow(obs), function(ii) {
    row <- matrix(as.numeric(derivH_betaW[ii,]), ncol=1)
    mat <- ( row %*% t(row) ) * weights[ii] * SIGMATAB[ii]
  })
  
  
  sigma2 <- Reduce("+", H.list)/sum(weights)
  setSigma2(this, sigma2)
})

setMethodS3("updateSigma", "NPVI", function(this, dev, exact=TRUE, effICW, ...) {
  updateSigmaNonTab(this, dev, exact=exact, effICW, ...)
  updateSigmaTab(this, dev, exact=exact, effICW, ...)
})

#<Cabral>
setMethodS3("updateSigmaNonTab", "NPVI", function(this, dev, exact=TRUE, effICW, ...) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'dev':
  if (mode(dev) != "function") {
    throw("Argument 'dev' should be a function, not a ", mode(dev));
  }
  ## Argument 'exact':
  exact <- Arguments$getLogical(exact);
  ## Argument 'effICW':
  if (exact) {
    if (missing(effICW)) {
      throw("Argument 'effICW' is required when 'exact' is TRUE");
    }
    if (mode(effICW) != "function") {
      throw("Argument 'dev' should be a function, not a ", mode(dev));
    }
  }

  fW <- getFW(this)
  eps <- getEpsilon(this)
  eps <- matrix( as.numeric(eps), ncol = 1) 
  
  sigma <- getSigma(this, tabulate=FALSE);
  g <- getG(this, tabulate=FALSE);
  
  if (!exact) { ## if do not use exact expression
    sigma1 <- function(W) {
        sigma(W) + t( eps %*% dev(W) ); 
    }
  } else { ## if use exact expression
    sigma1 <- function(W) {
      sigmaW <- sigma(W);
      theEffICW <- effICW(W)
      
      numerator <- sigmaW +   ( devW + sigmaW * theEffICW)  %*% eps ;         
      denominator <-  1 + theEffICW  %*% eps ;
        
      numerator/denominator;
    }
  }
  
  
  sigmaAux1 <- function(W) {
    sigma1(W)/(1-g(W))
  }
  setSigmaAux(this, sigmaAux1);
  setSigma(this, sigma1);
})


setMethodS3("updateSigmaTab", "NPVI", function(this, dev, exact=TRUE, effICW, ...) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'dev':
  if (mode(dev) != "function") {
    throw("Argument 'dev' should be a function, not a ", mode(dev));
  }
  ## Argument 'exact':
  exact <- Arguments$getLogical(exact);
  
  ## Argument 'effICW':
  if (exact) {
    if (missing(effICW)) {
      throw("Argument 'effICW' is required when 'exact' is TRUE");
    }
    if (mode(effICW) != "function") {
      throw("Argument 'dev' should be a function, not a ", mode(dev));
    }
  }

  fW <- getFW(this) 
  eps <- getEpsilon(this)
  eps <- matrix(as.numeric(eps), ncol=1)
  
  sigma <- getSigma(this, tabulate=TRUE)
  g <- getG(this, tabulate=TRUE)
  obs <- getObs(this, tabulate=TRUE);
  sigmaW <- sigma(obs[, "W"])
  gW <- g(obs[, "W"])
  
  rm(sigma, g, obs)

  obs <- getObs(this)
  devW <- dev(fW(obs))
  W <- obs[, "W"]
  
  
  if (!exact) { ## do not use exact expression
    sigma1W <- sigmaW + t(eps %*%  devW ) ;
  } else { ## use exact expression
    theEffICW <- effICW(W)
    numerator <- sigmaW + ( devW + sigmaW*theEffICW) %*% eps  ;    
    denominator <-  1 +  theEffICW %*% eps ;
      
    sigma1W <- numerator/denominator;
  }

  sigmaAux1W <- sigma1W/(1-gW)
  sigmaAux1tab <- function(ii) {
    if(is.matrix(sigmaAux1W)){
      sigmaAux1W[ii,]
    }
    else{
      sigmaAux1W[ii];
    }        
  }
  setSigmaAuxTab(this, sigmaAux1tab)

  sigma1tab <- function(ii) {
    if(is.matrix(sigma1W)){
      sigma1W[ii,]
    }
    else{
      sigma1W[ii];
    }    
  }
  setSigmaTab(this, sigma1tab)
  
  ## updating the value of the sigma2
  derivH_beta <- getDerivH_beta(this)
  weights <-getWeightsW(this)
  
  ## To do : to remove at a later stage.... just ensuring that the weights are all set to 1
  #weights <- rep(1, length(weights)) 
    
  W <- fW(obs)
  
  derivH_betaW <- derivH_beta(W) 
  
  H.list <- lapply(1:nrow(obs), function(ii){
    row <- matrix(as.numeric(derivH_betaW[ii,]), ncol=1)
    #mat <- ( row %*% t(row) ) * weights[ii] * sigma1W[ii]
    mat <- ( row %*% t(row) ) * sigma1W[ii]
  })
  
  #sigma2 <- Reduce("+", H.list)/sum(weights)
  sigma2 <- Reduce("+", H.list)/length(H.list)
  setSigma2(this, sigma2)
})

############################################################################
## HISTORY:
## 2015-09-10
## o Created.
############################################################################

