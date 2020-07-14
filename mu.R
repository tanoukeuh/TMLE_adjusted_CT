setMethodS3("getMumin", "NPVI", function(this, ...) {
  this$.mumin;
})

setMethodS3("getMumax", "NPVI", function(this, ...) {
  this$.mumax;
})

setMethodS3("getMu", "NPVI", function(this, tabulate, ...) {
  if (missing(tabulate)) {
    tabulate <- getTabulate(this);
  }
  if (!tabulate) {
    this$.mu;
  } else {
    this$.mutab;
  } 
})

setMethodS3("getPureMu", "NPVI", function(this, tabulate, ...) {
  if (missing(tabulate)) {
    tabulate <- getTabulate(this);
  }
  if (!tabulate) {
    this$.puremu;
  } else {
    this$.puremutab;
  } 
})

setMethodS3("getMuAux", "NPVI", function(this, tabulate, ...) {
  if (missing(tabulate)) {
    tabulate <- getTabulate(this);
  }
  if (!tabulate) {
    this$.muAux;
  } else {
    this$.muAuxtab;
  } 
})

#<Cabral>
setMethodS3("setMu", "NPVI", function(this, mu, ...) {
  ## Argument 'mu':
  if ((!is.null(mu))  && (mode(mu)!="function")) {
    throw("Argument \var{mu} should be of mode 'function', not ", mode(mu));
  }
  
  mumin <- getMumin(this)
  mumax <- getMumax(this)
  
  isBetaMult <- getIsBetaMult(this)
  derivH_beta <- getDerivH_beta(this)
  
  thresholdedMu <- function(W) {
    if(isBetaMult){
      result <- threshold(mu(W), min=mumin, max=mumax) * derivH_beta(fW(cbind(W=W,X=NA,Y=NA)))
    }
    else{
      result <- threshold(mu(W), min=mumin, max=mumax)
    }
    result
  }
  
  this$.mu <- thresholdedMu;
  
  if(isBetaMult){
    thresholdedPureMu <- function(W) {
      threshold(mu(W), min=mumin, max=mumax)
    }
    
    this$.puremu <- thresholdedPureMu;
  }
})
#<\Cabral>

setMethodS3("setMuAux", "NPVI", function(this, muAux, ...) {
  ## Argument 'muAux':
  if ((!is.null(muAux))  && (mode(muAux)!="function")) {
    throw("Argument \var{muAux} should be of mode 'function', not ", mode(muAux));
  }
  
  mumin <- getMumin(this)
  mumax <- getMumax(this)
  thresholdedMuAux <- function(W) {
    threshold(muAux(W), min=mumin, max=mumax)
  }
  this$.muAux <- thresholdedMuAux ;
})


setMethodS3("setMuTab", "NPVI", function(this, mu, ...) {
  ## Argument 'mu':
  if ((!is.null(mu))  && (mode(mu)!="function")) {
    throw("Argument \var{mu} should be of mode 'function', not ", mode(mu));
  }

  mumin <- getMumin(this)
  mumax <- getMumax(this)
  
  isBetaMult <- getIsBetaMult(this)
  derivH_betaW <- getDerivH_beta(this)
  fW <- getFW(this)
  
  thresholdedMu <- function(W){
    if(isBetaMult){
      result <- threshold(mu(W), min=mumin, max=mumax) * derivH_betaW(fW(cbind(W=W,X=NA,Y=NA)))
    }
    else{
      result <- threshold(mu(W), min=mumin, max=mumax)  
    }
    result
  }
  this$.mutab <- thresholdedMu
  
  if(isBetaMult){
    thresholdedPureMu <- function(W) {
      threshold(mu(W), min=mumin, max=mumax)
    }
  
    this$.puremutab <- thresholdedPureMu
  }
  
})

setMethodS3("setMuAuxTab", "NPVI", function(this, muAux, ...) {
  ## Argument 'muAux':
  if ((!is.null(muAux))  && (mode(muAux)!="function")) {
    throw("Argument \var{muAux} should be of mode 'function', not ", mode(muAux));
  }

  mumin <- getMumin(this)
  mumax <- getMumax(this)
  thresholdedMuAux <- function(W) {
    threshold(muAux(W), min=mumin, max=mumax)
  }

  this$.muAuxtab <- thresholdedMuAux
})


setMethodS3("initializeMu", "NPVI", function(this, muAux, g, ...) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'muAux':
  if (mode(muAux) != "function") {
    throw("Argument 'muAux' should be a function, not a ", mode(muAux));
  }
  ## Argument 'g':
  if (mode(g) != "function") {
    throw("Argument 'g' should be a function, not a ", mode(g));
  }

  setMuAux(this, muAux)

  #<Cabral>
  ## mu
  mu <- function(W) {
    muAux(W)* as.vector((1-g(W)))
  }
  #<\Cabral>
  setMu(this, mu)

  ## tabulated version of 'muAux'
  fW <- getFW(this);
  obs <- getObs(this);  
  MUAUXTAB <- muAux(fW(obs)); ## a *vector*, not a function
  
  muAuxtab <- function(ii) {
    MUAUXTAB[ii];
  }
  
  setMuAuxTab(this, muAuxtab)
  
  ## tabulated version of 'mu', as well as a tabulated version of pureMu
  MUTAB <- mu(fW(obs)); ## a *vector*, not a function
  
  mutab <- function(ii) {
      MUTAB[ii];
  }
  setMuTab(this, mutab)
})

setMethodS3("updateMu", "NPVI", function(this, dev, exact=TRUE, effICW, ...) {
  updateMuNonTab(this, dev, exact=exact, effICW, ...)
  updateMuTab(this, dev, exact=exact, effICW, ...)
})

#<Cabral>
setMethodS3("updateMuNonTab", "NPVI", function(this, dev, exact=TRUE, effICW, ...) {
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
  isBetaMult <- getIsBetaMult(this)
  if(isBetaMult){
    eps <- matrix( as.numeric(eps), ncol = 1) 
  }

  
  if(isBetaMult){
    puremu <- getPureMu(this, tabulate=FALSE);
  }
  else{
    puremu <- getMu(this, tabulate=FALSE);
  }
  
  g <- getG(this, tabulate=FALSE);
  
  if (!exact) { ## if do not use exact expression
    mu1 <- function(W) {
      if(isBetaMult){
        puremu(W) + t( eps %*% dev(W) ); 
      }else{
        puremu(W) + eps * dev(W);        
      }      
    }
  } else { ## if use exact expression
    mu1 <- function(W) {
      muW <- puremu(W);
      theEffICW <- effICW(W)
      
      if(isBetaMult){
        numerator <- muW +  ( devW + muW * theEffICW ) %*%  eps ;         
        denominator <- as.numeric( 1 + theEffICW  %*% eps );
      }
      else{      
        numerator <- muW + eps * (dev(W) + muW*theEffICW);
        denominator <- 1 + eps*theEffICW;      
      }
      
      numerator/denominator;
    }
  }
  
  
  muAux1 <- function(W) {
    mu1(W)/(1-g(W))
  }
  setMuAux(this, muAux1);
  setMu(this, mu1);
})


setMethodS3("updateMuTab", "NPVI", function(this, dev, exact=TRUE, effICW, ...) {
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
  isBetaMult <- getIsBetaMult(this)  
  eps <- getEpsilon(this)
  if(isBetaMult){
    eps <- matrix(as.numeric(eps), ncol=1)
  }

  if(isBetaMult){
    puremu <- getPureMu(this, tabulate=TRUE)
  }
  else{
    puremu <- getMu(this, tabulate=TRUE)
  }
  
  g <- getG(this, tabulate=TRUE)
  obs <- getObs(this, tabulate=TRUE);
  
  if(isBetaMult){
    muW <- puremu(obs[, "W"])
  }
  else{
    muW <- mu(obs[,"W"])
  }
  
  gW <- g(obs[, "W"])
  rm(puremu, g, obs)

  obs <- getObs(this)
  devW <- dev(fW(obs))
  W <- obs[, "W"]
  #rm(obs)
  
  if (!exact) { ## do not use exact expression
    if(isBetaMult){
      mu1W <- muW + t(eps %*%  devW ) ;
    }
    else{
      mu1W <- muW + eps * devW;
    }
  } else { ## use exact expression
    theEffICW <- effICW(W)
    
    if(isBetaMult){
      numerator <- muW + ( devW + muW * theEffICW) %*% eps;    
      denominator <- 1 +  theEffICW %*% eps  ;
    }else{     
      ## the above should use the tabulated or real versions of mu and theta0
      ## depending on tabulate, because effICW works on true values or
      ## indices depending on 'tabulate' (see how it is created in 'NPVI.update')
      numerator <- muW + eps * (devW + muW*theEffICW);
      denominator <- 1 + eps*theEffICW;
    }    
    mu1W <- numerator/denominator;
  }

  #browser()
  muAux1W <- mu1W/(1-gW)
  muAux1tab <- function(ii) {
    if(is.matrix(muAux1W)){
      muAux1W[ii,]
    }
    else{
      muAux1W[ii];
    }        
  }
  setMuAuxTab(this, muAux1tab)

  mu1tab <- function(ii) {
    if(is.matrix(mu1W)){
      mu1W[ii,]
    }
    else{
      mu1W[ii];
    }    
  }
  setMuTab(this, mu1tab)
})

#<\Cabral>

############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

