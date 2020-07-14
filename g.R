setMethodS3("getGmin", "NPVI", function(this, ...) {
  this$.gmin;
})

setMethodS3("getGmax", "NPVI", function(this, ...) {
  this$.gmax;
})

setMethodS3("getG", "NPVI", function(this, tabulate, ...) {
  if (missing(tabulate)) {
    tabulate <- getTabulate(this);
  }
  if (!tabulate) {
    this$.g;
  } else {
    this$.gtab;
  }
})

setMethodS3("setG", "NPVI", function(this, g, ...) {
  ## Argument 'g':
  if ((!is.null(g))  && (mode(g)!="function")) {
    throw("Argument \var{g} should be of mode 'function', not ", mode(g));
  }

  gmin <- getGmin(this)
  gmax <- getGmax(this)
  thresholdedG <- function(W) {
    threshold(g(W), min=gmin, max=gmax)
  }
  this$.g <- thresholdedG
})

setMethodS3("setGTab", "NPVI", function(this, g, ...) {
  ## Argument 'g':
  if ((!is.null(g))  && (mode(g)!="function")) {
    throw("Argument \var{g} should be of mode 'function', not ", mode(g));
  }

  gmin <- getGmin(this)
  gmax <- getGmax(this)
  thresholdedG <- function(W) {
    threshold(g(W), min=gmin, max=gmax)
  }

  this$.gtab <- thresholdedG
})

setMethodS3("initializeG", "NPVI", function(this, g, ...) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'g':
  if (mode(g) != "function") {
    throw("Argument 'g' should be a function, not a ", mode(g));
  }

  ## g
  setG(this, g)

  ## tabulated version of 'g'
  fW <- getFW(this);
  obs <- getObs(this);
  GTAB <- g(fW(obs)); ## a *vector*, not a function
  gtab <- function(ii) {
    GTAB[ii];
  }
  setGTab(this, gtab)
})

setMethodS3("updateG", "NPVI", function(this, dev, exact=TRUE, effICW, ...) {
  updateGNonTab(this, dev, exact=exact, effICW, ...)
  updateGTab(this, dev, exact=exact, effICW, ...)
})

#<Cabral>
setMethodS3("updateGNonTab", "NPVI", function(this, dev, exact=TRUE, effICW, ...) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'dev':
  
  #<Cabral>
  mode <- mode(dev)
  if (mode != "function") {
    throw("Argument 'dev' should be a function, not a ", mode);
  }
  #<\Cabral>

  ## Argument 'exact':
  exact <- Arguments$getLogical(exact);
  
  ## Argument 'effICW':
  if (exact) {
    if (missing(effICW)) {
      throw("Argument 'effICW' is required when 'exact' is TRUE");
    }
    if (mode(effICW) != "function") {
      throw("Argument 'effICW' should be a function, not a ", mode(dev));
    }
  }

  fW <- getFW(this)
  eps <- getEpsilon(this)
  tabulate <- getTabulate(this)
  isBetaMult <- getIsBetaMult(this)

  g <- getG(this, tabulate=FALSE)

  if (!exact) { ## do not use exact expression
    g1 <- function(W) {
      logit <- qlogis
      expit <- plogis
      
      gW <- g(W);      
      if(isBetaMult){
        res <- logit(gW) + eps %*% dev(W) * 1/(gW*(1-gW));        
      }
      else{
        res <- logit(gW) + eps * dev(W) * 1/(gW*(1-gW));
      }      
      expit(res);
    }
  } else { ## use exact expression
    g1 <- function(W) {
      gW <- g(W)
      theEffICW <- effICW(W)
      if(isBetaMult){
        numerator <- gW +  (dev(W) + gW*theEffICW) %*% eps ;
        denominator <- 1 + theEffICW %*%  eps;        
      }else{      
        numerator <- gW + eps * (dev(W) + gW*theEffICW);
        denominator <- 1 + eps*theEffICW;
      }
      out <- numerator/denominator;
      return(out)
    }
  }
  setG(this, g1);
})

setMethodS3("updateGTab", "NPVI", function(this, dev, exact=TRUE, effICW, ...) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'dev':
  
  #<Cabral>
  mode <- mode(dev)
  if (mode != "function") {
    throw("Argument 'dev' should be a function, not a ", mode);
  }
  #<\Cabral>

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
  tabulate <- getTabulate(this)
  isBetaMult <- getIsBetaMult(this)

  g <- getG(this, tabulate=TRUE)
  obs <- getObs(this, tabulate=TRUE);
  gW <- g(obs[, "W"])
  rm(g, obs)

  obs <- getObs(this)
  devW <- dev(fW(obs))
  W <- obs[, "W"]
  rm(obs)
    
  if (!exact) { ## do not use exact expression
    logit <- qlogis
    expit <- plogis      
    if(isBetaMult){
      res <- logit(gW) +  ( devW %*% eps ) * 1/(gW*(1-gW));
    }else{    
      res <- logit(gW) + eps * devW * 1/(gW*(1-gW));
    }
    g1W <- expit(res);
  } else { ## use exact expression
    theEffICW <- effICW(W)
    ## the above should use the tabulated or real versions of mu and theta0
    ## depending on tabulate, because effICW works on true values or
    ## indices depending on 'tabulate' (see how it is created in 'NPVI.update')
    
    if(isBetaMult){
      numerator <- gW +  (devW + gW*theEffICW) %*% eps  ;
      denominator <- 1 +  theEffICW %*% eps ;
    }else{   
      numerator <- gW + eps * (devW + gW*theEffICW);
      denominator <- 1 + eps*theEffICW;
    }
    
    g1W <- numerator/denominator;
  }

  g1tab <- function(ii) {
    g1W[ii]
  }
  setGTab(this, g1tab)
})

#<\Cabral>

############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

