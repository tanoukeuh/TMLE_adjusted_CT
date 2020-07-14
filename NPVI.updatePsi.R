setMethodS3("updatePsi", "NPVI", function(this, B, ..., verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ## Argument 'B':
  B <- Arguments$getInteger(B);
  
  ## Argument 'verbose'
  verbose <- Arguments$getVerbose(verbose);

  ## Retrieve parameters
  family <- getFamily(this)
  fX <- getFX(this)
  
  obs <- getObs(this)
  W <- obs[, "W"]
  
  X <- fX(obs)
  Xq <- getXq(this)
  
  g <- getG(this);  
  mu <- getMu(this);
  #<Cabral>
  derivF_beta <- NULL
  #<\Cabral>
  muAux <- getMuAux(this);
  sigma2 <- getSigma2(this);
  weightsW <- getWeightsW(this);
  
  #<Cabral>
  #obsW <- fW(obs)
  isBetaMult <- getIsBetaMult(this);
  derivH_beta <- NULL
  obsFull <- NULL
  if(isBetaMult){
    derivF_beta <- getDerivF_beta(this)
    derivH_beta <- getDerivH_beta(this);
    obsFull <- getObs(this, tabulate=FALSE)
  }
  #<\Cabral>

  ## Perform 'B' simulations according to the estimated parameters
  verbose && enter(verbose, "Simulating ", B, " observations");  
  
  #<Cabral>
  sigmaAux <- NULL
  if(isBetaMult){
    ## computing sigma2_base which will be use for the estimation below
    sigma <- getSigma(this)
    sigmaW <- sigma(W)
    sigma2 <- mean(sigmaW)
    
    ## retreiving the value of sigmaAux
    sigmaAux <- getSigmaAux(this)
    ## changing the value of mu since we only want E[X | W]
    mu <- getPureMu(this)
  }
  
  obsB <- simulateData(B, W, X, Xq, g, mu, muAux, sigma2, weightsW=weightsW, family=family, verbose=verbose, sigmaAux=sigmaAux)    
  if( sum(is.na(obsB[,"X"])) != 0 ){
    #browser()
  }
  
  if(isBetaMult){
    sigma2 <- getSigma2(this)
  }
  #<\Cabral>
  
  verbose && str(verbose, obsB);
  verbose && exit(verbose);

  ## Calculate 'theta' and 'theta0' on these B samples
  theta <- getTheta(this)
  theta0 <- getTheta0(this)  

  ## Estimate psiPn:
  #<Cabral>
  realObs <- extractCovariatesData(this, obs)  ## To replace with the function getObs, using tabulate = FALSE
  psi1 <- estimatePsi(theta=theta, theta0=theta0, fX=fX, obs=obs, sigma2=sigma2, realObs=realObs, derivF_beta = derivF_beta, derivH_beta=derivH_beta, isBetaMult=isBetaMult, verbose=verbose) 
  this$.psiPn <- psi1$mean;
  this$.psiPn.sd <- psi1$sd;  
  
  ## Estimate psi:
  realObs <- extractCovariatesData(this, obsB)  ## To replace with the function getObs, using tabulate = FALSE
  psi0 <- estimatePsi(theta=theta, theta0=theta0, fX=fX, obs=obsB, sigma2=sigma2, realObs=realObs, derivF_beta=derivF_beta, derivH_beta=derivH_beta, isBetaMult=isBetaMult, verbose=verbose) 
  this$.psi <- psi0$mean;
    
  # estimating psi.sd
  if(isBetaMult){
    eic <- this$.efficientInfluenceCurve
    eic <- eic$eic
    
    if(!is.null(eic)){      
      varCovariance <- cov(eic)
      varCovarianceSd <- varCovariance/nrow(eic)
      #this$.psi.sd.norm2 <- norm(varCovarianceSd,"F")
      this$.psi.sd.norm2 <- varCovariance
      this$.psi.sd <- sqrt(diag(varCovarianceSd))
      #browser()
    }
    else{
      this$.psi.sd <- rep(0, length(psi0$mean))
    }
    # updated fBeta_fit
    updateFbeta_fit(this, psi0$mean)    
  }
  else{
    this$.psi.sd <- psi0$sd;
  }
  #<\Cabral>  
  
  rm(obsB, realObs)
})

############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

