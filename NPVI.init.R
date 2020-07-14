setMethodS3("init", "NPVI", function(this, flavor=c("learning", "superLearning"),
                                     cvControl=NULL,
                                     learnG=NULL,
                                     learnMuAux=NULL,
                                     learnTheta=NULL,
                                     bound=1e-1, B=1e4,
                                     light=TRUE, 
                                     trueGMu=NULL,
                                     trueTheta=NULL,
                                     SuperLearner.=NULL,
                                     ..., verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ## Argument 'flavor':
  flavor <- match.arg(flavor);
  learnMode <- switch(flavor,
                      learning="function",
                      superLearning="character");
  
  ## Argument 'learnG'
  mode <- mode(learnG);
  if (mode != learnMode) {
    throw("Argument 'learnG' should be of mode '", learnMode, "', not '", mode, "' for flavor: ", flavor);
  }

  ## Argument 'learnMuAux'
  mode <- mode(learnMuAux);
  if (mode != learnMode) {
    throw("Argument 'learnMuAux' should be of mode '", learnMode, "', not '", mode, "' for flavor: ", flavor);
  }

  ## Argument 'learnTheta'
  mode <- mode(learnTheta);
  if (mode != learnMode) {
    throw("Argument 'learnTheta' should be of mode '", learnMode, "', not '", mode, "' for flavor: ", flavor);
  }

  ## Argument 'bound':
  bound <- Arguments$getNumeric(bound);
  if (bound<=0) {
    throw("Argument 'bound' must be positive!\n")
  }
  
  ## Argument 'B':
  B <- Arguments$getInteger(B);
  
  ## Argument 'light'
  light <- Arguments$getLogical(light);

  ## Argument 'trueGMu'
  useTrueGMu <- (!is.null(trueGMu))
  if (useTrueGMu) {
    if (!is.list(trueGMu)) {
      throw("If not NULL, Argument 'trueGMu' should be a list")
    }
    trueG <- trueGMu[["g"]]
    if (mode(trueG) != "function") {
      throw("Argument 'trueGMu$g' should be a function, not a ", mode(trueG))
    }
    trueMuAux <- trueGMu[["muAux"]]
    if (mode(trueMuAux) != "function") {
      throw("Argument 'trueGMu$muAux' should be a function, not a ", mode(trueMuAux))
    }
  }
  
  #<Cabral>  
  ## Argument 'trueTheta'
  useTrueTheta <- (!is.null(trueTheta))
  
  if(useTrueTheta){
    if(!is.list(trueTheta)){
      throw("If not NULL, Argument 'trueTheta' should be a list")
    }
    
    trueThetaVal <- trueTheta[["theta"]]
    if(mode(trueThetaVal) != 'function'){
      throw("Argument 'trueTheta$theta' should be a function, not a ", mode(trueThetaVal))
    }
  }  
  #<\Cabral>

  ## Argument 'SuperLearner.'
  if (flavor=="superLearning") {
    if (is.null(SuperLearner.) || mode(SuperLearner.)!="function") {
      throw("Argument 'SuperLearner.' should be a function")
    }
  }
  
  ## Argument 'verbose'
  verbose <- Arguments$getVerbose(verbose);
  verbose <- less(verbose, 10);

  ## retrieving 'obs'
  obs <- getObs(this, tabulate=FALSE);
  
  ##<Cabral>
  ##------------------------------------------------------------------------------------------------------------
  ##
  ##   Initializing the functions f_beta, derivF_beta and sigma2(which will be a matrix)
  ##
  ##------------------------------------------------------------------------------------------------------------
  isBetaMult <- getIsBetaMult(this)
  if(isBetaMult){
    ff <- getFormula(this)    
    internalFunctions <- create_fBeta_And_DerivFBeta(this, ff,obs)
    
    fBeta_fit <- internalFunctions[[1]]
    f_beta <- internalFunctions[[2]]
    derivF_beta <- internalFunctions[[3]]
    h_beta <- internalFunctions[[4]]
    derivH_beta <- internalFunctions[[5]]
    
    setFbeta_fit(this,fBeta_fit);
    setF_beta(this,f_beta);
    setDerivF_beta(this, derivF_beta);
    setH_beta(this, h_beta);
    setDerivH_beta(this, derivH_beta);
    
    ## computing the value of Sigma = E_P[derivF_beta(X,W) * t(derivF_beta(X,W))]    
    #pseudoObs <- identifiablePseudoObs(ff, obs)
    #X <- obs[,"X"]
    #W <- extractW(obs)
    #derivFBetaXW <- derivF_beta(X,W)
    #H.list <- lapply(1:nrow(derivFBetaXW), function(ii) {
    #  row <- matrix(as.numeric(derivFBetaXW[ii, ]), ncol=1)
    #  mat <- row %*% t(row)
    #})
    
    
    ### to do :: To modify the below...
    #Sigma2 <- Reduce("+", H.list)/nrow(obs)
    #rm(X,W,derivFBetaXW,H.list)  
    
    #setSigma2(this,Sigma2);
  }
  ##<\Cabral>
  
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## learning
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  verbose && enter(verbose, "Estimating relevant features of the distribution");
  
  if (!useTrueGMu) {
    g <- estimateG(obs, flavor=flavor, learnG=learnG, light=light,
                   SuperLearner.=SuperLearner.,
                   ..., verbose=verbose);
    muAux <- estimateMuAux(obs, flavor=flavor, learnMuAux=learnMuAux, light=light,
                           SuperLearner.=SuperLearner.,
                           ..., verbose=verbose);
  } else {
    g <- trueG
    muAux <- trueMuAux
  } 
  
  #<Cabral>
  sigmaAux <- estimateSigmaAux(obs, flavor=flavor, learnMuAux=learnMuAux, light=light,
                         SuperLearner.=SuperLearner.,
                         ..., verbose=verbose);
  #<\Cabral>
  
  initializeG(this, g);
  initializeMu(this, muAux, g); 
  
  #<Cabral>
  if(isBetaMult){
    initializeSigma(this, sigmaAux,g);
  }
  #<\Cabral>
  
  #<Cabral>
  if(!useTrueTheta){
  theta <- estimateTheta(obs, flavor=flavor, learnTheta=learnTheta, light=light,
                         SuperLearner.=SuperLearner.,
                         ..., verbose=verbose);
  }
  else{
    theta <- trueThetaVal;
  }
  #<\Cabral>
  
  initializeTheta(this, theta);
  
  ##
  ##
  ##  Initializing 'psi'
  ##
  ##

  #<Cabral>
  if(!isBetaMult){
    sigma2 <- mean(obs[, "X"]^2);
    setSigma2(this, sigma2);
    verbose && exit(verbose);
  }
  #<\Cabral>

  verbose && enter(verbose, "Updating 'psi' accordingly");
  updatePsi(this, B, verbose=verbose);
  psi0 <- getPsi(this);
  verbose && str(verbose, psi0);
  
  verbose && enter(verbose, "Updating efficient influence curve and 'epsilon' accordinlgy");
  updateEfficientInfluenceCurve(this);

  ## Update history
  updateHistory(this);

  verbose && exit(verbose);
})

############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

