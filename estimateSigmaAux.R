##
## New file introduce for the estimation of Sigma
##

## Sigma = E[ X^2 | W]
## mu = E[X | W] - in the case of the parameter used by Pierre and Antoine
## given the above, we will keep the variable learnMuAux in the case the flavor is "learning" and simply insure that the values of the column
## corresponding to the values of X are squared before being used...

estimateSigmaAux <- function(obs, flavor=c("learning", "superLearning"), learnMuAux, isBetaMult = FALSE, 
                       light=TRUE, SuperLearner.=NULL, ..., verbose=FALSE) {

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'obs':
  obs <- validateArgumentObs(obs);
  
  ## Argument 'flavor':
  flavor <- match.arg(flavor);
  learnMode <- switch(flavor,
                      learning="function",
                      superLearning="character");

  ## Argument 'learnMuAux'
  mode <- mode(learnMuAux);
  if (mode != learnMode) {
    throw("Argument 'learnMuAux' should be of mode '", learnMode, "', not '", mode, "' for flavor: ", flavor);
  }

  ## Argument 'SuperLearner.'
  if (flavor=="superLearning") {
    if (is.null(SuperLearner.) || mode(SuperLearner.)!="function") {
      throw("Argument 'SuperLearner.' should be a function")
    }
  }
  
  ## Argument 'verbose'
  verbose <- Arguments$getVerbose(verbose);

  idx <- which(obs[, "X"] != 0);
  
  if (flavor=="learning") {
    obs[,"X"] <- obs[,"X"]^2
    sigmaAux <- learnMuAux(obs[idx, ], light=light, ...);
  } else if (flavor=="superLearning") {
    logSL <- as.logical(less(verbose, 10));  ## decrease verbosity in SuperLearner
    SL.library.muAux <- learnMuAux;
    obsD <- as.data.frame(obs)
    ## WW <- obsD[idx, "W", drop=FALSE]
    WW <- extractW(obsD[idx, ])

    Y <- obsD[idx,"X"]^2
    fitSigmaAux <- SuperLearner.(Y=Y, X=WW,
                              SL.library=SL.library.muAux, verbose=logSL,
                              family=gaussian(), ...);
    verbose && print(verbose, fitSigmaAux);
    sigmaAux <- function(W) {
      Wd <- as.data.frame(W)
      predict(fitSigmaAux, newdata=Wd)$pred;
    }
  }    
  
  verbose && cat(verbose, "sigma'(W):");
  verbose && print(verbose, summary(muAux(extractW(obs))));  
  
  sigmaAux
}

############################################################################
## HISTORY:
## 2015-09-10
## o Created.
############################################################################

