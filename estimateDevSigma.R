##
## New file introduce for the estimation of Sigma
##

## Sigma = E[ X^2 | W]
## mu = E[X | W] - in the case of the parameter used by Pierre and Antoine
## given the above, we will keep the variable learnDevMu in the case the flavor is "learning" and simply insur that the values of the column
## corresponding to the values of X are squared before being used...

estimateDevSigma <- function(sigmaW, obs, eic1, flavor=c("learning", "superLearning"), learnDevMu,
                          light=TRUE, SuperLearner.=NULL, ..., verbose=FALSE) {

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'sigmaW':
  sigmaW <- Arguments$getNumerics(sigmaW);  
  
  ## Argument 'obs':
  obs <- validateArgumentObs(obs, allowIntegers=TRUE);
  
  ## Argument 'eic1'
  eic1 <- Arguments$getNumerics(eic1);
  
  ## Argument 'flavor':
  flavor <- match.arg(flavor);
  learnDevMode <- switch(flavor,
                      learning="function",
                      superLearning="character");

  ## Argument 'learnDevMu'
  mode <- mode(learnDevMu);
  if (mode != learnDevMode) {
    throw("Argument 'learnDevMu' should be of mode '", learnDevMode, "', not '", mode, "' for flavor: ", flavor);
  }

  ## Argument 'SuperLearner.'
  if (flavor=="superLearning") {
    if (is.null(SuperLearner.) || mode(SuperLearner.)!="function") {
      throw("Argument 'SuperLearner.' should be a function")
    }
  }

  ## Argument 'verbose'
  verbose <- Arguments$getVerbose(verbose);
  
  d <- ncol(eic1)
  devSigma_fn <- list(NA,d)
  if(flavor != "learning"){
    fitDevSigma <- list(NA,d)
  }

  if (flavor=="learning") {
    for( i in 1:d){
      X2 <- obs[,"X"]^2
      pseudoObs <- obs
      pseudoObs[,"X"] <- X2
      devSigma_fn[[i]] <- learnDevMu(pseudoObs, eic1[,i], sigmaW, light=light, ...);
    }
  } else if (flavor=="superLearning") {
    logSL <- as.logical(less(verbose, 10));  ## decrease verbosity in superLearner
    SL.library.devMu <- learnDevMu
    obsD <- as.data.frame(obs)
    
    for(i in 1:d){
      ZdevMu <- (obsD[, "X"]^2 - sigmaW) * eic1[,i];
    
    
      fitDevSigma[[i]] <- SuperLearner.(Y=ZdevMu, X=extractW(obsD), ## obsD[, "W", drop=FALSE]
                                SL.library=SL.library.devMu, verbose=logSL,
                                family=gaussian(), ...);
      
      devSigma_fn[[i]] <- function(W, j) {
        Wd <- as.data.frame(W)
        predict(fitDevSigma[[j]], newdata=Wd)$pred;
      }
    }
  }
  
  devSigma <- function(W){
    n <- nrow(W)
    out <- matrix(NA,n,d)
    for(i in 1:d){
      if(flavor == "learning"){
        out[,i] <- devSigma_fn[[i]](W)
      }
      else{
        out[,i] <- devSigma_fn[[i]](W,i)
      }  
    }
    out
  }   
  
  verbose && cat(verbose, "devSigma(W):");
  verbose && print(verbose, summary(devSigma(extractW(obs))));
  
  devSigma
}


############################################################################
## HISTORY:
## 2015-09-10
## o Created.
############################################################################

