#<Cabral>
##
## it is important to not here that in the case where isBetaMult is true, muW, is not strictly speaking muW but rather a simplify version of it.
## To be more specific, we know that muW = E[ X h(W) | W]. In the case below, we called muW = E[X|W]
##
estimateDevMu <- function(muW, obs, eic1, flavor=c("learning", "superLearning"), learnDevMu, isBetaMult = FALSE, 
                          light=TRUE, SuperLearner.=NULL, ..., verbose=FALSE) {
#<\Cabral>
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'mu':
  muW <- Arguments$getNumerics(muW);  
  
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
  
  if(isBetaMult){
    ## devG will be function of the number of columns of the variable eic
    d <- ncol(eic1)
    devMu_fn <- list(NA,d)
    if(flavor != "learning"){
      fitDevMu <- list(NA,d)
    }
  }

  if (flavor=="learning") {
    #<Cabral>
    if(isBetaMult){
      for(i in 1:d){
        devMu_fn[[i]] <- learnDevMu(obs, eic1[,i], muW, light=light, ...);
      }
      
      devMu <- function(W){
        n <- nrow(W)
        out <- matrix(NA,n,d)
        for(i in 1:d){
          out[,i] <- devMu_fn[[i]](W)
        }
        out
      }
      
    }
    else{
      devMu <- learnDevMu(obs, eic1, muW, light=light, ...);
    }
    
    #browser()
    #verbose && cat(verbose, "devMu(W):");
    #verbose && print(verbose, summary(devMu(extractW(obs))));
    
  } else if (flavor=="superLearning") {
    logSL <- as.logical(less(verbose, 10));  ## decrease verbosity in superLearner
    SL.library.devMu <- learnDevMu
    obsD <- as.data.frame(obs)
    
    #<Cabral>
    
    if(isBetaMult){
      for(i in 1:d){
        ZdevMu <- (obsD[, "X"] - muW) * eic1[,i];
        fitDevMu[[i]] <- SuperLearner.(Y=ZdevMu, X=extractW(obsD), ## obsD[, "W", drop=FALSE]
                                  SL.library=SL.library.devMu, verbose=logSL,
                                  family=gaussian(), ...);
        devMu_fn[[i]] <- function(W,j) {
          Wd <- as.data.frame(W)
          predict(fitDevMu[[j]], newdata=Wd)$pred;
        }
      }
      
      devMu <- function(W){
        n <- nrow(W)
        out <- matrix(NA,n,d)
        for(i in 1:d){
          out[,i] <- devMu_fn[[i]](W,i)
        }
        out
      }
    }
    else{
      ZdevMu <- (obsD[, "X"] - muW) * eic1;
      fitDevMu <- SuperLearner.(Y=ZdevMu, X=extractW(obsD), ## obsD[, "W", drop=FALSE]
                                SL.library=SL.library.devMu, verbose=logSL,
                                family=gaussian(), ...);
      devMu <- function(W) {
        Wd <- as.data.frame(W)
        predict(fitDevMu, newdata=Wd)$pred;
      }
    }
    
    #verbose && cat(verbose, "devMu(W):");
    #verbose && print(verbose, summary(devMu(extractW(obs))));
  }
  
  
  verbose && cat(verbose, "devMu(W):");
  verbose && print(verbose, summary(devMu(extractW(obs))));
  
  #<\Cabral>
  devMu
}


############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

