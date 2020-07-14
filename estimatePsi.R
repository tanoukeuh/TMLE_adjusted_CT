#<Cabral>
estimatePsi <- function(theta, theta0, fX, obs, sigma2, realObs=NULL, isBetaMult=FALSE, derivF_beta = NULL, derivH_beta = NULL, isSimulationScheme=FALSE, ..., verbose=FALSE){
#<\Cabral>
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  ## Argument 'theta':
  mode <- mode(theta);
  if (mode != "function") {
    throw("Argument 'theta' should be of mode 'function', not '", mode);
  }

  ## Argument 'theta0':
  mode <- mode(theta0);
  if (mode != "function") {
    throw("Argument 'theta0' should be of mode 'function', not '", mode);
  }
  
  ## Argument 'derivF_Beta'
  #<Cabral>
  if(isBetaMult){
    mode <- mode(derivF_beta);  
    if( mode != "function"){
      throw("Argument 'derivF_beta' should be of mode 'function', not ", mode);
    }
  }
  #<\Cabral>
  
  ## Argument 'fX':
  mode <- mode(fX);
  if (mode != "function") {
    throw("Argument 'fX' should be of mode 'function', not '", mode);
  }
  
  ## Argument 'obs':
  obs <- validateArgumentObs(obs);
  
  ## Argument 'sigma2':
  #<Cabral>
  sigma2 <- Arguments$getNumerics(sigma2);
  #<\Cabral>
  
  ##--------------------------------------------------------
  ##    Computation of the variable 'psi'
  ##--------------------------------------------------------
  
  
  if( !isSimulationScheme ){
    T <- theta( obs[,c("X","W")]);    
    #T <- theta(obs[, c("X", "W")]);
  }
  else{
    #XW <- cbind( X = obs[,"X"], extractW(obs))
    XW <- cbind( X = obs[,"X"], W = obs[,"W"])
    T <- theta(XW);
  }  
  
  verbose && cat(verbose, "theta(X, W):");
  verbose && str(verbose, T);
  
  if( !isSimulationScheme ){
    T0 <- theta0( obs[,"W", drop= FALSE]);
    #T0 <- theta0(obs[, "W", drop=FALSE]);
  }
  else{
    #T0 <- theta0( extractW(obs));
    T0 <- theta0( obs[,"W"]);
  } 
  
  verbose && cat(verbose, "theta0(W):");
  verbose && str(verbose, T0);
  
  if(isBetaMult){
    X <- realObs[,"X"];
    W <- extractW(realObs);
    
    cleanIndices <- is.na(X) 
    
    derivF_betaXW <- derivF_beta(X,W)    
    
    ##
    ## colnames(derivF_betaXW) <- AdjustColumnsNamesForPsi(derivF_betaXW , W)
    ##    
    
    ### creating a formula in order to run glm fit  
    originalNames <- colnames(derivF_betaXW)
    lengthColNames <- length(colnames(derivF_betaXW))
    
    pseudoColNames <- paste( "V" ,1:lengthColNames,  sep="" )
    colnames(derivF_betaXW) <- pseudoColNames
    ffPsi <- as.formula( paste("Y~", paste( 
                                   pseudoColNames, collapse = "+"
                                   ), "-1", sep="" ) )
    
    # ffPsi <- as.formula(paste( paste("Y~", 
    #                               paste(colnames(derivF_betaXW), collapse="+"),
    #                               sep="") , "-1", sep=""))
        
    ##ffPsi <- as.formula(paste( paste("Y~", 
    ##                                 paste(colnames(derivF_betaXW), collapse="+"),
    ##                                 sep="") ,  sep=""))    
    
    Y <- T-T0  
    #if(!isSimulationScheme){
    #  Y <- as.matrix(Y)
    #  Y <- Y[!cleanIndices,]
    #}
    data <- cbind( derivF_betaXW , Y = Y)      
    ##data <- cbind(X, W, Y=Y)
    
    psiFit <- glm(ffPsi, data = as.data.frame(data) , family=  gaussian)
    out <- list(mean = psiFit$coefficients, sd = summary(psiFit)$coefficients[,2]) 
    names(out$mean) <- originalNames
    names(out$sd) <- originalNames
    
    ## compute theoretical value of Psi(P) when formula 'ff' is 'formula(Y~I(X*V) + I(X*W))'
    ## betaW <- mean(Y*X*W[,"W"])/mean((X*W[,"W"])^2)    
    #derivH_betaW <- derivH_beta(W)        
    #invSigma2 <- solve(sigma2)    
    #psiTheory <- invSigma2 %*% apply( X* Y * derivH_betaW , 2  , mean)
    
    rm(psiFit, X , W , derivF_betaXW , Y, T, T0, ffPsi, data)    
    out;
  }
  else{
    mean.psi1 <- mean(fX(obs) * (T - T0))/sigma2;
    sd.psi1 <- sd(fX(obs) * (T - T0))/sigma2;
    list(mean=mean.psi1, sd=sd.psi1);
  }      
  
}


############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

