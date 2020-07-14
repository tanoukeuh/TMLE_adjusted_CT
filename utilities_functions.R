##
## reformulation of the function learningLib$learnDevMu
##
learnDevMu2<- function (obs, derivFBeta, effIC1, muW, light = TRUE, verbose = FALSE, ...){
  
  W <- extractW(obs)  
  Z <- (derivFBeta - muW) 
  
  Z.list <- lapply(1:nrow(derivFBeta), function(ii) {
    rowZ <- matrix(as.numeric(Z[ii, ]), ncol=1)
    rowEic <- matrix(as.numeric(effIC1[ii, ]), ncol=1)
    mat <- t(rowZ) %*% rowEic     
  })  
  Z<- unlist(Z.list)
  
  obsZ <- cbind(obs, Z = Z)
  verbose && str(verbose, obsZ)
  varNames <- colnames(W)
  theFormula <- paste(varNames, collapse = " + ")
  theFormula2 <- paste("I(", varNames, "^2)", collapse = " + ", 
                       sep = "")
  theFormula <- paste("Z ~", theFormula, "+", theFormula2, 
                      sep = "")
  formula <- as.formula(theFormula)
  fit <- glm(formula, data = as.data.frame(obsZ), family = gaussian)
  rm(W, Z, obsZ)
  if (light) {
    fit <- getLightFit(fit)
  }
  foo <- function(W) {
    predict(fit, newdata = data.frame(W), type = "response")
  }
  attr(foo, "fit") <- fit
  return(foo)
}



## Returns an identifiable version of the formula ff to be used for the
## computation of the parameter of interest

identifiablePseudoObs <- function(ff, ## formula used to define the function f_beta
                                  obs ## observations 
                                  ){
          pseudoObs <- model.frame(ff , as.data.frame(obs))
          fit <- glm(ff, data= as.data.frame(obs), family= gaussian)
          idx.identif <- is.null(fit$coef)          
          
          if( !is.na(match("(Intercept)", names(fit$coef)))){
            pseudoObs <-cbind(intercept=1,
                              pseudoObs[, !idx.identif & colnames(pseudoObs)!="Y",drop=FALSE])
          }
          else{
            pseudoObs <- pseudoObs[, !idx.identif & colnames(pseudoObs)!="Y"]
          }          
          
          rm(fit);
          #pseudoObs <- obs[,"X"]*pseudoObs   # There is no need to multiply by obs[,"X]
          pseudoObs
  }


# returns a vector containing two functions and a glm result.
# the below code has been encapsulated so that he can be re-used in both the initialisation of the 
# the parameter TMLE_NPVI but also the function 'getSample'

create_fBeta_And_DerivFBeta <- function( this=NULL,
                                         ff , # formula used to define the function f_beta
                                         obs # initial data used for the computation
){  
  pseudoObs <- model.frame(ff, as.data.frame(obs))      
  fBeta_fit <- glm(ff, data = as.data.frame(obs), family = gaussian)  
  
  ##
  ##   f_beta(X,W) = X * h_beta(W) 
  ##  
  
  h_beta <- function( W , beta = NULL ){
    
    #varNames <- colnames(W)    
    #test <- match(  c("W") , varNames )
    #idxs <- which(is.na(test));
    #if(length(idxs)){
    #  throw("Missing column:", varNames[idxs]);
    #}
    
    YW <- cbind(Y=1, W)
    pseudoObs <- identifiablePseudoObs(ff, YW)
    
    pseudoBetaFit <- NULL
    if(!is.null(this)){
      pseudoBetaFit <- this$.fBeta_fit;
    }
    
    ## providing defaults values in case the above has not been initialize yet
    if(is.null(pseudoBetaFit)){
      pseudoBetaFit <- fBeta_fit
    }
    
    ## replacing beta with its approbriate values if it has been provided..
    if( !(missing(beta) || is.null(beta)) ){
      dimBeta <- pseudoBetaFit$coefficients
      if( length(dimBeta) != length(beta)  ){        
        throw("The parameter beta provided doesn't have the proper length.")
      }      
      pseudoBetaFit$coefficients <- beta
    }    
    
    #result <- predict(pseudoBetaFit , newdata = data.frame(W), type = "response")    
    #result <- predict(pseudoBetaFit , newdata = data.frame(pseudoObs), type = "response")
    result <- predict(pseudoBetaFit , newdata = data.frame(YW), type = "response")
    result
  }
  
  f_beta <- function(X ,W , beta = NULL ){    
    h_betaW <- h_beta(W, beta)
    result <- X * h_betaW
  }
  
  ## checking that the funciton f_beta verifies : f_beta(0,W) = 0 
  n = length(fBeta_fit$coefficients)
  X <- rep(0, 4)
  W <- extractW(obs)[1:4,]  
  if( sum( f_beta(X,W) != 0) != 0 ){
    throw("The formula used is incorrect as it will introduce a function f_beta such that f(X,W) != 0 ")
  }
  
  derivF_beta <- function(X,W){
    XW <- cbind(X,W)
    
    ## adding a dummy column Y which will be removed later on 
    XW <- cbind( XW , Y=apply(XW,1, mean))
    
    #varNames <- colnames(XW)
    
    #test <- match(c("X","W") , varNames )
    #idxs <- which(is.na(test));
    #if(length(idxs)){
    #  throw("Missing column:", varNames[idxs]);
    #}    
    
    obsD <- X * identifiablePseudoObs(ff, XW)      
    obsD1 <- as.matrix(obsD)      
    obsD1
  }
  
  derivH_beta <- function(W){
    W <- cbind(W)    
    ## adding a dummy column Y which will be removed later on 
    W <- cbind( W , Y=apply(W,1, mean))
    
    #varNames <- colnames(W)
    
    #test <- match(c("W") , varNames )
    #idxs <- which(is.na(test));
    #if(length(idxs)){
    #  throw("Missing column:", varNames[idxs]);
    #}    
    
    obsD <- identifiablePseudoObs(ff, W)      
    
    obsD1 <- as.matrix(obsD)      
    obsD1
  }
  
  out <- list(NA,4)  
  out[[1]] <- fBeta_fit
  out[[2]] <- f_beta
  out[[3]] <- derivF_beta
  out[[4]] <- h_beta
  out[[5]] <- derivH_beta
  out
}


###
### determine the proper association names between the formula names coming from derivF_beta and the extract names of W
###


AdjustColumnsNamesForPsi <- function(derivF_beta , W){
  colnamesDerivF <- colnames(derivF_beta)
  colnamesW <- colnames(W)
  
  n <- length(colnamesW)
  d <- length(colnamesDerivF)
  result <- rep(NA,d)
  
  for( i in 1:d){
    for( j in 1:n){
      if( grepl(colnamesW[j], colnamesDerivF[i]) == TRUE ){
        result[i] <- colnamesW[j]
      }
    }
  }
  
  result
}


###
###
### helper function to graph evolution of learning and superLearning method..
### It also compares results with the simulation scheme
###
###

TMLE_PlotHelp <- function( resultLearning , resultSuperLearning, resultSimulation , folderName){
  
  class <- class(resultLearning)[1]
  if( class != "NPVI_CA"){
    throw("The variable resultLearning should be a list not a ", class)
  }
  
  class <- class(resultSuperLearning)[1]
  if( class != "NPVI_CA"){
    throw("The variable resultSuperLearning should be a list not a ", class)
  }
  
  mode <- mode(resultSimulation)
  if( mode != "list"){
    throw("The variable resultSimulation should be a list not a ", mode)
  }
  
  colNames <- names(resultLearning$history$psi[[1]])
  
  ##
  ## creating folder name
  ##
  
  #mainPathName <- "/Desktop/phd_theis/covariate-adjusted variable with parametric function/graphs/"
  #pathName <- paste( mainPathName , folderName,"/", sep="")
  
  ## 
  ## Plotting the result of the learning method
  ##
  
  #pathFileName <- paste(pathName , "LearningMethodIterations.jpeg", sep="")
  pathFileName <- "LearningMethodIterations.jpeg"
  jpeg(file=pathFileName)
  UnlistedPsi <- unlist(resultLearning$history$psi)  
  matPsiLearning <- t( matrix( UnlistedPsi , nrow = length(colNames) , ncol=length(resultLearning$history$psi)  ) )  
  matplot(matPsiLearning , type="o", lty = 1:length(colNames) , pch = 1:length(colNames) , col=c(1:length(colNames)) , ylab="Values" , xlab = "Iterations" )
  legend("topleft", legend=colNames, lty = 1:length(colNames) , pch = 1:length(colNames) , cex = 0.8 , col=c(1:length(colNames)) )
  title(main="Learning Method")
  dev.off()
  
  
  ## 
  ## Plotting the result of the super Learning method
  ##
  
  #pathFileName <- paste(pathName , "LearningMethodIterations.jpeg")
  pathFileName <- "SuperLearningMethodIterations.jpeg"
  jpeg(file=pathFileName)
  UnlistedPsi <- unlist(resultSuperLearning$history$psi)  
  matPsiSuperLearning <- t( matrix( UnlistedPsi , nrow = length(colNames) , ncol=length(resultSuperLearning$history$psi)  ) )  
  matplot(matPsiSuperLearning , type="o", lty = 1:length(colNames) , pch = 1:length(colNames) , col=c(1:length(colNames)) , ylab="Values" , xlab = "Iterations" )
  legend("topleft", legend=colNames, lty = 1:length(colNames) , pch = 1:length(colNames) , cex = 0.8 , col=c(1:length(colNames)) )
  title(main="Super Learning Method")
  dev.off()
  
  
  ##
  ##
  ## comparing superLearner, learner and simulation results
  ##
  ##
  UnlistedPsiSd <- unlist(resultLearning$history$psi.sd)
  matPsiSdLearning <- t( matrix( UnlistedPsiSd , nrow = length(colNames) , ncol=length(resultLearning$history$psi)))
  
  UnlistedPsiSd <- unlist(resultSuperLearning$history$psi.sd)
  matPsiSdSuperLearning <- t( matrix( UnlistedPsiSd , nrow = length(colNames) , ncol=length(resultSuperLearning$history$psi)))
  
  simulationPsi <- resultSimulation$psi
  simulationSd <- resultSimulation$sd
  
  ## 
  ## Plotting the result of the learning method with error bars
  ##  
  
  pathFileName <- "LearningMethodIterations_withErrorBar.jpeg"
  jpeg(file=pathFileName)
  UnlistedPsi <- unlist(resultLearning$history$psi)  
  matPsiLearning <- t( matrix( UnlistedPsi , nrow = length(colNames) , ncol=length(resultLearning$history$psi)  ) ) 
  ylim <- range(matPsiLearning+matPsiSdLearning, 
                matPsiLearning-matPsiSdLearning,
                simulationPsi)
  matplot(matPsiLearning , type="o", lty = 1:length(colNames) , pch = 1:length(colNames) , col=c(1:length(colNames)) , ylab="Values" , xlab = "Iterations",
          xlim=c(1,nrow(matPsiLearning)), ylim=ylim)  
  plotCI(x=rep(1:nrow(matPsiLearning), ncol(matPsiLearning)), y=as.vector(matPsiLearning), uiw = as.vector(matPsiSdLearning) ,add=T)
  legend("topleft", legend=colNames, lty = 1:length(colNames) , pch = 1:length(colNames) , cex = 0.8 , col=c(1:length(colNames)) )
  title(main="Learning Method with Error Bars")
  dev.off()
  
  ## 
  ## Plotting the result of the super learning method with error bars
  ##
  
  pathFileName <- "SuperLearningMethodIterations_WithErrorBar.jpeg"
  jpeg(file=pathFileName)
  UnlistedPsi <- unlist(resultSuperLearning$history$psi)  
  matPsiSuperLearning <- t( matrix( UnlistedPsi , nrow = length(colNames) , ncol=length(resultSuperLearning$history$psi)  ) )  
  ylim <- range(matPsiSuperLearning+matPsiSdSuperLearning, 
                matPsiSuperLearning-matPsiSdSuperLearning,
                simulationPsi)
  matplot(matPsiSuperLearning , type="o", lty = 1:length(colNames) , pch = 1:length(colNames) , col=c(1:length(colNames)) , ylab="Values" , xlab = "Iterations" ,
          xlim=c(1,nrow(matPsiLearning)), ylim=ylim)  
  plotCI(x=rep(1:nrow(matPsiSuperLearning), ncol(matPsiSuperLearning)), y=as.vector(matPsiSuperLearning), uiw = as.vector(matPsiSdSuperLearning) ,add=T)
  legend("topleft", legend=colNames, lty = 1:length(colNames) , pch = 1:length(colNames) , cex = 0.8 , col=c(1:length(colNames)) )
  title(main="Super Learning Method with Error Bars")
  dev.off()
  
  
  for( i in 1:length(simulationPsi)){    
    ##
    ##  creating error bar for learning method
    ##
    
    fileName = paste(colNames[i] , "_Learning_VS_SIM.jpeg", sep="")
    jpeg(file=fileName)
    n = nrow(matPsiLearning) 
    ylim <- range(matPsiLearning+matPsiSdLearning, 
                  matPsiLearning-matPsiSdLearning,
                  simulationPsi)
    errbar(1:n, matPsiLearning[,i], matPsiLearning[,i] -  matPsiSdLearning[,i]  , matPsiLearning[,i] + matPsiSdLearning[,i] , type="o", lty = 1 , pch=1, col=1, ylab="Values" , xlab = "iterations", asp=2,
           ylim=ylim)    
    abline(a = simulationPsi[i] , b = 0 ,  col = 2,  lty = 2 , pch=2 )       
    legend("topleft", legend=c("Learner","Sim"), lty = 1:2 , pch = 1:2 , cex = 0.8 , col=c(1:2) )
    title(main= paste( colNames[i] , "- Learning Method", sep="")   ) 
    dev.off()
    
    ##
    ##  creating error bar for super learning method
    ##
    
    fileName = paste(colNames[i] , "_SuperLearning_VS_SIM.jpeg", sep="")
    jpeg(file=fileName)
    n = nrow(matPsiSuperLearning) 
    ylim <- range(matPsiSuperLearning+matPsiSdSuperLearning, 
                  matPsiSuperLearning-matPsiSdSuperLearning,
                  simulationPsi)
    errbar(1:n, matPsiSuperLearning[,i], matPsiSuperLearning[,i] -  matPsiSdSuperLearning[,i]  , matPsiSuperLearning[,i] + matPsiSdSuperLearning[,i] ,type="o", lty = 1 , pch=1, col=1, ylab="Values" , xlab = "iterations", asp=1,
           ylim=ylim)
    abline(a = simulationPsi[i] , b = 0 ,  col = 2,  lty = 2 , pch=2 )      
    legend("topleft", legend=c("SuperLearner","Sim"), lty = 1:2 , pch = 1:2 , cex = 0.8 , col=c(1:2) )
    title(main= paste( colNames[i] , " - Super Learning Method ", sep="")    )
    dev.off()
  }
  
}


###
###
### helper function to graph evolution of learning and superLearning method..
### It also compares results with the simulation scheme
###
###

TMLE_Result_Plotter <- function( result, resultSimulation=NULL , folderName=NULL){
  
  class <- class(result)[1]
  if( class != "NPVI"){
    throw("The variable result should be a NPVI and not a ", class)
  } 
  
  
  if(!is.null(resultSimulation)){
    mode <- mode(resultSimulation)
    if( mode != "list"){
      throw("The variable resultSimulation should be a list not a ", mode)
    }
  }
  
  resultType <- result$.flavor  
  colNames <- names(result$history$psi[[1]])  
  
  ## 
  ## Plotting the result of the method
  ##
  
  pathFileName <- paste( class , " - ", resultType , " - ", "MethodIterations.pdf", sep="")
  pdf(file=pathFileName)
  
  ## number of columns for the graphs
  nbColGraph <- 1
  if(!is.null(resultSimulation)){
    nbColGraph <- 2
  }
  par(mfrow=c(2,nbColGraph))
  
  UnlistedPsi <- unlist(result$history$psi)  
  matPsiLearning <- t( matrix( UnlistedPsi , nrow = length(colNames) , ncol=length(result$history$psi)  ) )  
  matplot(matPsiLearning , type="o", lty = 1:length(colNames) , pch = 1:length(colNames) , col=c(1:length(colNames)) , ylab="Values" , xlab = "Iterations" )
  legend("topleft", legend=colNames, lty = 1:length(colNames) , pch = 1:length(colNames) , cex = 0.8 , col=c(1:length(colNames)) )
  title(main= paste( resultType , " Method" , sep=""))
  
  #dev.off()  
  
  
  ##
  ## comparing superLearner, learner and simulation results
  ##  
  
  UnlistedPsiSd <- unlist(result$history$psi.sd)
  matPsiSd <- t( matrix( UnlistedPsiSd , nrow = length(colNames) , ncol=length(result$history$psi)))  
  
  if( ! is.null(resultSimulation)){
    simulationPsi <- resultSimulation$psi
    simulationSd <- resultSimulation$sd
  }
  
  ## 
  ## Plotting the result of the learning method with error bars
  ##  
  
  #pathFileName <- paste(  resultType ,"MethodIterations_withErrorBar.jpeg", sep="")
  #jpeg(file=pathFileName)
  
  UnlistedPsi <- unlist(result$history$psi)  
  matPsi <- t( matrix( UnlistedPsi , nrow = length(colNames) , ncol=length(result$history$psi)  ) ) 
  ylim <- range(matPsi + matPsiSd, 
                matPsi - matPsiSd,
                simulationPsi)
  matplot(matPsi , type="o", lty = 1:length(colNames) , pch = 1:length(colNames) , col=c(1:length(colNames)) , ylab="Values" , xlab = "Iterations",
          xlim=c(1,nrow(matPsi)), ylim=ylim)  
  plotCI(x=rep(1:nrow(matPsi), ncol(matPsi)), y=as.vector(matPsi), uiw = as.vector(matPsiSd) ,add=T)
  legend("topleft", legend=colNames, lty = 1:length(colNames) , pch = 1:length(colNames) , cex = 0.8 , col=c(1:length(colNames)) )
  
  title(main=paste( resultType, " Method with Error Bars", sep=""))
  
  
  #dev.off()
  
  
  if(!is.null(resultSimulation)){
    
    for( i in 1:length(simulationPsi)){    
      ##
      ##  creating error bar 
      ##
      
      #fileName = paste(colNames[i] , "_" , resultType, "_VS_SIM.jpeg", sep="")
      #jpeg(file=fileName)
      
      
      n = nrow(matPsi) 
      ylim <- range(matPsi+matPsiSd, 
                    matPsi-matPsiSd,
                    simulationPsi)
      errbar(1:n, matPsi[,i], matPsi[,i] -  matPsiSd[,i]  , matPsi[,i] + matPsiSd[,i] , type="o", lty = 1 , pch=1, col=1, ylab="Values" , xlab = "iterations", asp=2,
             ylim=ylim)    
      abline(a = simulationPsi[i] , b = 0 ,  col = 2,  lty = 2 , pch=2 )       
      legend("topleft", legend=c("Learner","Sim"), lty = 1:2 , pch = 1:2 , cex = 0.8 , col=c(1:2) )
      title(main= paste( colNames[i] , " - ", resultType, " Method", sep="")   ) 
      
      
      #dev.off()
      
    }
  }
  
  dev.off()
  
}


###
###    
###  Computing the plot related the different value 
###
###
plot_statisticalAnalysis_NPVICA <- function( result , simResult ,maxIter = NULL , partTitle = ""){  
  
  validINF <-  sapply(result, function(ll){!is.character(ll)}) 
  result <- result[validINF]
  
  mode <- mode(result)  
  if( mode != "list"){
    throw("The variable result should be a list not a ", mode)
  }
  
  resultType <- result[[1]]$.flavor
  
  nObs <- nrow(result[[1]]$.obs)
  nbIterations <- length(result)
  
  ##
  ## determinig the value of the max number of iterations in case it has not been fixed in the formula
  ##
  if(missing(maxIter)){
    maxIter = length(result[[1]]$history$psi)
  }
  
  
  ## The number of variables in our list
  nbElmts = length(result)
  
  ## nb of covariates in the analysis
  nbCovariates <- length(result[[1]]$history$psi[[1]])
  namesCovariates <- names(result[[1]]$history$psi[[1]])
  
  ##
  ##  storing the result of our analysis in a list where each element correspond to each of our iteration
  ##
  
  out <- replicate( nbCovariates, list(0) )
  for( i in 1:nbElmts){
    iResult <- result[[i]]
    matPsi <- t(matrix( unlist( iResult$history$psi), nrow=length(iResult$history$psi[[1]]), ncol=length(iResult$history$psi)) )
    
    for( j in 1 : length(iResult$history$psi) ){
      for( k in 1: nbCovariates){
        if(   i <= 1   ){          
          out[[k]] [[j]] <- list(0)          
          out[[k]] [[j]] <- matPsi[j,k]           
        }
        else{
          if( length(out[[k]]) >= j ){
            out[[k]] [[j]]   <- c(  out[[k]][[j]]  , matPsi[j,k] )    
          }
          else{
            out[[k]] [[j]] <- list(0)          
            out[[k]] [[j]] <- matPsi[j,k]
          }
        }
      }
    }    
  }  
  
  
  ##
  ##  creating a series of density points which will be needed for the graphs
  ##
  outDensity <- replicate( nbCovariates, list(0) )
  ylimDensity <- replicate( nbCovariates, list(0) )
  xlimDensity <- replicate( nbCovariates, list(0) )
  
  for(i in 1:nbCovariates){
    nbIter <- length(out[[i]])    
    
    for( j in 1:nbIter ){
      d <- density(out[[i]][[j]])
      outDensity[[i]][[j]] <- list(0)
      outDensity[[i]][[j]] <- d
      if( j==1){
        ylimDensity[[i]] <- range(d$y)
        xlimDensity[[i]] <- range(d$x)
      }else{
        ylimDensity[[i]] <- range( ylimDensity[[i]], d$y)
        xlimDensity[[i]] <- range( xlimDensity[[i]], d$x)
      }      
    }
  }
  
  ##
  ## plotting the graph showing the evolution of the empirical density at each iteration points
  ##  
  
  ## p represents the number of itereations steps that we would like to remove from the graph
  p <- 1
  
  ## m represents the starting iterations of our graph
  pStart <- 1
  
  for(i in 1:nbCovariates){
    fileName = paste(namesCovariates[i], " - Evolution Empirical Density - ", partTitle , " - ", resultType , ".pdf ", sep="")
    pdf(file=fileName)
    
    nbIter <- length(out[[i]])
    #d <- density(out[[i]][[1]])    
    d <- outDensity[[i]][[2]]
    plot( d  ,  type="o", lty = 1 , pch = 1 , col=1 , main=NA, xlab=NA, ylab=NA , xlim= xlimDensity[[i]], ylim=ylimDensity[[i]] )
    
    presence <- vector(length=nbIter-p)
    
    if(nbIter > 1 ){      
      for( j in pStart:(nbIter-p)){
        #d <- density(out[[i]][[j]])
        d <- outDensity[[i]][[j]]        
        presence[j-pStart+1] <- length(out[[i]][[j]])/nbIterations                
        lines( density(out[[i]][[j]])  ,  type="o", lty = j , pch = j , col=j )
      }
    }
    
    ## adding the result given by the simulation
    abline(v=simResult$psi[i], col=nbIter-p+1)    
    percent <- paste0( format( 100* presence, format="f", digits=3), "%" )    
    legend("topleft", legend=paste("k = ", (pStart):(nbIter-p) , "(" , percent , ")") , pch=(pStart):(nbIter-p), cex=0.8,  col=c((pStart):(nbIter-p)))
    title(main = paste( namesCovariates[i] , "- Evolution Density with ", partTitle , " ( N = ", nObs , " )" , " - ", "\n", nbIterations ," Iter. of " , resultType , sep=""))
    dev.off()
  }   
}




statisticalAnalysisProcedure_TMLE_NPVICA <- function(
  d,                ## number of cavariates to be considered
  N = 4,            ## number of procedure to be run      
  b = 300,          ## number of observations to be drawn in order to estimate our parameter of interest using TMLE_NPVI_CA 
  ff = NULL,        ## formula to be used for the formula      
  B = 5e4,         ## number of observations to be created in order to compute the value of the parameter of interest using the simulation scheme
  flavor="learning",## the learning procedure to use for some key variables of our procedures...
  maxIter = 5,      ## maximum number of iterations
  trueGMu = NULL,   ## true function for the varialbes "G" and "Mu"
  trueTheta=NULL,   ## true function for the varaibles "Theta"     
  dontComputeSim=FALSE  ## do not trigger the computation of the simulation values
) {
  
  ##
  ## compute the value of our parameter throught simulation scheme
  ##     
  
  if(!dontComputeSim){
    print("Computing parameter of interest value through simulation Scheme...")
    simulation <- getSample( B , O, lambda0,ff=ff, isBetaMult=TRUE, d=d, covariatesMean=rep(0,d), covariatesSigma=diag(d) , p=p, omega=omega, sigma2=1, Sigma3=S)  
  }
  
  ##
  ##  computing parameter of interest using TMLE_NPVI_CA
  ##  
  
  resultEstimate <- vector("list", N ) ;## initializing the variables which will old the entire set of result using "flavor" method    
  
  print("starting to compute estimate of parameter of interest using TMLE method...")
  for(i in 1:N){
    ##
    ## drawing the observation needed for our estimate scheme
    ##
    
    print( paste("Statistical Analysis - Step : ", i , sep="") )    
    sim1 <- getSample(b, O, lambda0,ff=ff, isBetaMult=TRUE, d=d, covariatesMean=rep(0,d), covariatesSigma=diag(d) , p=p, omega=omega, sigma2=1, Sigma3=S)
    obs <<- sim1$obs    
    
    tmle <- try(tmle.npvi(sim1$obs, f=identity, ff=ff, isBetaMult=TRUE, flavor=flavor, nMax=40, iter=maxIter, trueGMu=trueGMu, trueTheta=trueTheta) )
    failed <- inherits(tmle, "try-error")
    
    #if(failed){
    #  browser()
    #}
    
    if (flavor=="superLearning" & failed) {
      tmle <- tmle.npvi(sim1$obs, f=identity, ff=ff, isBetaMult=TRUE, flavor="learning", nMax=40, iter=maxIter, trueGMu=trueGMu, trueTheta=trueTheta)
      attr(tmle, "flag") <- "Flavor 'superLearning' failed, carried out flavor 'learning' instead."
    }  
    
    resultEstimate[[i]] <- tmle
    #resultEstimate[[i]] <- tmle.npvi(sim1$obs, f=identity, ff=ff, isBetaMult=TRUE, flavor=flavor, nMax=40, iter=maxIter, trueGMu=trueGMu, trueTheta=trueTheta)
  }
  
  if(!dontComputeSim){
    result <- list(sim = simulation, estimate = resultEstimate )    
  }
  else{
    result <- list(estimate = resultEstimate )
  }
}



statisticalAnalysisProcedure_TMLE_NPVI <- function(
  d,                ## number of cavariates to be considered
  N = 4,            ## number of procedure to be run      
  b = 300,          ## number of observations to be drawn in order to estimate our parameter of interest using TMLE_NPVI_CA 
  ff = NULL,        ## formula to be used for the formula      
  B = 5e4,         ## number of observations to be created in order to compute the value of the parameter of interest using the simulation scheme
  flavor="learning",## the learning procedure to use for some key variables of our procedures...
  maxIter = 5,      ## maximum number of iterations
  trueGMu = NULL,   ## true function for the varialbes "G" and "Mu"
  trueTheta=NULL,   ## true function for the varaibles "Theta" 
  useFourierData = TRUE
) {    
  
  ##
  ## compute the value of our parameter throught simulation scheme
  ##     
  
  print("Computing parameter of interest value through simulation Scheme...")
  simulation <- NULL
  ##simulation <- getSample( B , O, lambda0,ff=ff, isBetaMult=TRUE, d=d, covariatesMean=rep(0,d), covariatesSigma=diag(d) , p=p, omega=omega, sigma2=1, Sigma3=S)
  ##browser()
  ##
  ##  computing parameter of interest using TMLE_NPVI_CA
  ##
  
  #resultEstimate <- list(N,NA) ;## initializing the variables which will old the entire set of result using "flavor" method 
  resultEstimate <- vector("list", N );
  
  oneEstimationRun <- function(x, b=300, ff=NULL, d=1, maxIter=5 , trueGMu=NULL, trueTheta=NULL){# function of no arguments 
    if(useFourierData){
      typeFunctionF <- "sine_test"   ### option available here are : identity , square, sine
      typeFunctionG <- "identity"  ### option available here are : identity , square, cube
      distribW <- "normal"           ### option available here are : sd_normal, normal, beta, uniform, gamma
      nBasis <- 5
      
      base_fn <- basis_functions(b, lengthBasis = nBasis, d=1, computeFourierOfW=TRUE, distrib_Y = "imply", distrib_W = distribW ,
                                 distrib_X = "uniform", fn_typeF = typeFunctionF, fn_typeG = typeFunctionG)
      fourier_mat <- base_fn$fourier
      obs <- cbind( base_fn$Y, base_fn$X, fourier_mat)
      colnames(obs) <- c( "Y" , "X" , "W", paste("V",1:(nBasis-1),sep=""))
    }
    else{
      sim1 <- getSample(b, O, lambda0, ff=ff, isBetaMult=TRUE, d=d, covariatesMean=rep(0,d), covariatesSigma=diag(d) , p=p, omega=omega, sigma2=1, Sigma3=S)
      obs <- sim1$obs
    }
    
    tmle <- try(tmle.npvi(obs, f=identity, ff=ff, isBetaMult=TRUE, flavor=flavor, nMax=40, iter=maxIter, trueGMu=trueGMu, trueTheta=trueTheta) )
    failed <- inherits(tmle, "try-error")
    
    if (flavor=="superLearning" & failed) {
      tmle <- tmle.npvi(sim1$obs, f=identity, ff=ff, isBetaMult=TRUE, flavor="learning", nMax=40, iter=maxIter, trueGMu=trueGMu, trueTheta=trueTheta)
      attr(tmle, "flag") <- "Flavor 'superLearning' failed, carried out flavor 'learning' instead."
    }  
    
    print(tmle)
    #resultEstimate <- tmle.npvi(sim1$obs, f=identity, ff=ff, isBetaMult=TRUE, flavor="learning", nMax=30, iter=5, trueGMu=NULL, trueTheta=NULL)
    resultEstimate
  }
  
  estimationRuns <- function( nbRuns , nbCores ){
    require(parallel) 
    
    if(Sys.info()[1] == "Windows"){
      cl <- makeCluster(nbCores)
      result <- clusterApply(cl=cl, x =1:nbRuns, fun=oneEstimateRun , b=b, ff=ff, d=d, maxIter=maxIter,trueGMu=trueGMu,trueTheta=trueTheta )
      stopCluster(cl)
    }
    else{
      require(SnowballC)
      require(tm)
      result <- mclapply(X=1:nbRuns , FUN=oneEstimationRun, b=b, ff=ff, d=d, maxIter=maxIter,trueGMu=trueGMu,trueTheta=trueTheta , mc.cores=nbCores)
    }
    result
  }
  
  print("starting to compute estimate of parameter of interest using TMLE method...")
  resultEstimate <- estimationRuns(nbRuns = N, nbCores = 1)  ### should be 3 or 4 !!!!!! 
  
  result <- list(sim = simulation, estimate = resultEstimate )    
  #result <- list( estimate = resultEstimate )  
}

