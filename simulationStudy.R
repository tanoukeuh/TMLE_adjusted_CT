## updating workspace directory
currentDirectory <- "~/Desktop/phd_thesis/project_with_cristina/problem/code/"
setwd(currentDirectory)
#sourceDirectory(getwd())



## Loading essential packages
library(abind)
library(acepack)
library(arm)
library(class)
library(cluster)
library(datasets)
library(e1071)
library(fBasics)
library(fMultivar)
library(foreign)
library(Formula)
library(ggplot2)
library(graphics)
library(grDevices)
library(grid)
library(Hmisc)
library(lattice)
library(lme4)
library(lpSolve)
library(magic)
library(MASS)
library(Matrix)
library(methods)
library(mnormt)
library(mvtnorm)
library(nlme)
library(nloptr)
library(nnls)
library(parallel)
library(R.methodsS3)
library(R.oo)
library(R.utils)
library(Rcpp)
library(SparseM)
library(stats)
library(survival)
library(SuperLearner)
library(timeSeries)
library(utils)
library(timeDate)
library(geometry)  
library(gdata)


## defining a formula
ff <- formula(Y~I(W)+I(V1)+I(V2)+I(V3)+I(V4)-1)
##ff <- formula(Y~I(W)+I(W^2))

## function to use for antoine in order to run 1000 replications of the learning and super learning Method
AntoineSimulationRuns <- function( case = 1){
  if(case==1){
    result <- StartProcedure(N=10, b=1000, ff=ff, maxIter=4, parallel = TRUE)    ### learning case with no trueGMu or trueTheta
  }
  else if(case==2){
    result <- StartProcedure(N=1000, b=300, ff=ff, maxIter=4, trueGMu=trueGMu) ### learning case with trueGMu
  }
  else if(case==3){
    result <- StartProcedure(N=1000, b=300, ff=ff, maxIter=4, trueTheta=trueTheta) ### learning case with trueTheta
  }
  else if(case == 4){
    result <- StartProcedure(N=1000, b=300, ff=ff, maxIter=3, flavor="superLearning") ### superLearning case with no trueGMu or trueTheta
  }
  else if(case == 5){
    result <- StartProcedure(N=1000, b=300, ff=ff, maxIter=3 , flavor="superLearning" , trueGMu=trueGMu) ### superLearning case with trueGMu
  }
  else if(case == 6){
    result <- StartProcedure(N=1000, b=300, ff=ff, maxIter=3 , flavor="superLearning" , trueTheta = trueTheta) ### superLearning case with trueTheta
  }

  result
}

### File to run in order to set properly the entire environment
StartProcedure <- function(N = 4, b=300, ff=ff, trueGMu=NULL, trueTheta=NULL, flavor="learning", maxIter = 3, parallel=FALSE){  
  
  
  ## Parameters for the simulation (case 'f=identity')
  O <<- cbind(W=c(0.05218652, 0.01113460),
             X=c(2.722713, 9.362432),
             Y=c(-0.4569579, 1.2470822))
  O <<- rbind(NA, O)
  lambda0 <<- function(W) {-W}
  p <<- c(0, 1/2, 1/2)
  omega <<- c(0, 3, 3)
  S <<- matrix(c(10, 1, 1, 0.5), 2 ,2)
  f <<- identity
  Sigma1 <<-matrix(c(1, 1/sqrt(2), 1/sqrt(2), 1), 2, 2)
  Sigma3 <<- Sigma1
  
  ### proceduce
  if(parallel){
    SIM <- statisticalAnalysisProcedure_TMLE_NPVI(d=1, b=b, N = N, ff=ff, trueGMu=trueGMu, trueTheta=trueTheta , maxIter = maxIter, flavor=flavor)  
  }
  else{    
    SIM <- statisticalAnalysisProcedure_TMLE_NPVICA(d=1, b=b, N = N, ff=ff, trueGMu=trueGMu, trueTheta=trueTheta , maxIter = maxIter, flavor=flavor, dontComputeSim=FALSE)  
  }
  ##  title definition
  partTitle <- ""
  if(!is.null(trueGMu)){
    partTitle <- "trueGMu"
  }  
  else if(!is.null(trueTheta)){
    partTitle <- "trueTheta"
  }
  
  tradesInIC <-  function(SIM){
    truth <- SIM$sim$psi
    d <- length(truth)
    INF <- SIM$estimate
    validINF <-  sapply(INF, function(ll){!is.character(ll)}) 
    PSI <- sapply(INF, function(ll){ll$psi})    
    
    ##
    ## initializing list needed for confidence interval
    ##
    nameCovariates <- names(SIM$sim$psi)
    interval <- list(NA,NA)
    test <- list(NA,NA)
    result <- matrix(NA, nrow=1, ncol=d)
    colnames(result) <- names(SIM$estimate[[1]]$history$psi[[1]])
    
    ##
    ## construction of interval and test
    ##    
    
    nObs <- nrow(getObs(SIM$estimate[[1]]))
    nbIterations <- length(INF)
    resultType <- getFlavor(SIM$estimate[[1]])
    
    for(i in  1:d){
      alpha <- 1-getConfLevel(SIM$estimate[[i]])
      interval[[i]] <- sapply(INF, function(ll, idx=i){ll$.psi[idx]  +  c(-1,1) * qnorm(1 - alpha/2) * ll$.psi.sd[idx] })
      test[[i]] <- sapply(1:ncol(interval[[i]]) , function(ii){interval[[i]][1,ii]<= truth[i] & interval[[i]][2,ii] >= truth[i] })
      result[1,i] <- paste0( format( 100* sum(test[[i]])/length(test[[i]]), format="f", digits=3), "%" )
      
      
      ##
      ## plotting graph corresponding to the final values of the estimation regardless of the number of iterations steps
      ##      
      finalResult <- sapply(INF, function(ll, idx=i){ll$psi[idx]})
      d <- density(unlist(finalResult))
      ylimDensity <- range(d$y)
      xlimDensity <- range(d$x)
      
      fileName = paste(nameCovariates[i], " - Final Empirical Density - ", resultType , " - ", partTitle   ,".pdf ", sep="")
      pdf(file=fileName)
      plot( d ,  type="o", lty = 1 , pch = 1 , col=1 , main=NA, xlab=NA, ylab=NA , xlim= xlimDensity, ylim=ylimDensity)      
      
      ## adding the result given by the simulation
      abline(v=truth[i], col = 5)          
      title(main = paste( nameCovariates[i] , "- Final Density ", partTitle, " ( N = ", nObs , " )" , " - ", "\n", nbIterations ," Iter. of " , resultType , "\n" , result[1,i] , " of CIs contain the true value.", sep=""))
      dev.off()
    }
    
    result <- list(interval=interval, test=test, result=result)
    result
  }
  
  ## plotting the evolution of the density per iterations step.
  graph <- plot_statisticalAnalysis_NPVICA(SIM$estimate, SIM$sim, maxIter= maxIter , partTitle=partTitle)
  
  
  IC <- tradesInIC(SIM)
  
  SIM <- list(sim=SIM, IC=IC)
  SIM
}

##  ------ in Urgence ---------####


## SUGGESTIONS - TO DO ##

### CODE
### vérifier la façon dont la mémoire est gérée
### plutôt: ne sauver que ce qui est pertinent dans 'statisticalAnalysisProcedure_TMLE_NPVICA'
###
### SIMULATIONS
### 1)
### jouer sur les formules: formula(Y~W), formula(Y~W+I(W^2))
### 2)
### jouer sur lambda0
### 3)
### estimation avec 'g' et 'mu' correctement spécifiés
### estimation sans aucune spécification correcte
### 4)
### jeu sur la taille de l'échantillon d'apprentissage (2e2, 1e3)
### 5)
### "learning" et "super-learning
###
### APPLICATION AUX DONNÉES RÉELLES!

##  ------ in Urgence ---------####

####
####
#tradesInIC2 <-  function(SIM){
#  truth <- SIM$sim$psi
#  d <- length(truth)
#  INF <- SIM$estimate
#  validINF <-  sapply(INF, function(ll){!is.character(ll)}) 
#  INF <- INF[validINF]
#  
#  PSI <- sapply(INF, function(ll){ll$psi}) 
#  
#  ##
#  ## initializing list needed for confidence interval
#  ##
#  nameCovariates <- names(SIM$sim$psi)
#  interval <- list(NA,NA)
#  test <- list(NA,NA)
#  result <- matrix(NA, nrow=1, ncol=d)
#  colnames(result) <- names(SIM$estimate[[1]]$history$psi[[1]])
#  
#  ##
#  ## construction of interval and test
#  ##    
#  nObs <- nrow(getObs(SIM$estimate[[1]]))
#  nbIterations <- length(INF)
#  resultType <- getFlavor(SIM$estimate[[1]])
#  
#  for(i in  1:d){
#    alpha <- 1-getConfLevel(SIM$estimate[[i]])
#    interval[[i]] <- sapply( INF, function(ll, idx=i){ll$.psi[idx]  +  c(-1,1) * qnorm(1 - alpha/2) * ll$.psi.sd[idx] })
#    test[[i]] <- sapply(1:ncol(interval[[i]]) , function(ii){interval[[i]][1,ii]<= truth[i] & interval[[i]][2,ii] >= truth[i] })
#    result[1,i] <- paste0( format( 100* sum(test[[i]])/length(test[[i]]), format="f", digits=3), "%" )
#    
#    ##
#    ## plotting graph corresponding to the final values of the estimation regardless of the number of iterations steps
#    ##      
#    finalResult <- lapply( INF,  function(ll, idx=i){ll$.psi[idx]})
#    
#    d <- density(unlist(finalResult))
#    ylimDensity <- range(d$y)
#    xlimDensity <- range(d$x)
#    
#    fileName = paste(nameCovariates[i], " - Final Empirical Density - TrueGMu.pdf ", sep="")
#    pdf(file=fileName)
#    plot( d ,  type="o", lty = 1 , pch = 1 , col=1 , main=NA, xlab=NA, ylab=NA , xlim= xlimDensity, ylim=ylimDensity)      
#    
#    ## adding the result given by the simulation
#    abline(v=truth[i], col = 5)          
#    title(main = paste( nameCovariates[i] , "- Final Density with trueTheta ( N = ", nObs , " )" , " - ", "\n", nbIterations ," Iter. of " , resultType , sep=""))
#    dev.off()
#  }
#  
#  result <- list( interval=interval, test=test, result=result)
#  result
#}
####
####



### Plotting Antoine Results 
plotAntoineResult <- function(  RES , SIM_RESULT , charact, confLevel =0.95){  
    
    ##  title definition
    partTitle <- ""  
    if(!is.null(charact$trueGMu)){
      partTitle <- "trueGMu"
    }  
    else if(!is.null(charact$trueTheta)){
      partTitle <- "trueTheta"
    }
      
    truth <- SIM_RESULT$psi
    d <- length(truth)
    
    #filtering out all the bad results
    INF <- RES
    INF <- INF[lapply(INF, is.character)==FALSE]    
     
    ##
    ## initializing list needed for confidence interval
    ##
    nameCovariates <- names(SIM_RESULT$psi)
    interval <- list(NA,NA)
    test <- list(NA,NA)
    result <- matrix(NA, nrow=1, ncol=d)  
    
    colnames(result) <- nameCovariates
    
    ##
    ## construction of interval and test
    ##    
    
    nObs <- charact$nObs
    nbIterations <- length(RES)
    resultType <- charact$flavor
    
    for(i in  1:d){
      alpha <- 1-confLevel      
      interval[[i]] <- sapply(INF, function(ll, idx=i){ll$psi[idx]  +  c(-1,1) * qnorm(1 - alpha/2) * sqrt(diag(ll$varCovar/nObs))[idx] }) 
      test[[i]] <- sapply(1:ncol(interval[[i]]) , function(ii){interval[[i]][1,ii]<= truth[i] & interval[[i]][2,ii] >= truth[i] }) 
      result[1,i] <- paste0( format( 100* sum(test[[i]])/length(test[[i]]), format="f", digits=3), "%" ) 
      
      #interval[[i]] <- sapply(INF, function(ll, idx=i){ll[idx]  +  c(-1,1) * qnorm(1 - alpha/2) * ll[length(ll)/2 + idx] }) 
      #interval[[i]] <- sapply(INF, function(ll, idx=i){ll[idx]  +  c(-1,1) *  sqrt( qchisq( 1-alpha/2 , df=d) * ll[length(ll)]) })  
      #test[[i]] <- sapply(1:ncol(interval[[i]]) , function(ii){interval[[i]][1,ii]<= truth[i] & interval[[i]][2,ii] >= truth[i] }) 
      #b <- 10
      #test[[i]] <- sapply(1:ncol(interval[[i]]) , function(ii){   round( interval[[i]][1,ii] , digits = b)   <= round( truth[i], digits = b) &  round( interval[[i]][2,ii]  , digits=b   ) >=   round( truth[i] , digits=b)    })
      
      ##
      ## plotting graph corresponding to the final values of the estimation regardless of the number of iterations steps
      ##      
      finalResult <- sapply(INF, function(ll, idx=i){ll$psi[idx]})
      den <- density(unlist(finalResult))
      ylimDensity <- range(den$y)
      xlimDensity <- range(den$x)
      
      fileName = paste(nameCovariates[i], " - Final Empirical Density - ", resultType , " - ", partTitle   ,".pdf ", sep="")
      pdf(file=fileName)
      plot( den ,  type="o", lty = 1 , pch = 1 , col=1 , main=NA, xlab=NA, ylab=NA , xlim= xlimDensity, ylim=ylimDensity)      
      
      ## adding the result given by the simulation
      abline(v=truth[i], col = 5)          
      title(main = paste( nameCovariates[i] , "- Final Density ", partTitle, " ( N = ", nObs , " )" , " - ", "\n", nbIterations ," Iter. of " , resultType , "\n" , result[1,i] , " of CIs contain the true value.", sep=""))
      dev.off()
    }   
}






### Plotting Antoine Results 
plotAntoineResult2 <- function(  RES , SIM_RESULT , charact, confLevel =0.95){  
  
  ##  title definition
  partTitle <- ""  
  if(!is.null(charact$trueGMu)){
    partTitle <- "trueGMu"
  }  
  else if(!is.null(charact$trueTheta)){
    partTitle <- "trueTheta"
  }
  
  truth <- SIM_RESULT$psi
  d <- length(truth)
  
  #filtering out all the bad results
  INF <- RES
  INF <- INF[lapply(INF, is.character)==FALSE]    
  
  if( length(RES) != length(INF)){
    browser()
  }
  
  ##
  ## initializing list needed for confidence interval
  ##
  nameCovariates <- names(SIM_RESULT$psi)
  interval <- list(NA,NA)
  test <- list(NA,NA)
  result <- matrix(NA, nrow=1, ncol=d)  
  
  colnames(result) <- nameCovariates
  
  ##
  ## construction of interval and test
  ##    
  
  nObs <- charact$nObs
  nbIterations <- length(RES)
  resultType <- charact$flavor
  
  
  ### Constructing the different confidence intervals.
  if(TRUE){
    
    ### construct a single confidence interval
    ### Sigma corresponds to a single line of the matrix of results constructed.
    constructSingleConfidenceInterval <- function(sigmaAndMu, confLevel, truth , n){
      
      d <- length(truth)
      
      ## Extracting and constructin the true sigma matrix
      EndIndexSigma <- length(sigmaAndMu)
      StartIndexSigma = length(sigmaAndMu) - d^2 + 1
      trueSigma <- matrix( sigmaAndMu[StartIndexSigma:EndIndexSigma], ncol=d, nrow=d,   byrow = TRUE )
      Mu <- sigmaAndMu[1:d]
      computeConfidenceInterval(trueSigma, Mu, n , confLevel)
    }    
    
    interval <- sapply(INF, function(ll){   constructSingleConfidenceInterval(ll, confLevel, truth, nObs)      })   
    
  }
  
  for(i in  1:d){
    
    #alpha <- 1-confLevel 
    
    #interval[[i]] <- sapply(INF, function(ll, idx=i){ll[idx]  +  c(-1,1) * qnorm(1 - alpha/2) * ll[length(ll)/2 + idx] })   
        
    ## test[[i]] <- sapply(1:ncol(interval[[i]]) , function(ii){interval[[i]][1,ii]<= truth[i] & interval[[i]][2,ii] >= truth[i]})   
    pseudoInterval <- t(interval)
    
    test[[i]] <- sapply(  1:length(pseudoInterval[,1]) , function(ii){      pseudoInterval[  ii , 2*i - 1 ]  <=  truth[i]        &          pseudoInterval[ ii, 2*i  ] >= truth[i]        }      )       
    
    result[1,i] <- paste0( format( 100* sum(test[[i]])/length(test[[i]]), format="f", digits=3), "%" )            
    
    ##
    ## plotting graph corresponding to the final values of the estimation regardless of the number of iterations steps
    ##      
    finalResult <- sapply(INF, function(ll, idx=i){ll[idx]})
    den <- density(unlist(finalResult))
    ylimDensity <- range(den$y)
    xlimDensity <- range(den$x)
    
    fileName = paste(nameCovariates[i], " - Final Empirical Density - ", resultType , " - ", partTitle   ,".pdf ", sep="")
    pdf(file=fileName)
    plot( den ,  type="o", lty = 1 , pch = 1 , col=1 , main=NA, xlab=NA, ylab=NA , xlim= xlimDensity, ylim=ylimDensity)      
    
    ## adding the result given by the simulation
    abline(v=truth[i], col = 5)          
    title(main = paste( nameCovariates[i] , "- Final Density ", partTitle, " ( N = ", nObs , " )" , " - ", "\n", nbIterations ," Iter. of " , resultType , "\n" , result[1,i] , " of CIs contain the true value.", sep=""))
    dev.off()
  }   
}





### Plotting Antoine Results 
plotAntoineResult3 <- function(  RES , SIM_RESULT , charact, confLevel =0.95){  
  
  ##  title definition
  partTitle <- ""  
  if(!is.null(charact$trueGMu)){
    partTitle <- "trueGMu"
  }  
  else if(!is.null(charact$trueTheta)){
    partTitle <- "trueTheta"
  }
  
  truth <- SIM_RESULT$psi
  d <- length(truth)
  
  #filtering out all the bad results
  INF <- RES
  INF <- INF[lapply(INF, is.character)==FALSE]    
  
  if( length(RES) != length(INF)){
    browser()
  }
  
  ##
  ## initializing list needed for confidence interval
  ##
  nameCovariates <- names(SIM_RESULT$psi)
  interval <- list(NA,NA)
  test <- list(NA,NA)
  result <- matrix(NA, nrow=1, ncol=d)  
  
  colnames(result) <- nameCovariates
  
  ##
  ## construction of interval and test
  ##    
  
  nObs <- charact$nObs
  nbIterations <- length(RES)
  resultType <- charact$flavor
  
  
  ### Constructing the different confidence intervals.
  if(TRUE){
    
    ##
    ## Test if multi dimensional point is in the confidence region
    ##
    
    constructTestConfidenceRegion <- function(  sigmaAndMu, confLevel, n, truePsi){      
      result <- testIfInRegion(  sigmaAndMu$psi, truePsi ,  sigmaAndMu$varCovar, n, confLevel)
    }    
    
    testOfInclusion <- sapply(INF, function(ll){ constructTestConfidenceRegion( ll, confLevel, nObs , truth ) })
    result_test <- sum(testOfInclusion)/length(testOfInclusion)
    print(result_test)
    browser()
    
    
    ### construct a single confidence interval
    ### Sigma corresponds to a single line of the matrix of results constructed.
    constructSingleConfidenceInterval <- function(sigmaAndMu, confLevel , n){
          
      ## Extracting and constructin the true sigma matrix      
      trueSigma <- sigmaAndMu$varCovar
      Mu <- sigmaAndMu$psi      
      computeConfidenceInterval(trueSigma, Mu, n , confLevel)
    }    
    
    interval <- sapply(INF, function(ll){   constructSingleConfidenceInterval(ll, confLevel, nObs)      })   
    
  }
  
  for(i in  1:d){
    
    #alpha <- 1-confLevel 
    
    #interval[[i]] <- sapply(INF, function(ll, idx=i){ll[idx]  +  c(-1,1) * qnorm(1 - alpha/2) * ll[length(ll)/2 + idx] })   
    
    ## test[[i]] <- sapply(1:ncol(interval[[i]]) , function(ii){interval[[i]][1,ii]<= truth[i] & interval[[i]][2,ii] >= truth[i]})   
    pseudoInterval <- t(interval)
    
    test[[i]] <- sapply(  1:length(pseudoInterval[,1]) , function(ii){      pseudoInterval[  ii , 2*i - 1 ]  <=  truth[i]        &          pseudoInterval[ ii, 2*i  ] >= truth[i]        }      )       
    
    result[1,i] <- paste0( format( 100* sum(test[[i]])/length(test[[i]]), format="f", digits=3), "%" )            
    
    ##
    ## plotting graph corresponding to the final values of the estimation regardless of the number of iterations steps
    ##      
    
    finalResult <- sapply(INF, function(ll, idx=i){ll$psi[idx]})
    den <- density(unlist(finalResult))
    ylimDensity <- range(den$y)
    xlimDensity <- range(den$x)
    
    fileName = paste(nameCovariates[i], " - Final Empirical Density - ", resultType , " - ", partTitle   ,".pdf ", sep="")
    pdf(file=fileName)
    plot( den ,  type="o", lty = 1 , pch = 1 , col=1 , main=NA, xlab=NA, ylab=NA , xlim= xlimDensity, ylim=ylimDensity)      
    
    ## adding the result given by the simulation
    abline(v=truth[i], col = 5)          
    title(main = paste( nameCovariates[i] , "- Final Density ", partTitle, " ( N = ", nObs , " )" , " - ", "\n", nbIterations ," Iter. of " , resultType , "\n" , result[1,i] , " of CIs contain the true value.", sep=""))
    dev.off()
  }   
}



###
### build cluster from the TCGA simulation study results
###


buildCluster <- function(RES, iClusterSize = 2){
  set.seed(12345)
  ##
  ## determine the length of each variable 'psi' computed for the simulation
  ##
  
  sizeEachPsi <- sapply(RES , function(ll){  length(ll$psi)   })
  nameEachPsi <- sapply(RES , function(ll){  ll$name          })
  
  ##
  ## clustering only genes having the same size
  ##
  
  minSize <- min(sizeEachPsi)
  maxSize <- max(sizeEachPsi)
  clusterSize <- iClusterSize  ## arbitrary number of the cluster we are trying to set up....
  result <- list()
  result2 <- list()
  iterator <- 1
  
  for(i in minSize:maxSize){
    indices <- which(sizeEachPsi==i)
    if(length(indices) > 1){
      clusterList <- list()
      clusterList[["position"]] <- paste("size - ", i, sep="")
        
      psi_list <- sapply(indices, function(ii){ RES[[ii]]$psi})
      
      genes_names <- sapply(indices, function(ii){ RES[[ii]]$name })
      psiNames <- names(RES[[indices[1]]]$psi)
      
      # matrixPsi <- matrix(unlist(psi_list), ncol=length(psiNames), byrow=TRUE)
      # colnames(matrixPsi) <- psiNames
      
      matrixPsi <- t(psi_list)
      rownames(matrixPsi) <- genes_names
      clusterList[["genesNames"]] <- genes_names
      
      for( j in 1:i ){
        X <- matrixPsi[,j]
        cl <- kmeans(X, clusterSize)
        clusterList[[psiNames[j]]] <- cl
        plot(X, col = cl$cluster)    
      }
      
      commonNames <- matchGenesNames(clusterList, clusterSize, psiNames)      
      iterator <- iterator + 1
    }
    
    namesIterator <- paste( iterator , "-pos", sep="")
    #result[[namesIterator]] <- clusterList
    result[[namesIterator]] <- commonNames
  }
  
  result  
}



##
## return a list, if any, of mathcing genes within the same cluster area
##


matchGenesNames <- function( clusterList, sizeCluster , clusterNames ){
  
  nbOfParameters <- length(clusterNames)
  results <- list()
  
  for( i in 1:sizeCluster){
    finalNames <- names( which(clusterList[[ clusterNames[1]  ]]$cluster==i) )
    for ( j in 2:nbOfParameters){
      
      intermNames <- names( which(clusterList[[clusterNames[j]]]$cluster==i))
      finalIndices <- match(  intermNames , finalNames)
      finalIndices <- finalIndices[!is.na(finalIndices)]
      finalNames <- finalNames[finalIndices]      
    }
    results[[paste("cluster number ", i, sep="")]] <- finalNames        
  }  
  
  results
}

##
## creating a plot positioning the the parameter of interest on a graph
##

plotParameterOfInterestReduction <- function(RES  , option = 1){
  
  parOfInterestVal <- t(sapply(RES, function(ll){ ll$psi}))
  posGeneOnChromosome <- sapply(RES, function(ll){ as.double( strsplit(ll$name , ",")[[1]][2] )  })
  
  if( option ==1 ){
    dataEllipse(parOfInterestVal[,1], parOfInterestVal[,2], levels=c(0.95), ylab="W", xlab="intercept" )
  }
  else if(option == 2){
    plot(parOfInterestVal , ylab="W" )
  }
  else if(option == 3){
    plot( x = posGeneOnChromosome , y = parOfInterestVal[,1] , ylab="Intercept" , xlab = "position of Gene")
  }
  else if(option == 4){
    plot( x = posGeneOnChromosome , y = parOfInterestVal[,2] , ylab="W" , xlab = "position of Gene")
  }
  else if(option == 5){
    #plot(0,0, , xlim=range(parOfInterestVal[,1]) , ylim=range(parOfInterestVal[,2]) ) ## initialization of the plot window - not very pretty - to change eventually
    plot(0.5,0.5, ylim=c(-10,10), xlim=c(0,1) )
    browser()
    nIter = nrow(parOfInterestVal)
    for( i in 1:nIter){
      abline(parOfInterestVal[i,1]  ,  parOfInterestVal[i,2] ) 
    }
  }
  #title("confidence region (50% and 95%) of the parameter of interest points in our TCGA dataset")
  #plot(parOfInterestVal , ylab="W" )
  #title("Parameter of Interest of Genes in our TCGA Dataset")
}


processResultBAAC <- function(RES){
  ##
  ## creating a matrix containing all our estimators
  ##
  dummyVal <- c(0,0,0)
  psi <- sapply(RES, function(ll){ll$psi})
  result <- t(sapply(1:length(psi), function(ll){ if(is.null(psi[[ll]])){psi[[ll]] <- dummyVal }else{psi[[ll]]}}) )
  colnames(result) <- names(RES[[3]]$psi)
  ##
  ## Creating a column containing the birth year of the driver and the number of elements used for the estimation
  ##
  
  dummy_birthYear <- 0
  birth_year <- sapply(RES, function(ll){ if( is.null(ll$birth_year)){ dummy_birthYear }else{ll$birth_year}})
  
  dummy_dataSize <- 0
  data_size <- sapply(RES, function(ll){ if( is.null(ll$dataSize)){ dummy_dataSize }else{ll$dataSize}})
  
  result <- cbind( result , "birthYear" = birth_year)
  result <- cbind( result , "dataSize" = data_size)
  
  ##
  ## clearning data and removing all which covariates estimates have been set to zero
  ##

  index <- which( result[,"V1"] != 0 )
  result <- result[index,]
  result
  
  
  ##
  ## plotting a graph per covariates
  ##
  covariatesNames <- c("Sexe", "Weather Condition", "Road Quality")
  birthYear <- result[,"birthYear"]
  
  for( i in 1:length(covariatesNames)){
    fileName = paste(covariatesNames[i], " - Scatter plot and linear regression", ".pdf", sep="")
    pdf(file=fileName)
    
    covariateData <- result[,i]
    plot( birthYear , covariateData, col= 1 , xlab="driver birth year" , ylab = "Causal Value")
    
    ##
    ## linear regression
    ##
    linearRegression <- lm(formula = covariateData ~ birthYear)
    ## adding linear regression to the above plot
    abline(linearRegression, col = 4)
    
    title(main = paste( covariatesNames[i] ))
    dev.off()
  }
  
}




### Plotting Antoine Results 
plotClimateChangeResult <- function(  RES , trueCovariatesNames = NULL ,confLevel =0.95 , startingPoint = 1){  
  
  ## covariateNames <- c("Under-five mortality rate", "Population growth","Urban population growth (annual)","CO2 emissions per capita (metric tons)", "CO2 emissions per units of GDP (kg/$1,000 of 2005 PPP $)","Energy use per capita (kilograms of oil equivalent)","Energy use per units of GDP (kg oil eq./$1,000 of 2005 PPP $)" )
  ## {Y,X}  <- c("Projected annual precipitation change (2045-2065, mm)","CO2 emissions per capita (metric tons)")
  
  #filtering out all the bad results
  INF <- RES
  INF <- INF[unlist(lapply(INF, function(ll){is.double(ll$psi[1])}))]
  #INF <- INF[lapply(INF, is.character)==FALSE]    
  
  ##
  ## initializing list needed for confidence interval
  ##
  if(!is.null(trueCovariatesNames)){
    nameCovariates <- trueCovariatesNames
  }
  else{
    nameCovariates <- names(INF[[1]]$psi)
  }

  interval <- list(NA,NA)
  #result <- matrix(NA, nrow=1, ncol=d)  
  #browser()
  #colnames(result) <- nameCovariates
  
  ##
  ## construction of interval and test
  ##    
  nbIterations <- length(nameCovariates)
  nbCovariates <- length(INF)
  
  for(i in  1:nbIterations){
    alpha <- 1-confLevel   
    bonferroniCorrection <- length(startingPoint:nbCovariates)
    #ll$dataSize[idx]
    #dataSize <- 160
    #interval <- sapply(INF, function(ll, idx=i){ll$psi[idx]  +  c(-1,1) * qnorm(1 - alpha/2) * sqrt(diag(ll$varCovar)/dataSize)[idx] }) 
    interval <- sapply(INF, function(ll, idx=i){ll$psi[idx]  +  c(-1,1) * qnorm(1 - alpha/(2*bonferroniCorrection)) * ll$sd[idx] }) 
    psiCov <- sapply(INF, function(ll, idx=i){ll$psi[i]})
    years <- sapply(INF, function(ll){ll$year})
    
    ## Only selecting the variables i need for the graphs    
    interval <- interval[,startingPoint:nbCovariates]
    psiCov <- psiCov[startingPoint:nbCovariates]
    years <- years[startingPoint:nbCovariates]
    
    #years <- c(1990,1991,1992,1993,1994,1995,1996,1997,1998)
    #years <- c(1990,1991)
    
    ##
    ## plotting the results by year with their corresponding confidence interval
    ##
    plotType = ".pdf"
    fileName = paste(nameCovariates[i], " - Evolution with confidence Interval",plotType, sep="")
    
    pdf(file=fileName)
    n = length(psiCov) 
    ylim <- range(interval[2,], interval[1,])
    #ylim <- c(-.1,.3)
    #browser()
    #errbar(years, psiCov, interval[1,], interval[2,] , type="o", lty = 1 , pch=1, col=1, ylab="Values" , xlab = "iterations", asp=2,
    #       ylim=ylim)    
    errbar(years, psiCov, interval[1,], interval[2,] , type="o",  ylab="Value" , xlab = "year of the analysis", ylim=ylim)    
    title(main= paste( nameCovariates[i] , " - ", " with Confidence Interval", sep="")   ) 
    dev.off()
    
    ##
    ## plotting graph corresponding to the final values of the estimation regardless of the number of iterations steps
    ##     
    finalResult <- sapply(INF, function(ll, idx=i){ll$psi[idx]})
    den <- density(unlist(finalResult))
    ylimDensity <- range(den$y)
    xlimDensity <- range(den$x)
    
    fileName = paste(nameCovariates[i], " - Final Empirical Density" ,plotType, sep="")
    pdf(file=fileName)
    plot( den ,  type="o", lty = 1 , pch = 1 , col=1 , main=NA, xlab=NA, ylab=NA , xlim= xlimDensity, ylim=ylimDensity)      
    
    ## adding the result given by the simulation
    title(main = paste( nameCovariates[i] , "- Final Density", sep=""))
    dev.off()
  }   
}
