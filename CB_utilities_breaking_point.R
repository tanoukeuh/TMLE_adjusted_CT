#####################################################################################################################
##
##  Cabral CHANANG - May 2017
## 
#####################################################################################################################
##  Description of the file :
##
##  1. Will be used to find the breaking point of a component of the function f.
##  2. By "breaking point", we refer here to the value attached to a given component which makes our estimator
##     incapable of capturing the signal coming from that component.
##
#####################################################################################################################



run_breakPoint <- function(ff,gg){
  
  ptm <- proc.time()
  
  ## break point position
  bkPtMax <- 1.5
  bkPtMin <- 0
  bkPtSq <- .02
  bkPtPos <- 14

  ## few more parameters
  nbCases <- 6
  step <- 75
  method <- 2
  lengthBaseX <- 25
  lengthBaseW <- 25
  typeFunctionF <- ff   ### option available here are : identity , square, sine, sine_test, sine_test_small, breaking_point
  typeFunctionG <- gg   ### option available here are : identity , square, cube , zero, g_sine, g_sine_small
  distribW <- "beta"    ### option available here are : sd_normal, normal, beta, uniform, gamma
  nbW <- 1
  mcCores <- 4
  subFnMcCores <<- 4
  
  
  mainResult_bkPoint <- startComputation_forBreakingPoint(bkPtMax, bkPtMin, bkPtSq, bkPtPos, nbCases, step, lengthBaseX, lengthBaseW, method, typeFunctionF, typeFunctionG, distribW, nbW, mcCores );
  
  
  ## adding element of the run to overall result
  mainResult_bkPoint$bkPtMax <- bkPtMax
  mainResult_bkPoint$bkPtMin <- bkPtMin
  mainResult_bkPoint$bkPtSq <- bkPtSq
  mainResult_bkPoint$bkPtPos <- bkPtPos
  
  mainResult_bkPoint$nbCase <- nbCases
  mainResult_bkPoint$step <- step
  mainResult_bkPoint$lengthBaseX <- lengthBaseX
  mainResult_bkPoint$lengthBaseW <- lengthBaseW
  mainResult_bkPoint$method <- method
  mainResult_bkPoint$typeFunctionF <- typeFunctionF
  mainResult_bkPoint$typeFunctionG <- typeFunctionG
  mainResult_bkPoint$distribW <- distribW
  mainResult_bkPoint$nbW <- nbW
  
  # saving the file containing all results
  sRunDescription = paste("method=", method , "_bkPtMax=", bkPtMax , "_bkPtMin=", bkPtMin , "_bkPtPos=", bkPtPos ,"_bkPtSq=", bkPtSq, "_Step=",step, "_lengthBaseX=", lengthBaseX,
                          "_lengthBaseW=",lengthBaseW, "_typeFunctionF=",typeFunctionF, "_typeFunctionG=",typeFunctionG,
                          "_distribW=", distribW, sep="")
  filePath = paste(sRunDescription, ".RData")
  save(mainResult_bkPoint, file=filePath)
  
  ## plotting sparsity marker
  breakPointSeq <- seq(bkPtMin , bkPtMax , bkPtSq)
  plotBreakingPointAnalysis(mainResult_bkPoint$matrix_sparsity, sRunDescription, breakPointSeq, step)
  
  ## test
  print( proc.time() - ptm )
  ##test
  
  ##
  ## removing all the local variables
  ##
  rm(ptm, mainResult_bkPoint, breakPointSeq)
}


##
## 1. we only focus here on the method 2, which has produced strong result up to this point
## 2. The number of repetitions will be replaced by the number of breaking points that we attempts to calculate. We will not need
##    multiple repetitions per breaking point because for a single repetition, the simulation will be done n time where n is our sample size
##
startComputation_forBreakingPoint <- function(bkPtMax, bkPtMin, bkPtSq, bkPtPos, nbCases, step, lengthBaseX, lengthBaseW, 
                             method , typeFunctionF, typeFunctionG, distribW, nbW, mcCores){
  
  ## few checks
  if( bkPtMax <= bkPtMin){
    throw("the max braeking point is smaller or equal to the min breaking point !")
  }
  
  if(bkPtSq >= bkPtMax){
    throw("The sequence of breaking point is greater or equal to the max breaking point value")
  }
  
  breakPointValues <- seq(bkPtMin, bkPtMax, bkPtSq)
  nbRepetitions <- length(breakPointValues)
  
  result <- list()
  result_sparsity <- matrix(0, ncol = nbRepetitions, nrow = nbCases)
  
  for(i in 1:nbCases){
    print(paste("I am running the case numero " , i, sep=""))
    nBasis <- i * step
    
    values <- 1:nbRepetitions
    fun <- function(i){
      if(method == 1){
        throw("The first method is not considered for this exercise !!! ")
      }
      else{
        intermediaryTotalResult <- estimation_Of_function_fW_2_Prediction_error(nBasis, lengthBaseX, lengthBaseW, typeFunctionF, typeFunctionG, distribW, nbW = nbW, mcCores=mcCores, breakingFactor = breakPointValues[i])
      }
      return(intermediaryTotalResult)
    }
    
    
    partialRes <- mclapply(values, fun, mc.cores = mcCores)
    
    
    for(j in 1:nbRepetitions){
      interm <- partialRes[[j]]
      matrix_theta <- interm$estimate_theta
      result_sparsity[i,j] <- length(which(matrix_theta[,bkPtPos]==0))/nrow(matrix_theta)
    }
    
    argName <- paste("case-",i,sep="")
    result[[argName]] <- partialRes
  }
  result$matrix_sparsity <- result_sparsit
  
  ##
  ## removing all the local variables
  ##
  rm(result_sparsity, partialRes)
  
  ##
  ## returning final result
  ##
  y
  return(result)
}




## plotting the matrix result of the breaking point above.
plotBreakingPointAnalysis <- function(result_matrix, sRunDescription, bkPtSequence, initialSampleSize){
  
  ## Plotting the result of the method
  pathFileName <- paste( "Evolution_of_BreakPoint_" , sRunDescription, "_file.pdf", sep="")
  pdf(file=pathFileName)
  
  nCases <- ncol(result_matrix)
  nSampleSize <- nrow(result_matrix)
  values <- 1:nSampleSize
  
  ## first iteration
  plot(x = bkPtSequence, y=result_matrix[1,], type="o", col=1, xlab="break point sequence", ylab="sparsity percentage")
  
  ## plotting the rest of the sequence
  for(i in 2:nSampleSize){
    lines(bkPtSequence, result_matrix[i,], type="o", col=i)
  }
  
  ## creating title and legend
  legend("topleft", paste(" Sample Size : ", initialSampleSize*values , sep=""), col=c(1:nSampleSize), cex=.5,lwd=3)
  title(main= paste("Evol. estimate of function of Sparsity given sample size"))
  
  dev.off()
}
