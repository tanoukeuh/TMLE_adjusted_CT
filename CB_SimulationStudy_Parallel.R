#############################################################################################################
##
##  Cabral CHANANG - December 2016
## 
#############################################################################################################
##  Description of the file :
##
##  will be used as a triggering mechanism for our simulation.
##
#############################################################################################################


##
## updating workspace directory
##

currentDirectory <- "~/Desktop/phd_thesis/project_with_cristina/problem/code/"
setwd(currentDirectory)
#sourceDirectory(getwd())

## Loading essential packages
library(abind)
library(class)
library(datasets)
library(fBasics)
library(fMultivar)
library(ggplot2)
library(graphics)
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
library(R.oo)
library(R.utils)
library(SparseM)
library(stats)
library(survival)
library(SuperLearner)
library(timeSeries)
library(utils)
library(timeDate)
library(geometry)  
library(gdata)
library(fda)
library(splines)
library(Hmisc)
library(lars)
library(glmnet)
library(parallel)
library(SnowballC)
library(doMC)
library(flare)

registerDoMC(cores=4)


##
## simple test of my procedure
##
simple_test <- function(){
  
  test <- proc.time()
  
  run_tmle_algo <- TRUE
  runLassoWithError <- TRUE
  ff <- "sine"
  gg <- "identity"
  lengthBaseW <- 20
  lengthBaseX <- 20
  nbCores <- 4
  
  # with error
  mu <- 0.01
  run_parallel(ff,gg,runLassoWithError=runLassoWithError, mu=mu, nbCores = nbCores, lengthBaseW = lengthBaseW, lengthBaseX = lengthBaseX)
  
  #tmle procedure
  #run_parallel(ff,gg,run_tmle_algo=run_tmle_algo, nbCores = nbCores, lengthBaseW = lengthBaseW, lengthBaseX = lengthBaseX)
  
  print(proc.time() - test)
}

##
## running all the simulation
##

run_all_procedure <- function(){
  nbCores <- 4
  
  ## with error simulation
  run_parallel_lassoWithError_lowDimension(nbCores)
  run_parallel_lassoWithError_highDimension(nbCores)
  
  ## TMLE comparison
  run_parallel_tmle_lowDimension(nbCores)
  run_parallel_tmle_highDimension(nbCores)
}


##
## run the tmle simulation
##
run_parallel_tmle_lowDimension <- function(nbCores){
  run_tmle_algo <- TRUE
  ff <- "sine"
  gg <- "identity"
  lengthBaseW <- 20
  lengthBaseX <- 20
  nbCores <- nbCores
  
  run_parallel(ff,gg,run_tmle_algo=run_tmle_algo, nbCores = nbCores, lengthBaseW = lengthBaseW, lengthBaseX = lengthBaseX)
}

run_parallel_tmle_highDimension <- function(nbCores){
  run_tmle_algo <- TRUE
  ff <- "sine"
  gg <- "identity"
  lengthBaseW <- 200
  lengthBaseX <- 50
  nbCores <- nbCores
  
  run_parallel(ff,gg,run_tmle_algo=run_tmle_algo, nbCores = nbCores, lengthBaseW = lengthBaseW, lengthBaseX = lengthBaseX)
}

##
## run the simulation case with error
##

run_parallel_lassoWithError_lowDimension <- function(nbCores){
  
  test <- proc.time()
  
  run_tmle_algo <- FALSE
  runLassoWithError <- TRUE
  ff <- "sine_test"
  gg <- "identity"
  lengthBaseW <- 20
  lengthBaseX <- 20
  nbCores <- nbCores
  listMu <- seq(0.05,.15,.05)
  
  for(i in 1:length(listMu)){
    mu <- listMu[i]
    run_parallel(ff,gg,runLassoWithError=runLassoWithError, mu=mu, nbCores = nbCores, lengthBaseW = lengthBaseW, lengthBaseX = lengthBaseX)
  }
  
  print(proc.time() - test)
}

run_parallel_lassoWithError_highDimension <- function(nbCores){
  
  test <- proc.time()
  
  run_tmle_algo <- FALSE
  runLassoWithError <- TRUE
  ff <- "sine_test"
  gg <- "identity"
  lengthBaseW <- 200
  lengthBaseX <- 50
  nbCores <- nbCores
  listMu <- seq(0.00,.15,.05)
  
  for(i in 1:length(listMu)){
    mu <- listMu[i]
    run_parallel(ff,gg,runLassoWithError=runLassoWithError, mu=mu, nbCores = nbCores, lengthBaseW = lengthBaseW, lengthBaseX = lengthBaseX)
  }
  
  
  print(proc.time() - test)
}

##
## start of the script
##
run_parallel <- function(ff, gg, run_tmle_algo = FALSE, runLassoWithError = FALSE, mu=NULL , nbCores = 4, lengthBaseW = 5, lengthBaseX = 5){
  
  ### test
  #ptm <- proc.time()
  ###test
  
  ##
  ##  validity of parameters
  ##
  if(is.null(mu) && runLassoWithError == TRUE){
    throw("The parameter mu can not be set to TRUE and the lasso with Error set to TRUE as well.")
  }

  ###
  ###   Key Parameters of the simulation
  ###
  ## TMLE and Coordinate descent addition
  run_tmle_algo <- run_tmle_algo
  runLassoWithError <- runLassoWithError
  mu <- mu
  mcCores <- nbCores
  
  isLasso <- TRUE
  isLassoPenalized <- FALSE
  nbRepetitions <- 20
  nbCases <- 4
  step <- 75
  method <- 2
  lengthBaseX <- lengthBaseX
  lengthBaseW <- lengthBaseW
  typeFunctionF <- ff   ### option available here are : identity , square, sine, sine_test, sine_test_small
  typeFunctionG <- gg   ### option available here are : identity , square, cube , zero, g_sine, g_sine_small
  distribW <- "normal"    ### option available here are : sd_normal, normal, beta, uniform, gamma
  nbW <- 1
  subFnMcCores <<- 4
  maxFourier_basis <<- 4
  
  mainResult <- startComputation(nbRepetitions,nbCases,step, lengthBaseX, lengthBaseW, method, typeFunctionF, typeFunctionG, distribW, nbW, mcCores, isLasso, isLassoPenalized, run_tmle_algo, runLassoWithError, mu );
  
  ## adding element of the run to overall result
  mainResult$nbRepetitions <- nbRepetitions
  mainResult$nbCase <- nbCases
  mainResult$step <- step
  mainResult$lengthBaseX <- lengthBaseX
  mainResult$lengthBaseW <- lengthBaseW
  mainResult$method <- method
  mainResult$typeFunctionF <- typeFunctionF
  mainResult$typeFunctionG <- typeFunctionG
  mainResult$distribW <- distribW
  mainResult$nbW <- nbW
  mainResult$isLasso <- isLasso
  mainResult$isLassoPenalized <- isLassoPenalized
  
  ## TMLE data
  mainResult$run_tmle_algo <- run_tmle_algo
  ## LASSO with in error variable
  mainResult$runLassoWithError <- runLassoWithError
  mainResult$mu <- mu
  
  # plotting the evolution of error
  error <- mainResult$error
  if(run_tmle_algo){
    error_tmle <- mainResult$error_tmle
    #riskPrediction_comparison_Result_Plotter(error, error_tmle, step, lengthBaseX, lengthBaseW, method, nbCases, nbRepetitions, typeFunctionF, typeFunctionG, distribW, nbW);
  }
  else{
    #riskPrediction_Result_Plotter(error, step, lengthBaseX, lengthBaseW, method, nbCases, nbRepetitions, typeFunctionF, typeFunctionG, distribW, nbW);
  }
  
  ## plotting the evolution of sparsity
  #sparsity <- mainResult$sparsity
  #riskPrediction_Result_Plotter(sparsity, step, lengthBaseX, lengthBaseW, method, nbCases, nbRepetitions, typeFunctionF, typeFunctionG, distribW, nbW, isSparsity=TRUE);
    
  ## plotting the evolution of error of the function G
  #errorG <- mainResult$errorG
  #riskPrediction_Result_Plotter(errorG, step, lengthBaseX, lengthBaseW, method, nbCases, nbRepetitions, typeFunctionF, typeFunctionG, distribW, nbW, isFunctionG=TRUE);
  
  # saving the file containing all results
  sTMLE <- ""
  if(run_tmle_algo){
    sTMLE <- "_with_TMLE"
  }
  
  sWithError <- ""
  if(runLassoWithError){
    sWithError <- paste("_with_Error_", mu*100, "%", sep="")
  }
  
  sRunDescription = paste("method=", method ,"_nbRepetitions=",nbRepetitions,"_Step=",step, "_lengthBaseX=", lengthBaseX,
                          "_lengthBaseW=",lengthBaseW, "_typeFunctionF=",typeFunctionF, "_typeFunctionG=",typeFunctionG,
                          "_distribW=", distribW, "_nbW=", nbW, sTMLE, sWithError, sep="")
  filePath = paste(sRunDescription, ".RData")
  save(mainResult, file=filePath)
  
  ### test
  # print( proc.time() - ptm )
  ###test
  
  ##
  ## removing all local variables
  ##
  rm(mainResult)
}

startComputation <- function(nbRepetitions, nbCases, step, lengthBaseX, lengthBaseW, 
                             method , typeFunctionF, typeFunctionG, distribW, nbW, mcCores,isLasso = TRUE, isLassoPenalized = FALSE,
                             run_tmle_algo = FALSE, runLassoWithError = FALSE, mu=NULL){
  error <- matrix(0, ncol=2, nrow =nbCases)
  error_tmle <- matrix(0, ncol=2, nrow =nbCases)
  errorG <- matrix(0, ncol=2, nrow=nbCases)
  sparsity <- matrix(0, ncol=2, nrow =nbCases)
  result <- list()
  
  for(i in 1:nbCases){
    print(paste("I am running the case numero " , i, "out of ", nbCases, sep=""))
    intermediary_error <- rep(0, nbRepetitions)
    intermediary_error_tmle <- rep(0, nbRepetitions)
    intermediary_errorG <- rep(0, nbRepetitions)
    intermediary_sparsity <- rep(0, nbRepetitions)
    nBasis <- i * step
    
    values <- 1:nbRepetitions
    fun <- function(i){
      if(method == 1){
        intermediaryTotalResult <- predictionRiskCalculator(nBasis, lengthBaseX, lengthBaseW, typeFunctionF, typeFunctionG, distribW, nbW=nbW, mcCores=mcCores)
      }
      else{
        intermediaryTotalResult <- estimation_Of_function_fW_2_Prediction_error(nBasis, lengthBaseX, lengthBaseW, typeFunctionF, typeFunctionG, distribW, nbW = nbW, mcCores=mcCores, isLasso= isLasso, isLassoPenalized=isLassoPenalized, runTMLE_Algo = run_tmle_algo, runLassoWithError=runLassoWithError, mu=mu)
      }
      return(intermediaryTotalResult)
    }

    partialRes <- mclapply(values, fun, mc.cores = mcCores)
    
    for(j in 1:nbRepetitions){
      interm <- partialRes[[j]]
      intermediary_error[j] <- interm$result
      intermediary_errorG[j] <- interm$resultG
      intermediary_sparsity[j] <- interm$sparsity
      if(run_tmle_algo){
        intermediary_error_tmle[j] <- intermediaryTotalResult$result_tmle
      }
    }

    error[i,1] <- mean(intermediary_error)
    error[i,2] <- sd(intermediary_error)
    
    if(run_tmle_algo){
      error_tmle[i,1] <- mean(intermediary_error_tmle)
      error_tmle[i,2] <- sd(intermediary_error_tmle)
    }
    
    errorG[i,1] <- mean(intermediary_errorG)
    errorG[i,2] <- sd(intermediary_errorG)
    
    sparsity[i,1] <- mean(intermediary_sparsity)
    sparsity[i,2] <- sd(intermediary_sparsity)
    
    argName <- paste("case-",i,sep="")
    result[[argName]] <- partialRes
  }

  #result$main <- main
  result$error <- error
  if(run_tmle_algo){
    result$error_tmle <- error_tmle
  }
  result$errorG <- errorG
  result$sparsity <- sparsity
  
  ##
  ## removing all local variables
  ##
  rm(partialRes, intermediary_sparsity, intermediary_error, intermediary_errorG)
  
  ## return final variables
  return(result)
}


##
## Computing the error corresponding to the approximation of the function f using the brut force methodology
##
predictionRiskCalculator <- function(nBasis, lengthBaseX, lengthBaseW, typeFunctionF, typeFunctionG, distribW, nbW , mcCores ){

  base_fn <- basis_functions(nBasis, d=nbW, distrib_Y = "imply", distrib_W = distribW, 
                             distrib_X = "uniform", fn_typeF = typeFunctionF, fn_typeG = typeFunctionG)
  fourier_mat <- base_fn$fourier
  spline_mat <- base_fn$spline
  Y <- base_fn$Y
  W <- base_fn$W
  X <- base_fn$X
  
  result_error <- rep(0, nBasis)
  result_sparsity <- rep(0,nBasis)
  values <- 1:nBasis
  
  ## result <- least_square_solution(Y, fourier_mat, spline_mat)
  iterX <- lengthBaseX      #25
  iterW <- lengthBaseW      #25
  
  fun_list_Of_ThetaMatrix <- function(i){
    Y_minus_i <- Y[-i,]
    W_minus_i <- W[-i,]
    X_minus_i <- X[-i,]
    fourier_mat_minus_i <- fourier_mat[-i,]
    spline_mat_minus_i <- spline_mat[-i,]
    result <- least_square_solution_Theta(Y_minus_i, fourier_mat_minus_i, spline_mat_minus_i, iterX, iterW)
    return(result)
  }
  
  # list of theta matrix
  listOfThetaMatrix <- mclapply(values, fun_list_Of_ThetaMatrix, mc.cores = mcCores)
  
  # list of sparsity
  fun_list_of_sparsity <- function(i){
    sparse <- length(which(listOfThetaMatrix[[i]]==0))/length(listOfThetaMatrix[[i]])
    return(sparse)
  }
  listOfSparsity <- mclapply(values, fun_list_of_sparsity, mc.cores = mcCores)
  
  fun_list_of_result <- function(i){
    fourier_mat_minus_i <- fourier_mat[-i,]
    spline_mat_minus_i <- spline_mat[-i,]
    
    ## computing the estimate of f_w using the ith values
    deriv_fourier_basis_i <- derivative_fourier_extract(X[i,], nBasis)
    usedSpline_mat_i <- spline_mat[i,1:iterW]
    result <- listOfThetaMatrix[[i]]
    estimate_fW_i <- computeEstimateY(result, deriv_fourier_basis_i[1,1:iterX], usedSpline_mat_i)
    
    ## computing result error
    W_i <- matrix(W[i,], nrow=1, ncol = length(W[i,]))
    result_error_fun <- (f_W(W_i, typeFunctionF) - estimate_fW_i)^2
  
    return(result_error_fun)
  }
  listOfError <- mclapply(values, fun_list_of_result, mc.cores = mcCores)
  
  #final result
  listResult <- list()
  listResult$base_fn <- base_fn
  listResult$thetaMatrix <- listOfThetaMatrix
  listResult$result <- sum(unlist(listOfError))/nBasis
  listResult$sparsity <- mean(unlist(listOfSparsity))
  return(listResult)
}


##
## computation of the estimate of the function f_W, using the first methodology of our paper
##

## When user is trying to run LASSO procedure with error-in-variable, the value for "mu" must be provided.

## Computing the error corresponding to the approximation of the function f 
estimation_Of_function_fW_2_Prediction_error <- function(nBasis,lengthBaseX, lengthBaseW, typeFunctionF, typeFunctionG, distribW, nbW, mcCores, 
                                                         breakingFactor = 1, isLasso = TRUE, isLassoPenalized = FALSE,
                                                         runTMLE_Algo = FALSE, runLassoWithError = FALSE, mu = NULL){ # not sure yet how to acheive fourier for W in  R^2
  
  ## checking validity of some variables.
  if(runLassoWithError){
    if(is.null(mu)){
      throw("[estimation_Of_function_fW_2_Prediction_error]:the variable mu must be provided if we want to run LASSO with in error variable.")
    }
  }
  
  ## Instead of relying on the number of element in the sample, we will rely on the user defined basis for our simulation
  ## hence nBasis will be different..
  lengthBasis <- lengthBaseX + lengthBaseW
  
  ## add an option here such that the number of basis is reduce significantly when looking at dimension greater than 2 to avoid speed issue
  base_fn <- basis_functions(nBasis, lengthBasis, d=nbW, computeFourierOfW=TRUE, distrib_Y = "imply", distrib_W = distribW ,
                             distrib_X = "uniform", fn_typeF = typeFunctionF, fn_typeG = typeFunctionG, breakingFactor = breakingFactor, isWithError = runLassoWithError, mu = mu)
  
  
  fourier_mat <- base_fn$fourier
  if(nbW==1){
    fourier_mat <- fourier_mat[,1:lengthBasis]
  }
  
  Y <- base_fn$Y
  W <- base_fn$W
  X <- base_fn$X
  
  if(runLassoWithError){
    Z <- base_fn$Z
  }
  
  iterW1 <- lengthBaseX   #30
  iterW2 <- lengthBaseW   #30
  
  result_matrix_theta <- matrix(0, nrow = nBasis, ncol = iterW2)
  result_matrix_eta <- matrix(0, nrow = nBasis, ncol = iterW1)
  result_matrix_beta <- matrix(0, nrow = nBasis, ncol = iterW2+iterW1)
  values <- 1:nBasis
  
  # compute the theta
  fun_listOfBeta <- function(i){
    Y_minus_i <- Y[-i,]
    X_minus_i <- X[-i,]
    fourier_mat_minus_i <- fourier_mat[-i,]
    
    if(runLassoWithError){
      Z_minus_i <- Z[-i,]
    }
    
    if(runLassoWithError){
      # estimate_beta <- cv.lasso_coordinate_descent_cpp(Y=Y_minus_i, z=Z_minus_i, mu=mu, Phi_mat=fourier_mat_minus_i, iterW1=iterW1, iterW2=iterW2, mcCores = mcCores)
      nbFolds <- 3
      estimate_beta <- rep(0,iterW1 + iterW2)
      indicesForTest <- sample(length(Y_minus_i))
      folds <- cut(seq(1,length(indicesForTest)), breaks=nbFolds, labels=FALSE)
      
      cv_lasso_coordinate_descent_all_cpp(estimate_beta , Y_minus_i , Z_minus_i, mu, fourier_mat_minus_i, indicesForTest, folds,
                                      iterW1,  iterW2, nbFolds);
      
      
      
    }
    else{
      estimate_beta <- least_square_solution_beta(Y_minus_i, X_minus_i , fourier_mat_minus_i, iterW1 , iterW2 , isLasso, isLassoPenalized)
    }
    
    # we will return the full vector beta and post process it
    # estimate_theta <- estimate_beta[(iterW1+1):(length(estimate_beta))]
    estimate_beta <- estimate_beta[1:(length(estimate_beta))]
    return(estimate_beta)
  }
  list_Of_Beta <- mclapply(values, fun_listOfBeta, mc.cores = mcCores)
  

  
  #compute the matrix of beta (corresponds to function f) and eta (corresponds to function g)
  result_matrix_beta <- matrix(unlist(list_Of_Beta), ncol=iterW1 + iterW2, byrow = TRUE)
  
  ## the matrix of eta corresponds to the (iterW1 + 1) first columns of the matrix above
  result_matrix_eta <- result_matrix_beta[,1:iterW1]
  
  ## the matrix of theta corresponds to the iterW2 last elements of hte matrix of beta
  length_beta <- ncol(result_matrix_beta)
  result_matrix_theta <- result_matrix_beta[,(iterW1+1):length_beta]
  
  #compute the prediction error
  fun_listOfError <- function(i){
    #browser()
    estimate_theta <- result_matrix_theta[i,]
    #estimate_fW_i <- X[i,]*( fourier_mat[i,(iterW1+1):(iterW1+iterW2)] %*% estimate_theta)
    estimate_fW_i <- fourier_mat[i,1:iterW2] %*% estimate_theta
    
    if(nbW == 1){
      result_error_fun <- (estimate_fW_i - f_W(W[i,], typeFunctionF, breakingFactor))^2
    }
    else{
      result_error_fun <- (estimate_fW_i - f_W(   matrix( W[i,], ncol=length(W[i,])) , typeFunctionF, breakingFactor))^2
    }
    return(result_error_fun)
  }
  list_Of_Error <- mclapply(values, fun_listOfError, mc.cores = mcCores)
  
  
  fun_listOfError_G <- function(i){
    estimate_eta <- result_matrix_eta[i,]
    estimate_gW_i <- fourier_mat[i,1:iterW1] %*% estimate_eta
    result_error_fun <- (estimate_gW_i - g_W(W[i,], typeFunctionG))^2
  }
  list_Of_Error_g <- mclapply(values, fun_listOfError_G, mc.cores = mcCores)
  
  
  #compute the sparsity
  fun_listOfSparsity <- function(i){
    sparse <- length(which(result_matrix_theta[i,]==0))/length(result_matrix_theta[i,])
    return(sparse)
  }
  list_Of_Sparse <- mclapply(values, fun_listOfSparsity, mc.cores = mcCores)
  
  
  #browser()
  #final result
  result <- list()
  result$base_fn <- base_fn
  result$estimate_theta <- result_matrix_theta
  result$estimate_eta <- result_matrix_eta
  result$result <- sum(unlist(list_Of_Error))/nBasis
  result$resultG <- sum(unlist(list_Of_Error_g))/nBasis
  result$sparsity <- mean(unlist(list_Of_Sparse))
  
  ##
  ## removing some local variables
  ##
  rm(list_Of_Beta, list_Of_Sparse, list_Of_Error_g, list_Of_Error, result_matrix_beta, result_matrix_eta, result_matrix_theta)
  
  ##
  ## running the TMLE based algorithm if requested by the user
  ##
  if(runTMLE_Algo){
    fourier_mat_tmle <- fourier_mat[,1:iterW2]
    maxIter_tmle <- 5  ### to increase to 5 at a latter
    flavor <- "learning"
    result_error_tmle <- rep(0, nBasis)
    result_matrix_tmle <- matrix(0, nrow = nBasis, ncol = iterW2)
    B <- 1e4   ### parameter use to generate estimate using the TMLE method.
    ## computation of the formula ff
    ff_string <- paste( "Y~I(W)+" , paste( paste("I(V",1:(iterW2-1), ")",sep=""), sep="", collapse="+"),"-1", sep="", collapse="")
    ff <- formula(ff_string)
    
    oneEstimationRun <- function(increment){# function of no arguments
      Y_minus_i <- Y[-increment,]
      X_minus_i <- X[-increment,]
      fourier_mat_minus_i <- fourier_mat_tmle[-increment,]
      
      obs <- cbind(Y_minus_i, X_minus_i, fourier_mat_minus_i)
      colnames(obs) <- c( "Y" , "X" , "W", paste("V",1:(iterW2-1),sep=""))
      estimate <- tmle.npvi(obs, f=identity, ff=ff, isBetaMult=TRUE, flavor=flavor, nMax=40, iter=maxIter_tmle, B = B, trueGMu=NULL, trueTheta=NULL)
      return(estimate)
    }
    
    require(parallel)
    require(SnowballC)
    require(tm)
    list_of_result <- mclapply(X=1:nBasis , FUN=oneEstimationRun , mc.cores=mcCores)  ### reduced to  1 for a test... should be 4
    
    
    
    ### computing the error associated with the TMLE estimate
    fun_listTMLE_ERROR <- function(i){
      estimate_theta_tmle <- list_of_result[[i]]$psi
      
      ## calculating the error
      estimate_fW_i <-  fourier_mat_tmle[i,]%*%estimate_theta_tmle
      result_error_tmle <- (estimate_fW_i - f_W(W[i,], typeFunctionF))^2
      
      ## assignment variables
      result_matrix_tmle[i,] <- estimate_theta_tmle
      
      ## return the error value
      return(result_error_tmle)
    }
    
    result_error_tmle_list <- mclapply(X=1:nBasis, FUN=fun_listTMLE_ERROR, mc.cores = mcCores)
    
    result$result_tmle <- sum(unlist(result_error_tmle_list))/nBasis
    result$result_matrix_tmle <- result_matrix_tmle
    
    ## removing all local variables
    rm(result_error_tmle_list, result_matrix_tmle, list_of_result)
  }
  
  return(result)
}
