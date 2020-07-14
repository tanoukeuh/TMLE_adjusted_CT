#############################################################################################################
##
##  Cabral CHANANG - October 2017
## 
#############################################################################################################
##  Description of the file :
##
##  Implementation of a bespoke coordinate descent in order to implement to solve our optimization problem
##  with error-in-variable
##
##  url: http://jocelynchi.com/a-coordinate-descent-algorithm-for-the-lasso-problem
##
#############################################################################################################

##
## soft thresholding function
##

##
## if we consider the equation V(B_l) = a*B_l^2 + b*B_l  + \lambda * |B_l|
## then through the resolution of an unidimensional LASSO equation using "sub-differential", we can conclude that
## B_l  =  \frac{-1}{2a} *  (\lambda + b) if \lambda < -b
## B_l  =  \frac{1}{2a} *  (\lambda + b)  if \lambda < b
## B_l  = 0                               if  b \in [ -\lambda, \lambda   ]

soft_thresholding <- function(a , b , lambda){
  ## checking the validity of some variables
  if(lambda < 0){
    print("i am here")
    throw("[soft_thresholding]:the variable lambda must be positive")
  }
  
  if(a == 0){
    throw("[soft_thresholding]:the variable a must be different from zero")
  }
  
  beta = 0
  if( lambda < -b){
    beta = -(lambda + b)/(2*a)
  }
  else if( lambda < b){
    beta = (lambda - b)/(2*a)
  }
  else if( -lambda <= b && b <= lambda){
    beta = 0
  }
  else{ ## should never get here -- just being bizarre.
    throw("[soft_thresholding]:the variable provided are not conformed...")
  }
  
  return(beta)
}

##
## checking if the solution in hand has converged. To achieve that, we check if it is within the sub differential space
##

##
## let's consider the equation \beta = \argmin \| Y - Z * \beta \|_2^2  +  \mu^2 * \| H * \beta \|_2^2 + \lambda \| \beta \|_1
## we want to check whether a set \hat{\beta} is an optimal solution of the equation above
## define by g(\beta) = 2*Z^t*( Y - Z * \beta ) - 2*\mu^2*H^t*(H * \beta)
## hence \beta is a solution when : forall k, g_k \in [-\lambda, \lambda] if \beta_k = 0 and g_k = \lambda * sign(\beta_k) if \beta_k \neq 0
##

subDifferentialConditionForBeta <- function( Y , Z , beta, mu , H , lambda, tolerance = 1e-3){
  ## checking the validity of some variables
  
  if(lambda<0){
    throw("[subDifferentialConditionForBeta]:The variable lambda must be positive")
  }
  
  if(ncol(Z) != length(beta)){
    throw("[subDifferentialConditionForBeta]:The number of column of the matrix Z must be equal to the length of beta")
  }
  
  if(ncol(H) != length(beta)){
    throw("[subDifferentialConditionForBeta]:The number of colums of the matrix H must be equal to the length of beta")
  }
  
  ## computing the vector g
  g <- 2* ( t(Z)%*%( Y - Z%*%beta) ) - 2*(mu^2)* ( t(H)%*%(H%*%beta) )
  
  ## checking if conditions are respected or not.
  singleConditionCheck <- function(i){
    if(i > length(g)){
      throw("[subDifferentialConditionForBeta]:the variable i must be lower than the length of g")
    }
    
    g_i <- g[i]
    beta_i <- beta[i]
    result <- FALSE
    if( abs(beta_i) > tolerance ){ ## beta_i is assumed to be different to zero here
      if( abs( g_i - lambda*sign(beta_i)) < tolerance ){
        result <- TRUE
      }
    }
    else{
      if( (g_i + lambda) > -1*tolerance && ( lambda - g_i) > -1*tolerance ){
        result <- TRUE
      }
    }
    return(result)
  }

  isBetaSolution <- FALSE
  result <- mclapply(X = 1:length(g), singleConditionCheck, mc.cores = 4)   ## to replace by a variable.
  
  if(sum(unlist(result)) == length(g)){
    isBetaSolution <- TRUE
  }
  
  ##
  ## removing local variables
  ##
  rm(result)
  
  ## returning final value
  return(isBetaSolution)
}

##
##  writing the COORDINATE DESCENT algorithm which will help us solve our optimization problem
##

##
## we are trying to solve the equation \beta = \argmin \| Y - Z * \beta \|_2^2  +  \mu^2 * \| H * \beta \|_2^2 + \lambda \| \beta \|_1
##
lasso_coordinate_descent_givenLambda <- function(Y, Z, mu, H, lambda, beta0 = NULL, maxIter = 100, tolerance = 1e-3){
  ##
  ## checking the validity of the input variable
  ##
  ## beta0 is the starting point of the algorithm. When not specified we will assumed that all of its parameter are equal to zero
  if(is.null(beta0)){
    beta0 <- rep(0, ncol(Z))
  }
  else{
    if( length(beta0) != ncol(Z)){
      throw("[lasso_coordinate_descent_givenLambda]: the length of the variable beta0 should have the equal to the number of columns of matrix Z")
    }
    
    if(length(beta0) != ncol(H) ){
      throw("[lasso_coordinate_descent_givenLambda]: the length of the variable beta0 should be equal to the number of columns of the matrix H")
    }
    
    if( length(Y) != nrow(H) || length(Y) != nrow(Z)){
      throw("[lasso_coordinate_descent_givenLambda]: the variable Y, Z and H should have the same number of rows")
    }
  }
  
  ##
  ## starting the coordinate descent algorithm
  ##
  increment <- 1
  p <- length(beta0)
  beta <- beta0
  result <- beta0
  
  while(increment < maxIter){
    for(j in 1:p){
      R_j = Y - ( Z %*% beta  - Z[,j]*beta[j])
      S_j = H %*% beta - H[,j]*beta[j]
      
      vector_A <- Z[,j]^2 + mu^2 * H[,j]^2
      vector_B <- 2 * ( mu^2 * S_j * H[,j]   - R_j * Z[,j])
      
      a <- sum(vector_A)
      b <- sum(vector_B)
      
      newBeta_j <- soft_thresholding( a, b, lambda)
      ## updating beta_j
      beta[j] <- newBeta_j
    }
    
    
    ## storing the last known values as the result
    result <- beta
    
    ## let's check if the new vector beta is a solution
    if(subDifferentialConditionForBeta(Y , Z , beta, mu , H , lambda, tolerance = tolerance)){
      break;
    }
    
    ## incrementation
    increment <- increment + 1
  }
  
  ## if maxIter was reached then returned the last solved solution
  if(increment == maxIter){
    warning("[subDifferentialConditionForBeta]: the maximum iterations value was reached without finding a suitable solution.. returning the last computed value beta")
  }
  
  return(result)
}


##
## coordinate descent algorithm based on cross validation
##

##
## given a dataset, an a sequence for lambda, i would like to determine the lambda which will provide me the best fit possible for my optimization problem
##
cv.lasso_coordinate_descent <- function(Y, z, mu, Phi_mat, iterW1, iterW2, lambda_max = 10, lambda_min = 0, lambda_seq = 1, maxIter = 100, tolerance = 1e-5, nbFolds = 10, mcCores = 1){
  #
  # checking the validity of few variables
  #
  if(lambda_max <0){
    throw("[cv.lasso_coordinate_descent]:the variable lambda_max must be positive")
  }
  
  if(lambda_seq < 0){
    throw("[cv.lasso_coordinate_descent]:the variable lambda_seq must be positive")
  }
  
  if( ncol(Phi_mat) != iterW1 + iterW2 ){
    throw("[cv.lasso_coordinate_descent]:the number of columns of the variable phi_mat must be equal to the sume of variables iterW1 and iterW2")
  }
  
  if( length(Y) != length(z) || length(z)!= nrow(Phi_mat)){
    throw("[cv.lasso_coordinate_descent]:the variable Y, z must have the same length as the number of row of the matrix Phi_mat")
  }
 
  ## computing key relevant variables
  Z <- cbind(Phi_mat[,1:iterW1], z*Phi_mat[,1:iterW2])
  H <- cbind(matrix(0,ncol=iterW1, nrow=length(Y)), Phi_mat[,1:iterW2])
  list_lambda <- seq(lambda_min,lambda_max, lambda_seq)     #### to change - start the sequence at 0
  
  MSE_output <- matrix(0, ncol=length(list_lambda), nrow=nbFolds)
  
  ## creating the different testing and training folds...  -- will be taking care of at a second stage if time allows.
  indicesForTest <- sample(length(Y))
  folds <- cut(seq(1,length(indicesForTest)), breaks=nbFolds, labels=FALSE)
  
  for(i in 1:length(list_lambda)){
    ## computing the relevant variable
    lambda <- list_lambda[i]
    
    computeSingleMSE <- function(j){
      testIndices <- which(folds==j, arr.ind = TRUE)
      testIndices <- indicesForTest[testIndices]
      
      ## defining the training set
      Z_training <- Z[-testIndices,]
      H_training <- H[-testIndices,]
      Y_training <- Y[-testIndices]
      
      #ptm <- proc.time()
      #print("i am here 2")
      beta <- lasso_coordinate_descent_givenLambda(Y_training, Z_training, mu, H_training, lambda, beta0 = NULL, maxIter = maxIter, tolerance = tolerance)
      #print(proc.time() - ptm)
      
      MSE <- (sum((Y[testIndices] - Z[testIndices,] %*%beta)^2)  + mu^2 * sum( (H[testIndices,]%*%beta)^2))/length(Y[testIndices])
      return(MSE)
    }
    
    lambda_based_MSE <- mclapply(X=1:nbFolds, computeSingleMSE,mc.cores = mcCores)
    
    MSE_output[,i] <- unlist(lambda_based_MSE)
  }
  
  
  ## finding the best lambda given the MSE results computed above
  ## it correspond to the lambda which minimize the biais-variance trade off : Baias^2 + Variance
  Biais_Variance_Calculator <- function(i){
    result <- mean(MSE_output[,i])^2 + sd(MSE_output[,i])^2
    return(result)
  }
  
  biaisVariance <- mclapply(X=1:ncol(MSE_output), Biais_Variance_Calculator, mc.cores = mcCores)
  biaisVariance <- unlist(biaisVariance)
  lambdaSolutionIndice <- which(biaisVariance==min(biaisVariance))
  
  ## computing beta using the best lambda we found
  bestLambda <- list_lambda[lambdaSolutionIndice]
  result <- lasso_coordinate_descent_givenLambda(Y, Z, mu, H, bestLambda, beta0 = NULL, maxIter = maxIter, tolerance = tolerance)
  
  ##
  ## removing local variables
  ##
  rm(MSE_output, list_lambda,lambda_based_MSE, biaisVariance)
  
  ## returning final values
  #browser()
  return(result)
}


##
## given a dataset, an a sequence for lambda, i would like to determine the lambda which will provide me the best fit possible for my optimization problem
##
cv.lasso_coordinate_descent_cpp <- function(Y, z, mu, Phi_mat, iterW1, iterW2, lambda_max = 10, lambda_min = 0, lambda_seq = 1, maxIter = 10, tolerance = 1e-1, nbFolds = 10, mcCores = 1){
  #
  # checking the validity of few variables
  #
  if(lambda_max <0){
    throw("[cv.lasso_coordinate_descent]:the variable lambda_max must be positive")
  }
  
  if(lambda_seq < 0){
    throw("[cv.lasso_coordinate_descent]:the variable lambda_seq must be positive")
  }
  
  if( ncol(Phi_mat) != iterW1 + iterW2 ){
    throw("[cv.lasso_coordinate_descent]:the number of columns of the variable phi_mat must be equal to the sume of variables iterW1 and iterW2")
  }
  
  if( length(Y) != length(z) || length(z)!= nrow(Phi_mat)){
    throw("[cv.lasso_coordinate_descent]:the variable Y, z must have the same length as the number of row of the matrix Phi_mat")
  }
  
  ## computing key relevant variables
  Z <- cbind(Phi_mat[,1:iterW1], z*Phi_mat[,1:iterW2])
  H <- cbind(matrix(0,ncol=iterW1, nrow=length(Y)), Phi_mat[,1:iterW2])
  list_lambda <- seq(lambda_min,lambda_max, lambda_seq)     #### to change - start the sequence at 0
  
  MSE_output <- matrix(0, ncol=length(list_lambda), nrow=nbFolds)
  
  ## creating the different testing and training folds...  -- will be taking care of at a second stage if time allows.
  indicesForTest <- sample(length(Y))
  folds <- cut(seq(1,length(indicesForTest)), breaks=nbFolds, labels=FALSE)
  
  for(i in 1:length(list_lambda)){
    ## computing the relevant variable
    lambda <- list_lambda[i]
    
    computeSingleMSE <- function(j){
      testIndices <- which(folds==j, arr.ind = TRUE)
      testIndices <- indicesForTest[testIndices]
      
      ## defining the training set
      Z_training <- Z[-testIndices,]
      H_training <- H[-testIndices,]
      Y_training <- Y[-testIndices]
      
      #ptm <- proc.time()
      print(paste("[Coordinate Descent - CPP]: I am running the case numero " , i, "out of ", length(list_lambda), sep=""))
      beta <- rep(0 , ncol(Z))
      lasso_coordinate_descent_givenLambda_cpp(Y_training, Z_training, mu, H_training, lambda, beta, maxIter = maxIter, tolerance = tolerance)
      #print(proc.time() - ptm)
      
      MSE <- (sum((Y[testIndices] - Z[testIndices,] %*%beta)^2)  + mu^2 * sum( (H[testIndices,]%*%beta)^2))/length(Y[testIndices])
      return(MSE)
    }
    
    lambda_based_MSE <- mclapply(X=1:nbFolds, computeSingleMSE,mc.cores = mcCores)
    
    MSE_output[,i] <- unlist(lambda_based_MSE)
  }
  
  
  ## finding the best lambda given the MSE results computed above
  ## it correspond to the lambda which minimize the biais-variance trade off : Baias^2 + Variance
  Biais_Variance_Calculator <- function(i){
    result <- mean(MSE_output[,i])^2 + sd(MSE_output[,i])^2
    return(result)
  }
  
  biaisVariance <- mclapply(X=1:ncol(MSE_output), Biais_Variance_Calculator, mc.cores = mcCores)
  biaisVariance <- unlist(biaisVariance)
  lambdaSolutionIndice <- which(biaisVariance==min(biaisVariance))
  
  ## computing beta using the best lambda we found
  bestLambda <- list_lambda[lambdaSolutionIndice]
  result <- rep(0 , ncol(Z))
  lasso_coordinate_descent_givenLambda_cpp(Y, Z, mu, H, bestLambda, result, maxIter = maxIter, tolerance = tolerance)
  
  ##
  ## removing local variables
  ##
  rm(MSE_output, list_lambda,lambda_based_MSE, biaisVariance)
  
  ## returning final values
  #browser()
  return(result)
}
