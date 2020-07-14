################################################################################################################
##
##  Cabral CHANANG - January 2017
## 
################################################################################################################
##  Description of the file :
##
##  The following file will contain all the necessary functions as well as processo to obtain solution of
##  of our problem using the third methodology defined in the working paper
##
################################################################################################################


n <- 50
#marginal_result <- marginal_based_theta(n)
#plot_marginal_variables(marginal_result$marginalVariables$eta)
#plot_error(marginal_result)


##
## marginal function in charge of computing and plotting the final result
##
marginal_based_theta <- function(n, bFourier=FALSE){
  result <- basis_functions(n, isComingFromMarginalComputation = TRUE)
  W <- result$W
  X <- result$X
  Y <- result$Y
  
  if(bFourier){
    Phi_W <- result$spline
  }
  else{
    Phi_W <- result$fourier
  }
  
  ## we need to attribute columns names to Phi_W
  n <- ncol(Phi_W)
  colnames_Phi <- paste("Phi",1:n, sep="")
  colnames(Phi_W) <- colnames_Phi
  
  result1 <- compute_theta_through_marginal_expectation(W,X,Y,Phi_W )
  
  mainResult <- list()
  mainResult$basisResult <- result
  mainResult$marginalVariables <- result1
  
  ##
  ## removing local variables
  ##
  rm(result, result1)
  
  ## returning final result
  return(mainResult)
}

##
##  Based on the third approach of our working document
##  E[Y|W] = Phi(W)   (   eta +   E[X|W] * theta) 
##  E[Y|X]  = E[Phi(W)|X] ( eta +  X * theta   )
##  the solution of the above set of equations will allow us to compute theta
## 
##  p: size of the function 
##
compute_theta_through_marginal_expectation <- function(W,X,Y, Phi_W){
  
  ## number of elements in our simulation
  n <- length(X)
  
  ## number of columns of our basis function will define the size of the variables \theta and \eta
  p <- ncol(Phi_W)
  
  ## we will use half of the elements for simulation and the other half for projection
  m <- n - floor(n/2)
  E_Phi_W_X <- matrix(0, ncol=p, nrow = m)
  E_Y_X <- rep(0, m)
  E_Y_W <- rep(0,m)
  E_X_W <- rep(0,m)
  
  
  W1 <- matrix( W[1:(n-m),] , ncol=ncol(W) , nrow = n-m )
  W2 <- matrix( W[(n-m+1):n,], ncol=ncol(W) , nrow = m )
  colnames(W1) <- colnames(W)
  colnames(W2) <- colnames(W)
  
  Y1 <- matrix( Y[1:(n-m),1], ncol=1, nrow= n-m)
  colnames(Y1) <- colnames(Y)
  
  X1 <- matrix( X[1:(n-m),1], ncol=1, nrow = n-m)
  colnames(X1) <- colnames(X)
  
  X2 <- matrix( X[(n-m+1):n,1] , ncol=1, nrow=m )
  colnames(X2) <- colnames(X)
  
  E_Y_W_fun <- learn_Esperance_Y_Knowing_X(Y1,W1)
  E_Y_X_fun <- learn_Esperance_Y_Knowing_X(Y1,X1)
  E_X_W_fun <- learn_Esperance_Y_Knowing_X(X1,W1)
  
  
  E_Y_W <- E_Y_W_fun(W2)
  E_Y_X <- E_Y_X_fun(X2)
  E_X_W <- E_X_W_fun(W2)
  
  ## starting calculation of different key elements
  for(i in 1:p){
    Phi_W_i_1 <- matrix( Phi_W[1:(n-m),i] , ncol=1, nrow = n-m)
    colnames(Phi_W_i_1) <- colnames(Phi_W)[i]
    E_Phi_W_X_fun <- learn_Esperance_Y_Knowing_X(Phi_W_i_1,X1)
    E_Phi_W_X[,i] <- E_Phi_W_X_fun(X2)
  }
  
  ## initialisation of the variables \eta and \theta and extracting proper values for Phi_W
  eta <- matrix(0, ncol=p, nrow = m)
  theta <- matrix(0,ncol=p, nrow = m)
  Phi_W2 <- Phi_W[(n-m+1):n,]
  
  for(i in 1:p){
    for(j in 1:m){
      theta[j,i] <- (1/( X2[j] - E_X_W[j])) * ((E_Y_X[j]/E_Phi_W_X[j,i]) - (E_Y_W[j]/Phi_W2[j,i])) 
      eta[j,i] <- (E_Y_W[j]/Phi_W2[j,i]) - E_X_W[j] * theta[j,i]
    }
  }
  
  result <- list()
  result$theta <- theta
  result$eta <- eta
  return(result)
}






##
## we will attempt to calculate below the measure E[Y|X] assuming here that both X and Y are vectors
## the variable "obs", will contain both the variables X and Y
##
learn_Esperance_Y_Knowing_X <- function (Y,X, ...){
  
  varNames <- colnames(X)
  varNamesY <- colnames(Y)
  
  theFormula <- paste(varNames, collapse = " + ")
  
  theFormula2 <- paste("I(", varNames, "^2)", collapse = " + ", 
                         sep = "")
  theFormula <- paste(varNamesY ," ~", theFormula, "+", theFormula2, 
                        sep = " ")
  
  #theFormula <- paste("Y ~", theFormula)
  
  formula <- as.formula(theFormula)
  
  obs <- cbind(Y,X)
  fit <- glm(formula, data = as.data.frame(obs), family = gaussian)
  
  foo <- function(W) {
    predict(fit, newdata = data.frame(W), type = "response")
  }
  
  
  ##
  ## removing local variables
  ##
  rm(fit, obs)
  
  ## returning final result
  return(foo)
}