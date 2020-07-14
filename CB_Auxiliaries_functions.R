#############################################################################################################
##
##  Cabral CHANANG - December 2016
## 
##  Implementation of several bases
##
#############################################################################################################
##  Description of the file :
##
##  We will implement two main generators of fourier and spline basis
##
##
#############################################################################################################


##
## return a vector containing spline and fourier basis based on a set of simulated vectors
## following the TMLE method develop by Pierre & Antoine.
##

## d represents the size of the covariate variables

basis_functions <-function(n, lengthBasis = n,  d=1, eigenValSplineBasis = TRUE, spline_df = 4, 
                           isComingFromMarginalComputation = FALSE , 
                           computeFourierOfW = FALSE,
                           distrib_Y = "identity",
                           distrib_W = "uniform",
                           distrib_X = "uniform",
                           fn_typeF = "identity",
                           fn_typeG = "identity",
                           breakingFactor = 1,
                           isWithError=FALSE,
                           mu=NULL){
  
  ## Simulating a data set of "n" i.i.d. observations
  sim <- sample_generator(n,d, distributionTypeY = distrib_Y, distributionTypeW = "uniform", typeFunctionF = fn_typeF, typefunctionG = fn_typeG, breakingFactor = breakingFactor, isWithError=isWithError, mu=mu)
 
  ## defining column names for the variables
  if(isWithError){
    colNamesW <- setdiff(colnames(sim), c("X","Y","Z"))
  }
  else{
    colNamesW <- setdiff(colnames(sim), c("X","Y"))
  }
  ## getting a fourier basis
  x <- matrix( sim[,"X"] , nrow = n)
  colnames(x) <- "X"
  
  y <- matrix( sim[,"Y"] , nrow=n )
  colnames(y) <- "Y"
  
  if(isWithError){
    z <- matrix(sim[,"Z"], nrow=n)
    colnames(z)<- "Z"
  }
  
  w <- matrix(sim[,colNamesW], nrow=n)
  colnames(w) <- colNamesW
  
  ## create a fourier basis functions {constant, sin, cos}
  ## maxFourier_basis <<- 30 ## arbirary set of maxFourierBasis which might be useful for 
  if(isComingFromMarginalComputation || computeFourierOfW){
    w_prime <- w
    if(ncol(w)==1){
      w_prime <- as.vector(w)
      maxFourier_basis <- lengthBasis
    }
    
    #fourier_basis <- fourier_extract(w_prime,n)
    fourier_basis <- fourier_extract(w_prime,maxFourier_basis)
  }
  else{
    fourier_basis <- fourier_extract(x,n)
  }
  
  
  ## create a B-spline basis functions - create a n X n smoothing matrix. If argument "eigenVal" is set to true
  ## the basis created will based of orthogonal vectors.
  
  ## temp change
  useFourierInsteadOfSpline <- TRUE
  if(useFourierInsteadOfSpline){
    #spline_basis <- smooth.matrix(w, df=spline_df, eigenValSplineBasis)
    w_prime <- w
    if(ncol(w)==1){
      w_prime <- as.vector(w)
      maxFourier_basis <- lengthBasis
    }
    ##spline_basis <- fourier_extract(w_prime,n)
    spline_basis <- fourier_extract(w_prime,maxFourier_basis)
  }
  else{
    spline_basis <- smooth.matrix(w, df=spline_df, eigenValSplineBasis)
  }
  
  ## storing final results
  result <- list()
  result$fourier <- fourier_basis
  result$spline <- spline_basis
  result$Y <- y
  result$W <- w
  result$X <- x
  if(isWithError){
    result$Z <- z
  }
  
  
  ##
  ## cleaning up local variables
  ##
  rm(sim, fourier_basis, spline_basis)
  
  ## returning variables
  return(result)
}


##
## compute the estimate of Y given the matrix Theta and both basis functions
##

computeEstimateY <- function( theta , fourierBasis, splineBasis){
  
  if(!is.matrix(fourierBasis)){
    fourierBasis <- matrix(fourierBasis, nrow=1,ncol=length(fourierBasis))
  }
  
  if(!is.matrix(splineBasis)){
    splineBasis <- matrix(splineBasis, nrow=1,ncol=length(splineBasis))
  }
  
  n <- nrow(fourierBasis)
  result <- rep(0,n)
  
  values <- 1:n
  fun <- function(i){
    answer <- t(fourierBasis[i,]) %*% theta %*% splineBasis[i,]
    return(answer)
  }
  
  result <- mclapply(values, fun, mc.cores = subFnMcCores)
  
  ## returning final variables
  return(unlist(result))
}


##
## we want here to plot the different elements of our matrix
##
## "theta" contains the different elements of our estimate calculation
##

plot_marginal_variables <- function(theta, plot_type="p"){
  n <- nrow(theta)
  p <- ncol(theta)
  x <- 1:n
  y <- theta[,1]
  plot(x,y, type=plot_type, col = 1, ylim=range(theta))
  
  for( i in 2:p){
    y <- theta[,i]
    lines(x , y , type=plot_type, col=i)
  }
  
  #legend("topleft", legend=paste("k = ",1:p) , pch=1:p, cex=0.8,  col=c(1:p))
}

plot_error <- function(marginal){
  Y <- marginal$basisResult$Y
  X <- marginal$basisResult$X
  Phi_W <- marginal$basisResult$spline
  
  n <- ncol(Phi_W)
  theta <- rep(0,n)
  eta <- rep(0,n)
  
  for(i in 1:n){
    theta[i] <- median(marginal$marginalVariables$theta[,i])
    eta[i] <- median(marginal$marginalVariables$eta[,i])
  }
  
  error <- rep(0, length(Y))
  for(i in 1:n){
    error[i] <- Y[i] - Phi_W[i,] %*% (eta + X[i]*theta)
  }
  
  base <- 1:n
  plot(base, error, type = "p", col="red", ylim=range(error))
}

