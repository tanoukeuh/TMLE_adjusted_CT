##########################################################################################################################
##
##  Cabral CHANANG - February 2017
## 
##########################################################################################################################
##  Description of the file :
##
##  The following file will contain all the necessary funcitons needed in order to produced pre-determine distribution
##  for the variables X,Y and W. We will use this base as a way to define the function "f" of our parameter
##
##########################################################################################################################


##
## we will introduce distributions of the variable "X" such that X \in [0,1]
## by default, we will assume it to be a uniform distribution
##

## n is the number of variable we want to get from the distribution
sample_X <- function(n, distribution_type = "uniform", distribution_typeY = "imply"){
  
  if( n <= 0){
    throw("the variable n should be a strictly positive number - it is currently equal to ", n )
  }
  
  result <- rep(0,n)
  if(distribution_type == "uniform"){
    result <- runif(n)  ## by default, runif has a max set a 1 and a min set at 0
  }
  else if(distribution_type == "beta"){
    ## these are arbitraries values we have choosen for the beta distribution
    ## might be change at a latter stage
    alpha <- 2
    beta <- 2
    result <- rbeta(n, alpha , beta)
  }
  else{
    throw("the distribution type ", distribution_type , " is not recognized.")
  }
  
  ##if(distribution_typeY == "imply"){ ## the 10% element of X will be set to Zero.
    if(TRUE){
      p <- round(n/10,0)
      setZero <- round( runif(p) * 10 , 0 ) +1
      result[setZero] <- 0
    }
  ##}
  
  return(result)
}

##
## distribution of the variable Z = Z + \nu, where \nu is a gaussian distribution of  mean 0 and variance \mu
##
sample_Z  <- function(X, mu){
  ## checking the validity of few variables
  if(mu<0){
    throw("[sample_Z]:The variable mu must be strictly positive")
  }
  
  ## sampling the value of Z
  listOfElementsOfNu <- rnorm(length(X), mean=0, sd=mu)
  Z <- X + listOfElementsOfNu
  
  return(Z)
}


##
## we will introduce distributions of the variable "Y"  - Y \in \R
## by default, we will assume it to be a standard normal distribution
##

## n is the number of variable we want to get from the distribution
sample_Y <- function(n, distribution_type = "sd_normal", X = NULL, W = NULL, isYRandom = TRUE,
                     fn_typeF = "identity", fn_typeG="identity", breakingFactor = 1){
  
  if( n <= 0){
    throw("the variable n should be a strictly positive number - it is currently equal to ", n )
  }
  
  result <- rep(0,n)
  if(distribution_type == "sd_normal"){
    result <- rnorm(n)
  }
  else if(distribution_type == "imply"){
    ## we want to derived the value of Y from the variables W, X , f_W and g_W
    # making sure that X and W are not null
    if(is.null(X)){
      throw("Can not sample Y. The variable X has not been specified.")
    }
    if(is.null(W)){
      throw("Can not sample Y. The variable W has not been specified.")
    }
    
    epsilon <- rep(0,n)
    if(isYRandom){
      epsilon <- rnorm(n)
    }
    
    result <- X * f_W(W, fn_typeF, breakingFactor) + g_W(W, fn_typeG) + epsilon
  }
  else if(distribution_type == "normal"){
    ## these are arbitraries values we have choosen for the normal distribution
    ## might be change at a latter stage
    meanVal <- 4
    sdVal <-  3
    result <- rnorm(n, mean = meanVal, sd=sdVal)
  }
  else if(distribution_type == "gamma"){
    shapeVal <- 4
    result <- rgamma(n, shape=shapeVal)
  }
  else{
    throw("the distribution type ", distribution_type , " is not recognized.")
  }
  
  return(result)
}



##
## we will introduce distributions of the variable "Y"  - Y \in \R
## by default, we will assume it to be a standard normal distribution
##

## n is the number of variable we want to get from the distribution
## d is the number of columns of the matrix W
sample_W <- function(n, d=1, distribution_type = "sd_normal"){
  
  if( n <= 0){
    throw("the variable n should be a strictly positive number - it is currently equal to ", n )
  }
  
  if( d <= 0){
    throw("the variable d should be a strictly positive number - it is currently equal to ", d )
  }
  
  result <- matrix(0, ncol=d, nrow = n) 
  if(distribution_type == "sd_normal"){
    for( i in 1:d){
      result[,i] <- rnorm(n)
    }
  }
  else if(distribution_type == "normal"){
    ## these are arbitraries values we have choosen for the normal distribution
    ## might be change at a latter stage
    meanVal <- 4
    sdVal <-  3
    for(i in 1:d){
      result[,i] <- rnorm(n, mean = meanVal, sd=sdVal)
    }
  }
  else if(distribution_type == "uniform"){
    for(i in 1:d){
      result[,i] <- runif(n)
    }
  }
  else if(distribution_type == "beta"){
    ## these are arbitraries values we have choosen for the beta distribution
    ## might be change at a latter stage
    alpha <- 2
    beta <- 2
    for(i in 1:d){
      result[,i] <- rbeta(n, alpha , beta)
    }
  }
  else if(distribution_type == "gamma"){
    shapeVal <- 4
    for(i in 1:d){
      result[,i] <- rgamma(n, shape=shapeVal)
    }
  }
  else{
    throw("the distribution type ", distribution_type , " is not recognized.")
  }
  
  return(result)
}


##
## defining the function h which is used to imply g
## in fact, we are looking for a function g such that Y = g(X,W) + \epsilon
## where g(X,W) = f(X,W) + E_P[Y | x=0,W] = f(X,W) + h(W)
##

g_W <- function(W, fn_type = "identity"){
  ncolW <- 1
  nrowW <- 1
  if(is.vector(W)){
    nrowW <- length(W)
  }
  else if(is.matrix(W)){
    nrowW <- nrow(W)
    ncolW <- ncol(W)
  }
  else{
    throw("the format of the variable W is not recognizable !")
  }
  
  matrixW <- matrix(W, ncol=ncolW, nrow = nrowW)
  result <- rep(0,nrowW)
  
  if(fn_type == "identity"){
    for(i in 1:nrowW){
      result[i] <- sum(matrixW[i,])
    }  
  }
  else if(fn_type == "square"){
    for(i in 1:nrowW){
      result[i] <- sum(matrixW[i,]^2)
    }
  }
  else if(fn_type == "cube"){
    for(i in 1:nrowW){
      result[i] <- sum(matrixW[i,]^3)
    }
  }
  else if(fn_type == "zero"){
    for(i in 1:nrowW){
      result[i] <- 0
    }
  }
  else if(fn_type == "g_sine"){
    for(i in 1:nrowW){
      result[i] <- sum( 9*cos(2*pi*3*matrixW[i,]) +  5*cos(2*pi*6*matrixW[i,]) + sqrt(3)*sin(2*pi*9*matrixW[i,]) - 7*sin(2*pi*11*matrixW[i,]))
    }
  }
  else if(fn_type == "g_sine_small"){
    for(i in 1:nrowW){
      result[i] <- sum( 4*cos(2*pi*3*matrixW[i,]) + cos(2*pi*6*matrixW[i,]) + sqrt(5)*sin(2*pi*2*matrixW[i,]) - 2*sin(2*pi*10*matrixW[i,]))
    }
  }
  else{
    throw("the fn_type is not recognized")
  }
  
  return(result)
}

##
## defining the function f which will be used in our cases
##

## the function W will be consider to be a matrix
f_W <- function(W, fn_type = "identity", breakPointFactor = 1){
  
  ncolW <- 1
  nrowW <- 1
  if(is.vector(W)){
    nrowW <- length(W)
    matrixW <- matrix(W, ncol = ncolW, nrow = nrowW)
  }
  else if(is.matrix(W)){
    nrowW <- nrow(W)
    ncolW <- ncol(W)
    matrixW <- W
  }
  else{
    throw("the format of the variable W is incorrect !")
  }
  
  result <- rep(0, nrowW)
  
  if(fn_type=="identity"){
      for(i in 1:nrowW){
        result[i] <- sum(matrixW[i,])
      }
  }
  else if(fn_type=="square"){
    for(i in 1:nrowW){
      result[i] <- sum(matrixW[i,]^2)
    }
  }
  else if(fn_type=="sine"){
    for(i in 1:nrowW){
      result[i] <- sum(sin(2*pi*4*matrixW[i,]))
    }
  }
  else if(fn_type=="sine_test"){
    for(i in 1:nrowW){
      # result[i] <- sum(  sin(2*pi*4*matrixW[i,])  +    cos(2*pi*7*matrixW[i,]) +     sin(2*pi*13*matrixW[i,]) +   sin(2*pi*17*matrixW[i,])     )
      # result[i] <- sum( 20*cos(2*pi*5*matrixW[i,]) + 20*sin(2*pi*5*matrixW[i,]) )
      result[i] <- sum( 20*cos(2*pi*5*matrixW[i,]) +  5*cos(2*pi*2*matrixW[i,]))
      #result[i] <- sum( 20*cos(2*pi*matrixW[i,]) +  5*sin(2*pi*2*matrixW[i,]))   ## for TMLE
      
      #result[i] <- sum( 20*cos(2*pi*5*matrixW[i,]) +  10*cos(2*pi*2*matrixW[i,]) + sqrt(5)*sin(2*pi*8*matrixW[i,]) - 7*sin(2*pi*11*matrixW[i,]))
    }
  }
  else if(fn_type=="sine_test_small"){
    for(i in 1:nrowW){
      result[i] <- sum( 20*cos(2*pi*5*matrixW[i,]) +  10*cos(2*pi*2*matrixW[i,])  + 0.05*sin(2*pi*7*matrixW[i,])    )
      # result[i] <- sum( 20*cos(2*pi*5*matrixW[i,]) +  10*cos(2*pi*2*matrixW[i,]) + sqrt(5)*sin(2*pi*8*matrixW[i,]) - 7*sin(2*pi*11*matrixW[i,]))
    }
  }
  else if(fn_type=="breaking_point"){
    for(i in 1:nrowW){
      result[i] <- sum( 20*cos(2*pi*5*matrixW[i,]) +  10*cos(2*pi*2*matrixW[i,])  + breakPointFactor * sin(2*pi*7*matrixW[i,])  - 7*sin(2*pi*11*matrixW[i,]))
    }
  }
  else if(fn_type=="sine_2d"){
    ##for(i in 1:nrowW){
    ##  result[i] <- sum( 20*cos(2*pi*5*matrixW[i,]) +  10*cos(2*pi*2*matrixW[i,])  + breakPointFactor * sin(2*pi*7*matrixW[i,])  - 7*sin(2*pi*11*matrixW[i,]))
    ##}
    
    result <- 10*cos(2*pi*2*matrixW[,1])*cos(2*pi*matrixW[,2]) + 4*sin(2*pi*matrixW[,1])*cos(2*pi*matrixW[,2])
  }
  else if(fn_type=="sine_3d"){
    ##for(i in 1:nrowW){
    ##  result[i] <- sum( 20*cos(2*pi*5*matrixW[i,]) +  10*cos(2*pi*2*matrixW[i,])  + breakPointFactor * sin(2*pi*7*matrixW[i,])  - 7*sin(2*pi*11*matrixW[i,]))
    ##}
    
    result <- 10*cos(2*pi*2*matrixW[,1])*sin(2*pi*matrixW[,2])*cos(2*pi*2*matrixW[,3]) + 4*sin(2*pi*matrixW[,1])*cos(2*pi*matrixW[,2])*cos(2*pi*2*matrixW[,3])
  }
  else{
    throw("the fn_type is not recognized")
  }
  
  return(result)
}


### generate the full scale of variable (W,X,Y) given the size of the sample and the dimension 
### of the covariates variables

## n is the size of the sample
## d is the number of column of the variables W
sample_generator <- function(n , d=1, distributionTypeW="sd_normal", distributionTypeX="uniform",
                             distributionTypeY="sd_normal",
                             typeFunctionF = "identity",
                             typefunctionG = "identity",
                             breakingFactor = 1,
                             isWithError = FALSE,
                             mu=NULL){
  
  ## checking the validity of few variables
  if(isWithError){
    if(is.null(mu)){
      throw("[sample_generator]:The value of the variable mu must be provided.")
    }
  }
  
  result <- matrix(0, ncol = d+2 , nrow=n)
  
  sColnames <- c("Y","X","W")
  if(d>1){
    sColnames <- c(sColnames,paste("W",1:(d-1),sep=""))
  }
  colnames(result) <- sColnames
  
  ## computing the values for X , Y and W
  X <- sample_X(n,distribution_type = distributionTypeX)
  
  ##X <- rep(1,n)
  W <- sample_W(n, d=d, distribution_type = distributionTypeW)
  
  if(distributionTypeY == "imply"){
    Y <- sample_Y(n, distribution_type = distributionTypeY , X = X, W = W, fn_typeF = typeFunctionF , fn_typeG = typefunctionG, breakingFactor = breakingFactor)
  }
  else{
    Y <- sample_Y(n, distribution_type = distributionTypeY, fn_typeF = typeFunctionF , fn_typeG = typefunctionG)  
  }
  
  result[,"X"] <- X
  result[,"Y"] <- Y
  colNamesW <- setdiff(sColnames, c("Y","X"))
  result[,colNamesW] <- W
  
  ## computing the value of Z = X + \nu, where \nu is a gaussian distribution with mean 0 and variable equal to mu^2
  if(isWithError){
    Z <- sample_Z(X, mu)
    result_Z <- matrix(Z, ncol=1,nrow=length(Z))
    colnames(result_Z) <- "Z"
    ## merging both results
    result <- cbind(result, result_Z)
  }
  
  return(result)
}
