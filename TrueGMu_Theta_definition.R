##
## registring the function that will be used to defined trueGMu
##


## Conditional expectation of X given W
#mu <- function(W) {
#  
#  vec <- f(O[, "X"] - O[2, "X"])
#  vec[is.na(vec)] <- 0 ## any real number would do here!
#  res <- condProbUW(W[,"W"]) %*% vec
#  dim(res) <- NULL
#  
#  res <- W * as.vector(res)
#  res
#}


mu <- function(W) {
  
  vec <- f(O[, "X"] - O[2, "X"])
  vec[is.na(vec)] <- 0 ## any real number would do here!
  res <- condProbUW(W[,"W"]) %*% vec
  dim(res) <- NULL
  
  res <- as.vector(res)
  res
}

muAux <- function(W) {
  mu(W)/(1-g(W))
}


##
## definition of condProbUW
##


condProbUW <- function(W) {
  u <- 1
  cpu1 <- rep(0, length(W))
  if (p[u]!=0) {
    cpu1 <- p[u] *
      dnorm(logit(W), mean=logit(O[u, "W"]), sd=omega[u])/(W*(1-W)) / omega[u]
  }
  
  u <- 2
  cpu2 <- rep(0, length(W))
  if (p[u]!=0) {
    cpu2 <- p[u] *
      dnorm(logit(W), mean=logit(O[u, "W"]), sd=omega[u])/(W*(1-W)) / omega[u]
  }
  
  u <- 3
  cpu3 <- rep(0, length(W))
  if (p[u]!=0) {
    cpu3 <- p[u] *
      dnorm(logit(W), mean=logit(O[u, "W"]), sd=omega[u])/(W*(1-W)) / omega[u]
  }
  out <- cbind(cpu1, cpu2, cpu3)
  ## make sure that they sum to 1
  Z <- apply(out, 1, sum, na.rm = TRUE)
  out <- out/Z
  out
}



##
##  definition of g
##
g <- function(W) {
  condProbUW(W[,"W"])[, 2]
}

trueGMu <- list(  g =g , mu = mu, muAux = muAux )

##
## registring function that will be used to define trueTheta
##

## Conditional expectation of Y given (W,X)
theta <- function(XW) {   
  #browser()
  XW_init <- XW
  XW <- cbind(XW[,c("X","W")])
  XW[, "X"] <- XW[, "X"] + O[2,"X"]
  X <- XW[, 1]        
  W <- XW[, "W"]
  
  dummy <- 1
  param <- c(Sigma1[1,2]/Sigma1[1,1],
             dummy,
             Sigma3[1,2]/Sigma3[1,1])
  ## if (X_1,X_2)~N(mu, S) then E(X_2|X_1)=mu_2+S_12/S_11*(X_1-mu_1)
  cpu <- condProbUXW(XW)
  
  res <- 0
  res1 <- (O[1, "Y"] + param[1] * (X-O[1, "X"])) * cpu[, 1]
  if (p[1] != 0) {
    res <- res + res1
  }
  
  res2 <- (O[2, "Y"] + lambda0(W)) * cpu[, 2]
  if (p[2] != 0) {
    res <- res+res2
  }
  
  res3 <- (O[3, "Y"] + param[3] * (X-O[3, "X"])) * cpu[, 3]
  if (p[3] != 0) {
    res <- res + res3
  }
  
  dim(res) <- NULL
  res
}

## function 'theta0'
theta0 <- function(W) {    
  #browser()
  XW <- cbind(X=0, W);
  theta(XW);
}

condProbUXW <- function(XW) {
  #X <- XW[, 1]
  #W <- XW[, 2]
  
  X <- XW[, "X"]
  W <- XW[, "W"]
  
  u <- 1
  cpu1 <- rep(0, nrow(XW))
  if (p[u]!=0) {
    cpu1 <- p[u] *
      dnorm(logit(W), mean=logit(O[u, "W"]), sd=omega[u])/(W*(1-W)) *
      dnorm(X, mean=O[u, "X"], sd=Sigma1[1, 1]) / omega[u]
  }
  
  u <- 2
  cpu2 <- rep(0, nrow(XW))
  if (p[u]!=0) {
    cpu2 <- p[u] *
      dnorm(logit(W), mean=logit(O[u, "W"]), sd=omega[u])/(W*(1-W)) *
      (X==O[u, "X"]) / omega[u]
  }
  
  u <- 3
  cpu3 <- rep(0, nrow(XW))
  if (p[u]!=0) {
    cpu3 <- p[u] *
      dnorm(logit(W), mean=logit(O[u, "W"]), sd=omega[u])/(W*(1-W)) *
      dnorm(X, mean=O[u, "X"], sd=Sigma3[1, 1]) / omega[u]
  }
  out <- cbind(cpu1, cpu2, cpu3)
  ## make sure that they sum to 1
  Z <- apply(out, 1, sum, na.rm = TRUE)
  out <- out/Z
  out
}

trueTheta <- list(theta = theta , theta0 = theta0)
