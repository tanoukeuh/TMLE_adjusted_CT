#############################################################################################################
##
##  Cabral CHANANG - December 2016
## 
#############################################################################################################
##  Description of the file :
##
##  The following file will contain some utilities functions.
##
###############################################################################################################
##
## The computation of the DR basis comes from the following URL :
## https://publish.illinois.edu/liangf/teaching/stat-424/rcode/rcode-smoothing-splines/
##
#############################################################################################################


## Here is how we obtain the Demmler & Reinsch (DR) Basis: we first obtain the smoother matrix S 
## (which is not returned by R, so we write our own script to compute it), and then the eigen-vectors of S are basically the DR basis functions.
smooth.matrix <- function(x, df, eigenVal = TRUE){
  ## return the smoother matrix with knots x and df
  ## this function is for x having unique values
  n = length(x);
  A = matrix(0, n, n);
  for(i in 1:n){
    y = rep(0, n); y[i]=1;
    ###yi = smooth.spline(x, y, df=df)$y;  ## to be used when we know for sure that the values of "x" are unique
    yi = predict(smooth.spline(x, y, df=df),x)$y;
    A[,i]= yi;
  }
  
  ## returning smoothing matrix
  ## either an orthogonal spline basis or a regular one
  if(eigenVal){
    eigenA <- eigen(A)
    return(eigenA$ve)
  }
  else{
    return(A)
  }
}


##
## based on the true fourier functions but give the option to include or not the constant based
##
fourier_extract <- function(x,n){
  
  ## defined m = n+1, because the base function "fourier" truncated the results. We ensure to have to exact number of columns
  ## once we are out of this function
  m <- n + 1
  ## need to check the dimension of x before moving on.
  if(is.vector(x)){
    result <- fourier(x , nbasis = m, period = 1)
    ## the fourier formula divide all the column by sqrt(2)... we then want to readjust the sin/cos base columns
    result[,2:m] <-  result[,2:m]/sqrt(2)
    result <- result[,1:n]
  }
  else if(is.matrix(x)){
    nColumns <- ncol(x)
    result <- fourier(x[,1], nbasis = m, period = 1)
    ## the fourier formula divide all the column by sqrt(2)... we then want to readjust the sin/cos base columns
    result[,2:m] <-  result[,2:m]/sqrt(2)
    
    for(i in 2:nColumns){
      interm_result <- fourier(x[,i], nbasis = m, period = 1)
      interm_result[,2:m] <-  interm_result[,2:m]/sqrt(2)
      result <- mergeBothMatrices(result, interm_result)
    }
  }
  
  ## return final result
  return(result)
}

##
## generate the legendre basis
## x : data points which can be either a vector o a matrix
## n : number of degree of libery
##
legendre_extract <- function(x,n){
  if(n <= 0){
    throw("[legendre_extract]: the degree of freedom n should be greater than 0")
  }
  
  
  if(is.vector(x)){
    result <- poly(x,n)
  }
  else{
    ## what we need here is the columns of the poly stack together and then we will add a combination of pairwise columns
    nbCols <- ncol(x)
    nbRows <- nrow(x)
    xColNames <- colnames(x)
    
    values <- 1:nbCols
    fun_poly <- function(i){
      midResult <- poly(x[,i],n, simple=TRUE)
      return(midResult)
    }
    
    fun_names <- function(i){
      if(n<2){
        midResult_name <- xColNames[i] 
      }
      else{
        midResult_name <- c( xColNames[i] , paste(xColNames[i],"^",2:n,sep=""))  
      }
      
      return(midResult_name)
    }
    
    ## this provide me with the base case
    poly_list <- mclapply(values, fun_poly, mc.cores = 4)
    names_list <- mclapply(values, fun_names, mc.cores = 4)
    
    poly_matrix <- matrix(unlist(poly_list), nrow=nbRows, ncol=nbCols*n)
    colnames(poly_matrix) <- unlist(names_list)
    
    
    if( n > 1){
      ## now let us construct the pair wise multiplication
      ## in order to echeive that we need to extract the first column of each base case mentioned above and then construct the pairwise multiplication
      extract_first_column <- function(i){
        return(poly_list[[i]][1])  ## it correspond to the base case - i.e the first coordinate of legendre basis
      }
      extract_first_column_name <- function(i){
        return(names_list[[i]][1])
      }
      extract_first_list <- mclapply(values, extract_first_column, mc.cores=4)
      extract_first_list_name <- mclapply(values, extract_first_column_name, mc.cores =4)
      
      extract_first_column_matrix <- matrix(unlist(extract_first_list), ncol=nbCols, nrow=nbRows)
      colnames(extract_first_column_matrix) <- unlist(extract_first_list_name)
      
      ## start building the pairwise multiplication
      build_pairwise_multi <- function(i){
        result <- extract_first_column_matrix[,i] * extract_first_column_matrix[,(i+1):nbCols]
        return(result)
      }
      
      build_pairwise_multi_name <- function(i){
        xColNames <- colnames(extract_first_column_matrix)
        result <- paste( xColNames[i] , xColNames[(i+1):nbCols], sep="*")
      }
      pairWiseValues <- 1:(nbCols-1)
      pairwise_matrix_list <- mclapply(pairWiseValues, build_pairwise_multi, mc.cores = 4)
      pairwise_matrix_list_names <- mclapply(pairWiseValues, build_pairwise_multi_name, mc.cores = 4)
      
      pairwise_matrix <- matrix(unlist(pairwise_matrix_list), nrow=nbRows, ncol=nbCols*(nbCols -1)/2)
      colnames(pairwise_matrix) <- unlist(pairwise_matrix_list_names)
    }
    
    
    ## creating the final matrix which is a merge between the one containing all the legendre basis and the one containing the cross of second order
    const_matrix <- matrix(1, nrow=nbRows, ncol=1)
    colnames(const_matrix) <- "cste"
    
    if(n>1){
      result <- cbind(const_matrix, poly_matrix, pairwise_matrix)
    }
    else{
      result <- cbind(const_matrix, poly_matrix)
    }
  }
  
  return(result)
}

##
## merged two matrices through the inner multiplication of both columns
##
mergeBothMatrices <- function(result, interm_result){
  ## we will assume here that both matrices have the same number of row but potentiallly different number of columns
  
  nbCol_res <- ncol(result)
  nbCol_int <- ncol(interm_result)
  nbRow_res <- nrow(result)
  #newMatrix <- matrix(0, ncol=nbCol_res*nbCol_int, nrow=nbRow_res)
  #colnames(newMatrix) <- rep(0, ncol(newMatrix))
  
  ## transform this into a mclapply functions for speed reason... can it be done ?
  #for( i in 1:nbCol_int){
  #  startCol <- (i-1) * nbCol_res + 1
  #  endCol <- i * nbCol_res
  #  colnames(newMatrix)[startCol:endCol] <- paste(colnames(result), colnames(interm_result)[i], sep="*")
  #  newMatrix[,startCol:endCol] <- result * interm_result[,i]
  #}
  
  fun_merge_matrix <- function(i){
    startCol <- (i-1) * nbCol_res + 1
    endCol <- i * nbCol_res
    fun_result <- result * interm_result[,i]
    colnames(fun_result) <- paste(colnames(result), colnames(interm_result)[i], sep="*")
    return(fun_result)
  }
  
  fun_merge_matrix_names <- function(i){
    colNames <- paste(colnames(result), colnames(interm_result)[i], sep="*")
    return(colNames)
  }
  
  values <- 1:nbCol_int
  fun_answer <- mclapply(values, fun_merge_matrix, mc.cores = 4)
  fun_answer_names <- mclapply(values, fun_merge_matrix_names, mc.cores = 4)
  
  newMatrix <- matrix(unlist(fun_answer), nrow = nbRow_res, ncol=nbCol_res*nbCol_int )
  colnames(newMatrix) <- unlist(fun_answer_names)
  
  ## returning final value
  return(newMatrix)
}

##
## retreive fourier basis for unidimensional fourier 
##
fourier_extract_recent_old <- function(x,n){
  ## need to check the dimension of x before moving on.
  result <- fourier(x , nbasis = n+1, period = 1)
  
  ## the fourier formula divide all the column by sqrt(2)... we then want to readjust the sin/cos base columns
  result[,2:(n+1)] <-  result[,2:(n+1)]/sqrt(2)
  
  ## return the n first columns
  result <- result[,1:n]
  return(result)
}

##
## derivative fourier extract  - to be changed since it does not cover the case where n is odd.
##

## we will assume here that the first column of the fourier basis is a constant
derivative_fourier_extract <- function(x , k){
  fourierBasis <- fourier_extract(x,k+1)
  if(!is.matrix(fourierBasis)){
    fourierNames <- names(fourierBasis)
    fourierBasis <- matrix(fourierBasis, nrow = 1, ncol=length(fourierBasis))
    colnames(fourierBasis) <- fourierNames
  }
  p <- ncol(fourierBasis)-1
  n <- nrow(fourierBasis)
  result <- matrix(0,ncol=p,nrow=n)
  
  values <- 1:p
  fun <- function(i){
    if(i==1){
      answer <- rep(0,n)
    }
    else{
      base <- round(i/2,0)
      if(i%%2 == 0){
        answer <- fourierBasis[,i+1] * base * 2 * pi 
      }
      else{
        answer <- fourierBasis[,i-1] * base * 2 * pi * (-1)
      }
    }
    return(answer)
  }
  fun_answer <- mclapply(values, fun, mc.cores = subFnMcCores)
  result <- matrix(unlist(fun_answer), ncol = p)
  
  ##
  ## removing all local variables
  ##
  rm(fn_answer, fourierBasis)
  
  ## returning final answer
  return(result)
}

fourier_extract_old <- function(x,n,include_constant=FALSE){
  
  base <- fourier(x , nbasis = 2*n+1, period = 1)
  nb_elts <- length(x) * n
  result <- matrix(rep(0,nb_elts), length(x) , n)
  start <- 1
  j <- 1
  if(!include_constant){
    start <- 2
    result[,1] <- base[,1]
    j <- 2
  }
  
  for(i in start:n){
    result[,i] <- base[,j]
    if(j%%2 == 0){
      j <- j + 1
    }
    else{
      j <- j + 3
    }
  }
  
  ##
  ## removing local variables 
  ##
  rm(base)
  
  ## returning final values
  return(result)
}



##
## computing the matrix \Theta, outlined in our second set of idea
##

## Y represents our effect variable
## phi_X represents the basis function representation of function based on X
## Psi_W represents the basis function representation of function based on W 
## iterX represents the size of Phi_X that we wish to use for the estimation exercise
## iterW represents the size of Psi_X that we wish to use for the estimation

## the equation linking all the variables is : Y = Phi_X * \Theta * Psi_W (1)
least_square_solution_Theta <- function(Y, Phi_X, Psi_W, iterX, iterW){
  ## checking that iterX and iterY have the right size
  if(ncol(Phi_X) < iterX){
    throw("the number of X basis to use must be smaller than ", ncol(Phi_X))
  }
  
  if(ncol(Psi_W) < iterW){
    throw("the number of W basis to use must be smaller than ", ncol(Psi_W))
  }
  
  ## the variable \theta will be of size iterX * iterW
  usedPhiX <- Phi_X[,1:iterX]
  usedPsiW <- Psi_W[,1:iterW]
  theta <- rep(0, iterX*iterW)
  
  ## the equation (1) above can be rewritten as Y = X * \theta_v (2) where \theta_v is the matrix \theta transformed into a vector
  ## hence let's construct the matrix X
  X <- matrix(0, nrow = length(Y), ncol = length(theta))
  
  values <- 1:length(Y)
  fun <- function(i){
    answer <- lineX_i(usedPhiX[i,], usedPsiW[i,])
  }
  
  fun_result <- mclapply(values, fun, mc.cores = subFnMcCores)
  X <- matrix(unlist(fun_result), ncol = length(theta), byrow = TRUE)
  
  if(TRUE){
    #lasso.fit <- cv.glmnet(X,Y, alpha = 1, intercept = TRUE, parallel = FALSE)
    lasso.fit <- cv.glmnet(X,Y, alpha = 1, intercept = FALSE, parallel = FALSE)
    lasso_coef <- coef(lasso.fit)
    theta_lasso <- lasso_coef[2:length(lasso_coef)]
    #browser()
    theta_mat <- matrix(theta_lasso, ncol=iterW, nrow=iterX, byrow=TRUE)
    
    ## removing some unused local variables
    rm(lasso.fit, lasso_coef)
  }
  else{
    ## transforming \theta into in a matrix
    theta_mat <- matrix(theta, ncol=iterW, nrow=iterX)
  }
  
  ##
  ## cleaning up local variables
  ##
  rm(fun_result, X)
  
  ## returning the final values
  return(theta_mat)
}


##
## computing the matrix \Theta, outlined in our second set of idea
##

## Y represents the response variable
## X represents the effect variable
## psi_W represents the basis function representation of function E[Y|X=0,W]
## Psi_W represents the basis function representation of function f_W 
## iterW1 represents the size of Psi_W1 that we wish to use for the estimation exercise
## iterW2 represents the size of Psi_W2 that we wish to use for the estimation

## the equation linking all the variables is : Y = [Psi_W1, X*PsiW2]%*%\beta
least_square_solution_beta <- function(Y, X, Psi_W, iterW1, iterW2, isLasso = TRUE, isLassoPenalize = FALSE){
  ## checking that iterX and iterY have the right size
  if(ncol(Psi_W) < iterW1){
    throw("the number of W1 basis to use must be smaller than ", ncol(Psi_W))
  }
  
  if(ncol(Psi_W) < iterW2){
    throw("the number of W2 basis to use must be smaller than ", ncol(Psi_W))
  }
  
  ## the variable \theta will be of size iterX * iterW
  usedPhiW1 <- Psi_W[,1:iterW1]
  usedPsiW2 <- Psi_W[,1:iterW2]
  usedPsiW2X <- X * usedPsiW2
  mat <- cbind(usedPhiW1, usedPsiW2X)
  
  ## removing the first column from the matrix
  mat <- mat[,2:ncol(mat)]
  ## beta <- rep(0, iterW1 + iterW2)
  
  if(isLasso){
    if(!isLassoPenalize){
      lasso.fit <- cv.glmnet(mat,Y, alpha = 1, intercept = TRUE, parallel = FALSE)
    }
    else{
      penalty <- rep(1,iterW1 + iterW2)
      ## only setting the function f variales to be set to 1
      penalty[1:iterW1] <- 0
      lasso.fit <- cv.glmnet(mat,Y, alpha = 1, intercept = TRUE, parallel = FALSE, penalty.factor = penalty)
    }
    lasso_coef <- coef(lasso.fit)
    ## theta_lasso <- lasso_coef[2:length(lasso_coef)]
    theta_lasso <- lasso_coef
    
    ## removing local variables
    rm(lasso.fit)
  }
  else{
    dantzig.fit <- slim(mat, Y, nlambda = 10,  method="dantzig", lambda.min.ratio=.01)
    theta_lasso <- rep(0,iterW1+iterW2)
    dantzigColLength <- ncol(dantzig.fit$beta)
    theta_lasso[0] <- dantzig.fit$intercept[,dantzigColLength]
    theta_lasso[1:(iterW1+iterW2-1)] <- dantzig.fit$beta[,dantzigColLength]
  }
  
  ## removing local variables
  rm(mat)
  
  ## final result
  return(theta_lasso)
}

##
## creating line of the matrix X
##

lineX_i <- function(PhiX_i , PsiW_i){
  result <- rep(0, length(PhiX_i)*length(PsiW_i))
  startIncrement <-1
  endIncrement <- length(PsiW_i)
  
  for(i in 1:length(PhiX_i)){
    result[startIncrement:endIncrement] <- PhiX_i[i] * PsiW_i
    startIncrement <- i * length(PsiW_i) + 1
    endIncrement <- startIncrement + length(PsiW_i) - 1
  }

  return(result)
}


##
##  compute the inverse "INV" of a diagonal vector "D" such that INV_i = 1/D_i.
##  if abs(D_i) < epsilon, then INV_i = D_i.
##  this function will be useful for pseudo inverse calculation
##

inverse_diag <- function(diag, epsilon = 10^-10){
  n <- length(diag)
  result <- rep(0,n)
  
  values <- 1:n
  fun <- function(i){
    if(abs(diag[i]) < epsilon){
      answer = diag[i]
    }
    else{
      answer = 1/diag[i]
    }
    return(answer)
  }
  
  result <- mclapply(values, fun, mc.cores = subFnMcCores)

  return(unlist(result))
}


##
## compute the pseudo inverse of a matrix throught the Singular Value Decomposition
##
##

pseudo_inverse <- function(mat){
  mat_svd <- svd(mat)
  inv_diag <- inverse_diag(mat_svd$d)
  result <- mat_svd$v %*% diag(inv_diag) %*% t(mat_svd$u)
  
  ## removing local variables
  rm(inv_diag, mat_svd)
  
  # final result
  return(result)
}



##
##  plotting the result of the prediction risk 
##

riskPrediction_Result_Plotter <- function( result, step, lengthBaseX, lengthBaseW, method, nbCases, nbRepetitions,
                                            typeFunctionF, typeFunctionG, distribW, nbW, isSparsity = FALSE, cf = 0.95, isFunctionG=FALSE){
  
  ##
  ## creating the graph title and name
  ##
  sRunDescription = paste("method=", method ,"_nbRepetitions=",nbRepetitions,"_Step=",step, "_lengthBaseX=", lengthBaseX,"_lengthBaseW=",
                            lengthBaseW, "_typeFunctionF=",typeFunctionF,"_typeFunctionG=",typeFunctionG,
                          "_distribW=",distribW,  "_nbW=",nbW,  "_conf.Inter=", cf, sep="")
  
  ## 
  ## Plotting the result of the method
  ##
  
  if(isSparsity){
    pathFileName <- paste( "Evolution of Sparsity_" , sRunDescription, "_file.pdf", sep="")
  }else{
    if(isFunctionG){
      pathFileName <- paste( "Evolution Prediction Error Function G_" , sRunDescription, "_file.pdf", sep="")
    }
    else{
      pathFileName <- paste( "Evolution Prediction Error_" , sRunDescription, "_file.pdf", sep="")
    }
  }
  pdf(file=pathFileName)
    
  
  ##  creating error bar 
  n = nrow(result) 
  x <- step *(1:n)
  q_cf <- qnorm( cf + (1-cf)/2)
  y <- result[,1]
  sd <- result[,2]
  shift <- q_cf * ( sd /sqrt(x) )
  
  y_up <- y + shift
  y_down <- y - shift
  
  y_lim <- c(min(y,y_up,y_down), max(y,y_up,y_down))
  x_lim <- c(min(x), max(x))
  #plot(x, y, type="o",  lty = 1 , pch=1, col=1, ylab="Values" , xlab = "Sample size", ylim= y_lim, xlim= x_lim)
  #error.bar(x, y , y_up, y_down)
  errbar(x , y, y_down  , y_up , type="o", lty = 1 , pch=1, col=1, ylab="Values" , xlab = "Sample size",
           ylim=y_lim)    
  #plot(x, y , ylim=y_lim, xlim=x_lim)
  #errbar(x , y, y_down  , y_up , type="o", lty = 1 , pch=1, col=1, ylab="Values" , xlab = "Sample size", asp=2)  
  
  if(isSparsity){
    title(main= paste("Evol. of Sparsity using method ", method," and ", nbRepetitions, " Repetitions. per step")) 
  }
  else{
    if(isFunctionG){
      title(main= paste("Evol. of Risk Pred. Fn G using method ", method," and ", nbRepetitions, " Reptitions. per step"))
    }
    else{
      title(main= paste("Evol. of Risk Pred. using method ", method," and ", nbRepetitions, " Reptitions. per step"))  
    }
    
  }
  dev.off()
}


error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  #arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ylim=ylim,...)
  arrows(x,upper, x, lower, angle=90, code=3, length=length,...)
}





##
## capacity to load a RData file and then plot key figures
## Here is the list of figures :
## 1. Evolution of sparsity
## 2. Evolution of prediction Risk
## 3. Evolution of estimated fourier coefficients (only available for method 2)
##


plotKeyFiguresOfSimulation <- function(fileName, plotTMLE = FALSE){
  #loading the main file
  load(paste(fileName,".RData",sep=""))
  
  # ploating the key figures of our simulation 
  ## adding element of the run to overall result
  nbRepetitions <- mainResult$nbRepetitions 
  nbCases <- mainResult$nbCase
  step <- mainResult$step 
  lengthBaseX <- mainResult$lengthBaseX
  lengthBaseW <- mainResult$lengthBaseW
  method <- mainResult$method 
  typeFunctionF <- mainResult$typeFunctionF
  typeFunctionG <- mainResult$typeFunctionG
  distribW <- mainResult$distribW
  nbW <- mainResult$nbW 
  isDataForErrorG_NotPresent<- is.null(mainResult$errorG)
  
  ## TMLE
  if(is.null(mainResult$run_tmle_algo)){
    run_algo_TMLE <- FALSE
  }
  else{
    run_algo_TMLE <- mainResult$run_tmle_algo
  }
  
  if(TRUE){
    # plotting the evolution of error
    error <- mainResult$error
    if(run_algo_TMLE){
      error_tmle <- mainResult$error_tmle
      riskPrediction_comparison_Result_Plotter(error, error_tmle, step, lengthBaseX, lengthBaseW, method, nbCases, nbIterations, typeFunctionF, typeFunctionG, distribW, nbW);
    }
    else{
      riskPrediction_Result_Plotter(error, step, lengthBaseX, lengthBaseW, method, nbCases, nbRepetitions, typeFunctionF, typeFunctionG, distribW, nbW);
    }
    
    # plotting the evolution of sparsity
    sparsity <- mainResult$sparsity
    riskPrediction_Result_Plotter(sparsity, step, lengthBaseX, lengthBaseW, method, nbCases, nbRepetitions, typeFunctionF, typeFunctionG, distribW, nbW, TRUE);
    
    ## plotting the evolution of error of function g
    if(!isDataForErrorG_NotPresent){
      errorG <- mainResult$errorG
      riskPrediction_Result_Plotter(error, step, lengthBaseX, lengthBaseW, method, nbCases, nbRepetitions, typeFunctionF, typeFunctionG, distribW, nbW, isFunctionG=TRUE);
    }
  }
  
  # plotting the evolution of the estimate of f_w
  matrixBeta <- matrix(0, ncol=nbCases, nrow=lengthBaseW)
  for(i in 1:nbCases){
    caseType <- paste("case-",i,sep="")
    iterBeta <- rep(0, lengthBaseW)
    for(j in 1:nbRepetitions){
      if(plotTMLE && run_algo_TMLE){
        iterBeta <- colMeans( rbind(colMeans(mainResult[[caseType]][[j]]$result_matrix_tmle) , iterBeta ))
      }
      else{
        iterBeta <- colMeans( rbind(colMeans(mainResult[[caseType]][[j]]$estimate_theta) , iterBeta ))
      }
    }
    matrixBeta[,i] <- iterBeta
  }
  riskEstimate <- list()
  riskEstimate$fw <- estimateFW_plotter(method, typeFunctionF, typeFunctionG, nbW, nbCases, lengthBaseX, lengthBaseW, matrixBeta)
  
  
  if(!isDataForErrorG_NotPresent){
    # plotting the evolution of the estimate of g_w
    matrixEta <- matrix(0,ncol=nbCases, nrow=lengthBaseX)
    for(i in 1:nbCases){
      caseType <- paste("case-",i,sep="")
      iterEta <- rep(0, lengthBaseX)
      for(j in 1:nbRepetitions){
        iterEta <- colMeans( rbind(colMeans(mainResult[[caseType]][[j]]$estimate_eta) , iterEta ))
      }
      matrixEta[,i] <- iterEta
    }
    riskEstimate$gw <- estimateFW_plotter(method, typeFunctionF, typeFunctionG, nbW, nbCases, lengthBaseX, lengthBaseW, matrixEta, isFunctionG=TRUE)
  }
  
  return(riskEstimate)
}

##
## plotting the evolution of the estimate of function f, relatively to its true values....
##
estimateFW_plotter <- function(method, typeFunctionF, typeFunctionG, nbW, nbCases, lengthBaseX, lengthBaseW, matrixBeta , isFunctionG=FALSE){
  
  # can only plot estimate if method is second one
  if(method!=2){
    throw("Can not yet plot the evolution of the estimate of function f_w for method different than method 2 !")
  }

  ## creating the graph title and name
  
  sRunDescription = paste("method=", method , "_lengthBaseX", lengthBaseX , "_lengthBaseW=", lengthBaseW, "_typeFunctionF=",typeFunctionF, "_typeFunctionG=", typeFunctionG , "_nbW=",nbW, sep="")

  ## Plotting the result of the method
  if(isFunctionG){
    pathFileName <- paste( "Evolution of estimate of g_W with_" , sRunDescription, "_file.pdf", sep="")
  }
  else{
    pathFileName <- paste( "Evolution of estimate of f_W with_" , sRunDescription, "_file.pdf", sep="")  
  }
  pdf(file=pathFileName)
  
  # values of W to be plotted - has to be between [0.1]
  sequence <- 0.005
  base_element <- seq(0,1,sequence)
  nBaseElement <- length(base_element)
  if(nbW > 1){
    list_base_element <- list()
    for(i in 1:nbW){
      list_base_element[[i]] <- base_element
    }
    base_element <- expand.grid(list_base_element)
    base_element_2d <- data.matrix(base_element)   # casting the variable into a matrix
    base_element <- seq(0,1,sequence)
  }
  
  # plotting the true function f_W/g_W
  if(isFunctionG){
    if(nbW==1){
      base_result <- g_W(base_element, typeFunctionG)
    }
    else if(nbW==2){
      base_result <- g_W(base_element_2d, typeFunctionG)
    }
    ylab_text <- "g(W) - values"
    if(nbW == 2){
      zlab_text <- ylab_text
    }
  }
  else{
    if(nbW==1){
      base_result <- f_W(base_element, typeFunctionF)
    }
    else if(nbW==2){
      base_result <- f_W(base_element_2d, typeFunctionF)
    }
    ylab_text <- "f(W) - values"
    if(nbW == 2){
      zlab_text <- ylab_text
    }
  }
  
  if(nbW == 1){
    plot(base_element, base_result, type="l", col=1, xlab="W - elements", ylab=ylab_text)
  }
  else if(nbW ==2){
    #####plot(base_element, base_result, type="l", col=1, xlab="W - elements", ylab=ylab_text)
    #####persp(base_element, base_element, base_result, type="l", col=1, xlab="W1 - elements", ylab="W2 - elements", zlab=zlab_text)
    persp(base_element, base_element, matrix(base_result, ncol=length(base_element)), col=1, xlab="W1 - elements", ylab="W2 - elements", zlab=zlab_text)
  }
  
  # adding all estimates of f_W/g_W
  if(isFunctionG){
    if(nbW==1){
      base_fourier_extract <- fourier_extract(base_element, lengthBaseX)
    }
    else if(nbW==2){
      base_fourier_extract <- fourier_extract(base_element_2d, maxFourier_basis)
      base_fourier_extract <- base_fourier_extract[,(1:lengthBaseX)]
    }
  }
  else{
    if(nbW==1){
      base_fourier_extract <- fourier_extract(base_element, lengthBaseW)  
    }
    else if(nbW==2){
      base_fourier_extract <- fourier_extract(base_element_2d, maxFourier_basis)
      base_fourier_extract <- base_fourier_extract[,(1:lengthBaseW)]
    }
  }
  estimate_result <- rep(0,nbCases)
  
  for(i in 1:nbCases){
    #print(paste("i am here : " , i))
    if(nbW==1){
      estimate_fW_element <- base_fourier_extract %*% matrixBeta[,i]
      lines(base_element, estimate_fW_element, type="l", col=i+1)
      estimate_result[i] <- sum((base_result - estimate_fW_element)^2)*sequence
    }
    else{
      estimate_fW_element <- base_fourier_extract %*% matrixBeta[,i]
      if(nbW==2){
        par(new=TRUE)
        persp(base_element, base_element, matrix(estimate_fW_element, ncol=length(base_element)), col=i+1, xlab="W1 - elements", ylab="W2 - elements", zlab=zlab_text)
      }
      estimate_result[i] <- sum((base_result - estimate_fW_element)^2)*sequence*sequence
    }
  }
  
  if(nbW <=2 ){
    legend("topleft", c("true function", paste(1:nbCases," Iteration", sep="")), col=c(1:(nbCases+1)), cex=.5,lwd=3)
    
    if(isFunctionG){
      title(main= paste("Evol. estimate of function g_W"))
    }
    else{
      title(main= paste("Evol. estimate of function f_W"))  
    }
  }
   
  dev.off()
  return(estimate_result)
}





##
##  plotting the comparison of risk prediction between "lasso" and "tmle" 
##
riskPrediction_comparison_Result_Plotter <- function( result_lasso, result_tmle, step, lengthBaseX, lengthBaseW, method, nbCases, nbRepetitions,
                                           typeFunctionF, typeFunctionG, distribW, nbW, cf = 0.95){
  
  ##
  ## creating the graph title and name
  ##
  sRunDescription = paste("comparison_between_lasso_tmle_" ,"_nbRepetitions=",nbRepetitions,"_Step=",step, "_lengthBaseX=", lengthBaseX,"_lengthBaseW=",
                          lengthBaseW, "_typeFunctionF=",typeFunctionF,"_typeFunctionG=",typeFunctionG,
                          "_distribW=",distribW,  "_nbW=",nbW,  "_conf.Inter=", cf, sep="")
  
  ## 
  ## Plotting the result of the method
  ##
  pathFileName <- paste("Comparison_Evolution_Error_with_TMLE_" , sRunDescription, "_file.pdf", sep="")
  pdf(file=pathFileName)
  
  ##  creating error bar 
  n = nrow(result_lasso) 
  x <- step *(1:n)

  y_lasso <- result_lasso[,1]
  y_tmle <- result_tmle[,1]
  
  y_lim <- range(y_lasso, y_tmle)
  x_lim <- range(x)
  
  
  ## lasso results
  ##plot(x, y_lasso, type="o",  lty = 1 , pch=1, col=1, ylab="Values" , xlab = "Sample size", ylim= y_lim, xlim= x_lim)
  ## tmle results
  ##plot(x, y_tmle, type="o",  lty = 1 , pch=1, col=2)
  
  colNames <- c("LASSO type","TMLE")
  matY <- cbind(y_lasso, y_tmle)  
  matplot( x ,  y = matY , type="o", lty = 1:2 , pch = 1:2 , col=c(1:2) , ylab="Values" , xlab = "Sample size" , ylim= y_lim, xlim= x_lim )
  legend("topleft", legend=colNames, lty = 1:2 , pch = 1:2 , cex = 0.8 , col=c(1:2) )
  
  title(main=paste("Evol. of Risk Prediction , using", nbRepetitions, " Repetitions per sample size"))
  dev.off()
}