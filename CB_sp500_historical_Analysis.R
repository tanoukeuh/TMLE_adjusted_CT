#############################################################################################################
##
##  Cabral CHANANG - December 2017
## 
#############################################################################################################
##  Description of the file :
##
##  Implementation of historical analysis based on the S&P 500
##
##
#############################################################################################################

## loading library
library(lubridate)


##
## calculate the presence of each GICS sectors in the final portfolio relative to the original one
##
presenceInGICSSectors <- function(coefResult){
  
  constituentsFullDescription <- loadObject("sp500_data_list.RData")
  listGicsSector <- unique(constituentsFullDescription[,"GICS Sector"])
  result <- matrix(1, nrow=1,ncol=length(listGicsSector))
  
  colnames(result) <- listGicsSector
  
  for(i in 1:length(listGicsSector)){
    result[i] <- length(which( coefResult$descriptionLightCoef[,"GICS Sector"] == listGicsSector[i]))/length(which( constituentsFullDescription[,"GICS Sector"] == listGicsSector[i]))
  }
  
  print(result)
}




##
## plotting the daily return corresponding to the SP500, as well as the difference with the replicating portfolio
##

plotDifferenceWithReplicatingPortfolio <- function(coefResult, yearProjection){
  dailyReturn <- coefResult$returnY
  dailyDifference <- coefResult$differenceReturn
  datesProjection <- as.Date(coefResult$projectionDates[1:length(dailyReturn)])
  
  
  pathFileName <- paste("Comparison_Evolution_DailyReturn_" , as.character(yearProjection), "_file.pdf", sep="")
  pdf(file=pathFileName)
  
  ##  creating the graph
  graphNames <- c("SP500 Daily return","Diff. of daily return betw. SP500 and replica. portf.")
  matY <- cbind(dailyReturn, dailyDifference)  
  
  x_lim <- range(datesProjection)
  y_lim <- range(matY)
  matplot( datesProjection ,  y = matY , type="l", lty = 1:2 , pch = 1:2 , col=c(1:2) , ylab="Return Values" , xlab = "Projection Dates" , ylim= y_lim, xlim= x_lim )
  legend("topleft", legend=graphNames, lty = 1:2 , pch = 1:2 , cex = 0.8 , col=c(1:2) )
  
  title(main=paste("Daily Risk Return vs Difference With replicating Portfolio"))
  dev.off()
}


##
## We would like here to use the results of our analysis in order to compare our derived portfolio with its synthetic or real one.
##
## isSPRealData = define if the base case is the real SP500 data or if it is a synthetic one
##
portfolio_projection <- function(coefResult, year_list_projection, isSPRealData = TRUE, basisType = "legendre", threshold = 0.2, freedom = 1){
  ## checking that the coefficient is not empty
  if(is.null(coefResult$lightCoef)){
    throw("[portfolio_projection]: the variable coefResult does not contain any data !")
  }
  
  ## extract the raw data from the files in our server
  constituentsFullDescription <- loadObject("sp500_data_list.RData")
  rownames(constituentsFullDescription) <- constituentsFullDescription[,1]
  constituentsDailyPrices <- read.csv("sp-500.csv")
  
  ## only keep the element that are within the coefficients results
  listOfComponents <- colnames(coefResult$lightCoef)
  listOfComponents_wo_intercept <- listOfComponents[2:length(listOfComponents)]
  constituentsDailyPrices_light <- constituentsDailyPrices[,c("X","Date", listOfComponents_wo_intercept)]
  
  ## extract the data corresponding the projection years
  mainDataXandY <- createTheMatrixXAndVectorY(constituentsDailyPrices_light, year_list_projection, freedom = freedom, basisType = basisType, threshold = threshold, useRealSPData = isSPRealData)
  
  if(isSPRealData){
    priceY <- mainDataXandY$SP500_prices
  }
  else{
    priceY <- mainDataXandY$Y$prices
  }
  
  pricesX <- mainDataXandY$keyData
  if(length(priceY) != nrow(pricesX)){
    throw("[portfolio_projection]: the size of both our instruments to be used for our projection is not valid!")
  }
  
  ## computing the weight corresponding to each instrument such that the sume is equal to 1
  ## normalizing the weight to 1
  coef <- coefResult$lightCoef
  percentage <- coef/sum(coef)
  
  
  ##
  ##  computing the return of both cases and then comparing the daily return
  ##
  returnY <- computePriceReturn(priceY)
  returnX <- computePriceReturn(pricesX)
  returnX <- cbind(matrix(1,nrow=nrow(returnX)), returnX)
  
  
  differenceReturn <- returnY - as.matrix(returnX) %*% t(percentage)
  result <- list()
  result$projectionDates <- mainDataXandY$keyDataFull[,"Date"]
  result$returnY <- returnY
  result$returnX <- returnX
  result$percentage <- percentage
  result$sd <- sd(differenceReturn)
  result$mean <- mean(differenceReturn)
  result$differenceReturn <- differenceReturn
  return(result)
}


## attempt to replicate the SP500 index return using a few amount of its constituents
run_sp500_index <- function(year_list, freedom = 1, basisType = "legendre", threshold = .2 , isContainError = FALSE, mu_error = 0.00, isCpp=FALSE){
  
  ## checking the validity of "years"
  if(!is.vector(year_list)){
    throw("[run_synthetic_example]: the list of years has to be a vector")
  }
  
  ## checking the validity of the elements in the list
  if(sum(is.na(year_list))!=0){
    throw("[run_synthetic_example]: the elements in the year list have to be double")
  }
  
  ## loading appropriate data - constituents prices and full name convention
  constituentsFullDescription <- loadObject("sp500_data_list.RData")
  rownames(constituentsFullDescription) <- constituentsFullDescription[,1]
  constituentsDailyPrices <- read.csv("sp-500.csv")
  
  ## creating the main matrix X to be use for our analysis
  mainDataXandY <- createTheMatrixXAndVectorY(constituentsDailyPrices, year_list, freedom = freedom, basisType = basisType, threshold = threshold, useRealSPData = TRUE, isContainError=isContainError, mu_error = mu_error)
  
  ## extracting the main matrix X and Y and computing the weight of our analysis
  returnY <- mainDataXandY$Y
  returnX <- mainDataXandY$X
  
  coefficientsLasso <- portfolioWeight(returnY, returnX, isContainError=isContainError, mu_error = mu_error, isCpp = isCpp)
  nonNullIndices <- which(coefficientsLasso != 0)
  coefficients_light <- matrix( coefficientsLasso[nonNullIndices] , nrow =1)
  
  if(isContainError){
    colnames(coefficients_light) <- names(coefficientsLasso)[nonNullIndices]
  }
  else{
    colnames(coefficients_light) <- dimnames(coefficientsLasso)[[1]][nonNullIndices]
  }
  coefLight_description <- constituentsFullDescription[colnames(coefficients_light)[2:length(nonNullIndices)],]
  
  
  ## storing all results
  result <-list()
  result$mainData <- mainDataXandY
  result$fullCoef <- coefficientsLasso
  result$lightCoef <- coefficients_light
  result$descriptionLightCoef <- coefLight_description
  result$yearList <- year_list
  return(result)
}

## creating  a synthetic portfolio from the available data and then running the procedure to see if we are able to reconstruct the initial portfolio
## year_list : is the years (a vector - may contain multiple data entry) on which we would want to run the procedure
run_synthetic_example <- function(year_list, freedom = 1, basisType = "legendre", threshold = .2, caseForY = "sector", nbStocksForY = 3){
  
  ## checking the validity of "years"
  if(!is.vector(year_list)){
    throw("[run_synthetic_example]: the list of years has to be a vector")
  }
  
  ## checking the validity of the elements in the list
  if(sum(is.na(year_list))!=0){
    throw("[run_synthetic_example]: the elements in the year list have to be double")
  }
  
  ## loading appropriate data - constituents prices and full name convention
  constituentsFullDescription <- loadObject("sp500_data_list.RData")
  rownames(constituentsFullDescription) <- constituentsFullDescription[,1]
  constituentsDailyPrices <- read.csv("sp-500.csv")
  
  ## creating the main matrix X to be use for our analysis
  mainDataXandY <- createTheMatrixXAndVectorY(constituentsDailyPrices, year_list, freedom = freedom, basisType = basisType, threshold = threshold, caseForY = caseForY, nbStocksForY = nbStocksForY)
  
  ## extracting the main matrix X and Y and computing the weight of our analysis
  returnY <- mainDataXandY$Y$return
  returnX <- mainDataXandY$X
  
  coefficientsLasso <- portfolioWeight(returnY, returnX)
  nonNullIndices <- which(coefficientsLasso != 0)
  coefficients_light <- matrix( coefficientsLasso[nonNullIndices] , nrow =1)
  
  colnames(coefficients_light) <- dimnames(coefficientsLasso)[[1]][nonNullIndices]
  coefLight_description <- constituentsFullDescription[colnames(coefficients_light)[2:length(nonNullIndices)],]
  
  ## extract the full description of the data used for this computation
  compositionY_description <- constituentsFullDescription[mainDataXandY$Y$composition,]
  
  ## storing all results
  result <-list()
  result$mainData <- mainDataXandY
  result$fullCoef <- coefficientsLasso
  result$lightCoef <- coefficients_light
  result$descriptionLightCoef <- coefLight_description
  result$compositionY <- compositionY_description
  result$yearList <- year_list
  return(result)
}


##
##
## Creating the vector Y (which is our synthetic portfolio) by either mixing up stock coming from different sectors
## or which have been choosen randomly
## Multiple cases are to be considered here :  
##
createTheVectorY <- function(keyData, case = "sector", nbStocks = 3){
  
  if(nbStocks > ncol(keyData) ){
    throw("[createTheVectorY]: the number of stocks to use is greater than the number of stocks available in the study")
  }
  
  ## initialization of final value
  result <- list()
  
  if( case == "random"){
    stockPosition <- as.integer(runif(nbStocks, min=1, max=ncol(keyData)))
  }
  else if(case =="sector"){  ## need to choose a sector randomly and then use it ot construct our synthetic portfolio
    constFullDescription <- loadObject("sp500_data_list.RData")
    
    ## only keep the ones which are within our own list of securities
    keyDataNames <- colnames(keyData)
    rownames(constFullDescription) <- constFullDescription[,1]
    constFullDescription_light <- constFullDescription[keyDataNames,]
  
    ## defining the list of sectors in our model
    sectorList <- unique(constFullDescription_light[,"GICS Sector"])
    ## Choosing the sector to pick
    sectorPick <- as.integer(runif(1, min=1, max=length(sectorList)))
    ## selecting the list of securities which are within the same sector
    sector_base_securities <- which(constFullDescription_light[,"GICS Sector"]==sectorList[sectorPick])
    
    ## choosing the final selection
    securityIndices <- as.integer(runif(nbStocks, min=1, max=length(sector_base_securities)))
    
    ## reusing definition used in the precedent case to create consistency and uniformity
    stockPosition <- sector_base_securities[securityIndices]
  }
  else{
    throw(paste("[createTheVectorY]:the sector specified is not currently available", case, sep=" "))
  }
  
  priceY = rowSums(keyData[,stockPosition])
  returnY = computePriceReturn(priceY)
  returnY_mat <- matrix(returnY, nrow=length(returnY), ncol=1)
  colnames(returnY_mat) <- paste( colnames(keyData)[stockPosition], collapse="+")
  
  ## storing the key values
  result$prices <- priceY
  result$return <- returnY_mat
  result$composition <-colnames(keyData)[stockPosition]
  return(result)
}

##
## creating the matrix X to be used for our LASSO computation
## we can use either the fourier or the legendre basis for our computation
## 1. threshold is needed in order to clean up the main data
##
createTheMatrixXAndVectorY <- function(mainData, yearList, freedom = 1, basisType = "legendre", threshold = .2, caseForY = "sector", nbStocksForY = 3, useRealSPData = FALSE, isContainError=FALSE, mu_error = 0.00){
  
  if(basisType != "legendre" && basisType != "fourier"){
    throw("[createTheMatrixX]: the basis type to used has to be either Fourier or Legendre")
  }
  
  ## filtering main data and only keeping the value of interest
  keyDataFull <- extract_data_givenListOfYears(mainData, yearList, threshold)
  nbCols <- ncol(keyDataFull)
  keyData <- keyDataFull[,3:nbCols]
  
  ## creating the vector of return linked to Y
  if(useRealSPData){
    realSP500Data <- read.csv("^GSPC.csv")
    rownames(keyDataFull) <- as.character(keyDataFull[,"Date"])
    rownames(realSP500Data) <- as.character(realSP500Data[,"Date"])
    realSP500Data_light_prices <- realSP500Data[as.character(keyDataFull[,"Date"]), "Close"]
    
    if(length(realSP500Data_light_prices) != length(keyDataFull[,"Date"])){
      throw("[createTheMatrixXAndVectorY]:they are some dates missing from the list of real SP500 data")
    }
    
    if(sum(is.na(realSP500Data_light_prices))!= 0){
      ## we need to adjust the keyDataFull and remove the dates which are not in both set
      indices_failures <- which(is.na(realSP500Data_light_prices))
      ## adjusting real SP500
      realSP500Data_light_prices <- realSP500Data_light_prices[-indices_failures]
      
      ## adjusting key data
      keyDataFull <- keyDataFull[-indices_failures,]
      keyData <- keyData[-indices_failures,]
    }
    
    ## computing the return given the list of prices
    resultY <- computePriceReturn(realSP500Data_light_prices)
  }
  else{
    resultY <- createTheVectorY(keyData, case=caseForY, nbStocks = nbStocksForY)
  }
  
  ## creating the matrix of return of X before constructing our basis function
  keyDataReturn <- computePriceReturn(keyData)
  
  ## adjusting the kayDataReturn 
  if(isContainError){
    nbRow <- nrow(keyDataReturn)
    error_vector <- rnorm(nbRow, mean=0, sd=mu_error)
    keyDataReturn <- keyDataReturn + error_vector
  }
  
  if(basisType=="legendre"){
    resultX <- legendre_extract(keyDataReturn, freedom)
  }
  else{
    resultX <- fourier_extract(keyDataReturn, freedom)
  }
  
  ## final return value
  result <- list()
  result$X <- resultX
  result$Y <- resultY
  result$keyData <- keyData
  result$keyDataFull <- keyDataFull
  if(useRealSPData){
    result$SP500_prices <- realSP500Data_light_prices
  }
  if(isContainError){
    result$error_vector <- error_vector
  }
  
  return(result)
}

## basic function in charge of applying either LASSO or Dantzip procedure in order to find the appropriate weight in a given portfolio
### Definition of varialbes
### Y == vector of daily return of either a synthetic portfolio or the SP500 itself
### X == the matrix containing the daily return per constituent in the portfolio - Each column corresponds to the return of a single constituent
### isLasso == defines whether we are using a LASSO based methdology or a DANTZIQ one.
### isLassoPenalize == defines whether we are penalizing the lasso coefficient or not
portfolioWeight <- function(Y, X, isLasso = TRUE, isLassoPenalize = FALSE , isContainError = FALSE, mu_error = 0.00, isCpp=FALSE){
  
  ## checking that there are no #NA values Y
  if(sum(is.na(Y)) != 0){
    throw("[portfolioWeight]:There is one(or multiple) values which are not double in the varialbe Y")
  }
  
  ## checking that htere are no #NA values in matrix X
  if(sum(is.na(X)) != 0){
    throw("[portfolioWeight]:There is one(or multiple) values which are not double in the varialbe X")
  }
  ## checking that X is a matrix
  if(!is.matrix(X)){
    throw("[portfolioWeight]:the variable X is not a matrix")
  }
  
  ## checking that X and Y have the same length
  if(nrow(X) != length(Y)){
    throw("[portfolioWeight]: the variable X and Y should have the same number of rows.")
  }
  
  #starting the calculation
  if(isLasso){
    if(!isLassoPenalize){
      if(isContainError){
        if(isCpp){
          lasso_coef <- cv.lasso_coordinate_descent_SP500_cpp(Y, X, mu_error)
        }else{
          lasso_coef <- cv.lasso_coordinate_descent_SP500(Y, X, mu_error)  
        }
        names(lasso_coef) <- colnames(X)
      }
      else{
        lasso.fit <- cv.glmnet(X ,Y, alpha = 1, intercept = TRUE, parallel = FALSE)
      }
    }
    else{
      penalty <- rep(1,ncol(X))
      ## only setting the function f variales to be set to 1
      penalty[1:iterW1] <- 0
      lasso.fit <- cv.glmnet(mat,Y, alpha = 1, intercept = TRUE, parallel = FALSE, penalty.factor = penalty)
    }
    
    if(!isContainError){
      lasso_coef <- coef(lasso.fit)
    }
    weights <- lasso_coef
    
    ## removing local variables
    if(!isContainError){
      rm(lasso.fit)
    }
  }
  else{
    dantzig.fit <- slim(X, Y, nlambda = 10,  method="dantzig", lambda.min.ratio=.01)
    dantzig_coeff <- rep(0,ncol(X))
    dantzigColLength <- ncol(dantzig.fit$beta)
    dantzig_coeff[0] <- dantzig.fit$intercept[,dantzigColLength]
    dantzig_coeff[1:ncol(X)] <- dantzig.fit$beta[,dantzigColLength]
    weights <- dantzig_coeff
  }

  ## final result
  return(weights)
}


## compute return given a matrix of prices
## returnType : two type will be available - discrete or log
## X : matrix/vector containing prices
computePriceReturn <- function( X , returnType = "discrete"){
  
  ## checking the matrix does not have any #NA values
  if(sum(is.na(X)) != 0){
    throw("[computePriceReturn]:The matrix contains values whihc are not double")
  }
  
  ## contain the value of the difference
  nbColumn <- 1
  nbRow <- 1
  isMatrix <- FALSE
  Xcolnames <- NULL
  
  if(is.matrix(X) || is.data.frame(X)){
    nbColumn <- ncol(X)
    nbRow <- nrow(X)
    isMatrix <- TRUE
    Xcolnames <- colnames(X)
  }
  else{
    nbRow <- length(X)
  }
  
  result <- matrix(0, ncol=nbColumn, nrow=nbRow-1)
  colnames(result) <- Xcolnames
  
  ### Starting the computation process...
  if(returnType == "discrete"){
    if(isMatrix){
      result <- ( X[2:nbRow,] - X[1:(nbRow-1),])/X[1:(nbRow-1),]
    }
    else{
      result <- ( X[2:nbRow] - X[1:(nbRow-1)])/X[1:(nbRow-1)]
    }
  }
  else if(returnType == "log"){
    if(isMatrix){
      result <- log( X[2:nbRow,]/X[1:(nbRow-1),])
    }
    else{
      result <- log( X[2:nbRow]/X[1:(nbRow-1)])
    }
  }
  else{
    throw(c("[computePriceReturn]: The return price provided is not available ", returnType))
  }
  
  return(result)
}

##
## extract data from main data set given a list of years
## and using a threshold to remove the column which have too many NA values
##
extract_data_givenListOfYears <- function(mainData, year_list, threshold = .5){
  #extracting date in the main data
  list_date_inData <- mainData[,"Date"]
  
  #extracting year in the main data
  list_date_year <- year(ymd(list_date_inData))
  
  ## indices corresponding to year in mainData
  indices <- is.element(list_date_year, year_list)
  
  ## removing all the rows not corresponding to the list of year provided
  intermData <- mainData[indices,]
  
  ## interpolating and removing all the columns which are not valid
  intermData <- interpolateMainDataAndRemoveNonValidColumns(intermData, threshold)
  
  finalData <- intermData[ rowSums(is.na(intermData)) == 0 , ]
  
  ## final values
  return(finalData)
}

##
## cleaning the main data frame by :
## 1. Removing the columns which have too many NA 
## 2. Replacing NA with the mean between two adjacents values.
##
interpolateMainDataAndRemoveNonValidColumns <- function(mainData, threshold = .5){
  ##checking the validity of mainData
  if(!is.matrix(matrix(mainData))){
    throw("[interpolateMainData]:the variable mainData should be a matrix !");
  }
  
  ## ensuring that the threshold is a proba value
  if(threshold >1 || threshold < 0){
    throw("[interpolateMainData]:the variable threshold should be between 0 and 1")
  }
  
  nbCols <- ncol(mainData)
  nbRows <- nrow(mainData)
  startOfMatrix <- 3
  #result <- matrix(mainData, ncol = nbCols, nrow = nbRows)
  result <- mainData
  isThereElementToRemove <- FALSE
  increment <- 0
  
  ## removing the columns which have too many na values
  for(i in startOfMatrix:nbCols){
    percentOfNA <- sum(is.na(result[,i]))/nbRows
    if(percentOfNA > threshold){
      if(increment == 0){
        indicesColToRemove <- -1 * i
        isThereElementToRemove <- TRUE
      }
      else{
        indicesColToRemove <- c( -1*i , indicesColToRemove)
      }
      increment <- increment + 1
    }
  }
  
  if(isThereElementToRemove){
    result <- result[,indicesColToRemove]
  }
  
  nbCols <- ncol(result)
  ## applying interpolation
  for(i in startOfMatrix:nbCols){
    indices <- is.na(result[,i])
    for(j in 1:nbRows){
      if(indices[j]){
        if(j!=1 && j!=nbRows){
          result[j,i] = ( result[j-1,i] + result[j+1,i])/2
        }
      }
    }
  }
  
  ## return final values
  return(result)
}


##
## given a dataset, an a sequence for lambda, i would like to determine the lambda which will provide me the best fit possible for my optimization problem
## the case develop here is specific to the SP500 case study here... Furthermore, it is extremely simplify as we assume that the legendre decomposition has a degree of freedom equal to 1
## as such, we can safely add the error term to the main matrix of return and then used the function below to infer the beta value...
##
cv.lasso_coordinate_descent_SP500 <- function(Y, returnMatrix, mu, lambda_max = .03, lambda_min = 0, lambda_seq = .005, maxIter = 100, tolerance = 1e-3, nbFolds = 10, mcCores = 4){
  #
  # checking the validity of few variables
  #
  if(lambda_max <0){
    throw("[cv.lasso_coordinate_descent_SP500]:the variable lambda_max must be positive")
  }
  
  if(lambda_seq < 0){
    throw("[cv.lasso_coordinate_descent_SP500]:the variable lambda_seq must be positive")
  }
  
  if( length(Y) != nrow(returnMatrix)){
    throw("[cv.lasso_coordinate_descent_SP500]:the variables Y must have the same length as the number of row of the matrix returnMatrix")
  }
  
  ## computing key relevant variables
  Z <- returnMatrix
  H <- matrix(1, ncol=ncol(Z), nrow = nrow(Z))
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
      
      beta <- lasso_coordinate_descent_givenLambda(Y_training, Z_training, mu, H_training, lambda, beta0 = NULL, maxIter = maxIter, tolerance = tolerance)
      
      MSE <- (sum((Y[testIndices] - Z[testIndices,] %*%beta)^2)  + mu^2 * sum( (H[testIndices,]%*%beta)^2))/length(Y[testIndices])
      return(MSE)
    }
    
    print(paste("I have done ", i, " out of ", length(list_lambda), sep=""))
    
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
  
  print(biaisVariance)
  
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


cv.lasso_coordinate_descent_SP500_cpp <- function(Y, returnMatrix, mu, lambda_max = .03, lambda_min = 0, lambda_seq = .005, maxIter = 100, tolerance = 1e-3, nbFolds = 10, mcCores = 1){
  #
  # checking the validity of few variables
  #
  if(lambda_max <0){
    throw("[cv.lasso_coordinate_descent_SP500_cpp]:the variable lambda_max must be positive")
  }
  
  if(lambda_seq < 0){
    throw("[cv.lasso_coordinate_descent_SP500_cpp]:the variable lambda_seq must be positive")
  }
  
  if( length(Y) != nrow(returnMatrix)){
    throw("[cv.lasso_coordinate_descent_SP500_cpp]:the variables Y must have the same length as the number of row of the matrix returnMatrix")
  }
  
  ## computing key relevant variables
  Z <- returnMatrix
  H <- matrix(1, ncol=ncol(Z), nrow = nrow(Z))
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
      
      beta <- rep(0, ncol(Z))
      lasso_coordinate_descent_givenLambda_cpp(Y_training, Z_training, mu, H_training, lambda, beta, maxIter = maxIter, tolerance = tolerance)
      
      MSE <- (sum((Y[testIndices] - Z[testIndices,] %*%beta)^2)  + mu^2 * sum( (H[testIndices,]%*%beta)^2))/length(Y[testIndices])
      return(MSE)
    }
    
    print(paste("CPP: I have done ", i, " out of ", length(list_lambda), sep=""))
    
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
  
  print(biaisVariance)
  
  ## computing beta using the best lambda we found
  bestLambda <- list_lambda[lambdaSolutionIndice]
  result <- rep(0, ncol(Z))
  lasso_coordinate_descent_givenLambda_cpp(Y, Z, mu, H, bestLambda, result, maxIter = maxIter, tolerance = tolerance)
  
  ##
  ## removing local variables
  ##
  rm(MSE_output, list_lambda,lambda_based_MSE, biaisVariance)
  
  ## returning final values
  #browser()
  return(result)
}


##
## plotting the daily return corresponding to the SP500, as well as the difference with the replicating portfolio from four different cases
##

plotDifferenceWithReplicatingPortfolioFromDifferentMu <- function(coefResult_list, yearProjection){
  
  nbOfCoef <- length(coefResult_list)
  nbItemsInPic <- nbOfCoef + 1
  result_projection <- list()
  dailyDifference <- list()
  
  for(i in 1:nbOfCoef){
    result_projection[[i]] <- portfolio_projection(coefResult_list[[i]], yearProjection)
    dailyDifference[[i]] <- result_projection[[i]]$differenceReturn
    print(paste("i am at stage ", i))
    
    
    ## need to merge dates and create the matrix Y
    if(i==1){
      dailyReturn <- result_projection[[i]]$returnY
      datesProjection <- as.Date(result_projection[[1]]$projectionDates[1:length(dailyReturn)])
      matY <- matrix(result_projection[[i]]$differenceReturn , ncol=1)
    }
    else{
      newDailyReturn <- result_projection[[i]]$returnY
      testIndices <- (datesProjection %in% as.Date(result_projection[[i]]$projectionDates[1:length(newDailyReturn)]))
      testIndicesOut <- ( as.Date(result_projection[[i]]$projectionDates[1:length(newDailyReturn)]) %in%  datesProjection )
      
      datesProjection <- datesProjection[testIndices]
      dailyReturn <- dailyReturn[testIndices]
      matY <- cbind(result_projection[[i]]$differenceReturn[testIndicesOut] , matY[testIndices,])
    }
  }
  
  pathFileName <- paste("Comparison_Evolution_DailyReturn_per_value_of_mu_" , as.character(yearProjection), "_file.pdf", sep="")
  pdf(file=pathFileName)
  
  ##  creating the graph
  #graphNames <- c("SP500 Daily return","Diff. of daily return betw. SP500 and replica. portf. mu = 0", "Diff. of daily return betw. SP500 and replica. portf. mu = 5e-3", "Diff. of daily return betw. SP500 and replica. portf. mu = 1e-2", "Diff. of daily return betw. SP500 and replica. portf. mu = 1.5e-2")
  graphNames <- c("SP500 Daily return","Diff. of daily return betw. SP500 and replica. portf. mu = 0", "Diff. of daily return betw. SP500 and replica. portf. mu = 5e-3")
  
  #matY <- matrix( unlist(dailyDifference) , ncol = nbOfCoef)
  matY <- cbind(dailyReturn, matY)  
  
  x_lim <- range(datesProjection)
  y_lim <- range(matY)
  
  matplot( datesProjection ,  y = matY , type="l", lty = 1:nbItemsInPic , pch = 1:nbItemsInPic , col=c(1:nbItemsInPic) , ylab="Return Values" , xlab = "Projection Dates" , ylim= y_lim, xlim= x_lim )
  legend("topleft", legend=graphNames, lty = 1:nbItemsInPic , pch = 1:nbItemsInPic , cex = 0.8 , col=c(1:nbItemsInPic) )
  
  title(main=paste("Daily Risk Return vs Difference With replicating Portfolio"))
  dev.off()
}
