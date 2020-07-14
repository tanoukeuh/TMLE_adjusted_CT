#############################################################################################################
##
##  Cabral CHANANG - December 2017
## 
#############################################################################################################
##  Description of the file :
##
##  will be used as a triggering mechanism for our simulation - updated for Antoine
##
#############################################################################################################


##
## running all the simulation
##
run_all_procedure_ForAntoine <- function(){
  nbCores <- 4
  listMu <- seq(0,.4,.05)
  
  ## with error simulation
  for(i in 1:length(listMu)){
    mu_i <- listMu[i]
    run_parallel_lassoWithError_lowDimension_SingleIteration(nbCores, mu_i)
    run_parallel_lassoWithError_highDimension_SingleIteration(nbCores, mu_i)
  }
  
  ## TMLE comparison
  run_parallel_tmle_lowDimension(nbCores)    ## only run a single simulation
  run_parallel_tmle_highDimension(nbCores)   ## only run a single simulation
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
run_parallel_lassoWithError_lowDimension_SingleIteration <- function(nbCores, mu_i){
  run_tmle_algo <- FALSE
  runLassoWithError <- TRUE
  ff <- "sine"
  gg <- "identity"
  lengthBaseW <- 20
  lengthBaseX <- 20
  nbCores <- nbCores
  mu <- mu_i
  
  run_parallel(ff,gg,runLassoWithError=runLassoWithError, mu=mu_i, nbCores = nbCores, lengthBaseW = lengthBaseW, lengthBaseX = lengthBaseX)
}

run_parallel_lassoWithError_highDimension_SingleIteration <- function(nbCores, mu_i){
  run_tmle_algo <- FALSE
  runLassoWithError <- TRUE
  ff <- "sine"
  gg <- "identity"
  lengthBaseW <- 200
  lengthBaseX <- 50
  nbCores <- nbCores
  
  mu <- mu_i
  run_parallel(ff,gg,runLassoWithError=runLassoWithError, mu=mu, nbCores = nbCores, lengthBaseW = lengthBaseW, lengthBaseX = lengthBaseX)
}