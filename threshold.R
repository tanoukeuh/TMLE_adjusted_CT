threshold <- function(xx, min=0, max=1-min) {
  
  #<Cabral>
  min <- Arguments$getNumerics(min);
  max <- Arguments$getNumerics(max);
  
  xx <- Arguments$getNumerics(xx);
  if(is.matrix(xx)){
    out <- sapply(1:nrow(xx), function(ii){pmin(max, pmax(min, xx[ii,]))})
    out <- t(out)
  }
  else{
    out <- pmin(max, pmax(min, xx))
  }  
  #<\Cabral>
}

############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

