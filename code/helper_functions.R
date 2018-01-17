library(quantreg)
library(dplyr)

## BLOM ESTIMATOR
## Equivalent to quantile(x, probs, type = 6)

#Calculate the Blom estimator for Q3 values for all length categories
BlomEstimator <- function(vec, percentile){
    #The Blom method here is as described from Gerow 2009
    vec <- sort(vec)
    #Identify which value in vec is at the percentile requested
    nDf <- percentile*(length(vec)+1)
    #round the number down to the nearest integer
    L <- floor(nDf)
    #Calculate difference
    s <- (nDf-L)
    #Calculate the percentile value
    return(vec[L]+(s*(vec[L+1]-vec[L])))

    }


## BOOTSTRAP ESTIMATES OF SLOPE AND INTERCEPT

EstLinQuantRegCoefDist <- function(data, n, numSamples, tau = 0.75){
  
    regCoefs <- 
      as.data.frame(matrix(NA, nrow = n, ncol=2)) %>%
      rename(intercept = V1, 
             slope = V2)      
  
    for(i in 1:n){
      rowVec <- sample(1:nrow(data), numSamples, replace=T)
      tmp <- data[rowVec, ]
      mod <- rq(log10(weight)~log10(length), data=tmp, tau = tau)
    
      regCoefs[i,1] <- mod$coefficients[1]
      regCoefs[i,2] <- mod$coefficients[2]
    }
    
    names(regCoefs) <- c("intercept", "slope")
    return(regCoefs)
  }

# Calculate standard error of a vector

std_error <- function(x){
  sd(x)/sqrt(length(x))
}

# Helper for assinging length classes

round_down <- function(x,to=10)
{
  to*(x %/% to)
}