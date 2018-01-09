  setwd('C:/Users/Steven Harris Ranney/Desktop/Manuscripts/Ws - Quantiles/code')

  wae <- read.table('WAE Clean.txt', header=T)

  dF <- wae
  dfNames <- c("lake", "length", "weight")
  n = 10000
  numSamples = 150
  tau = 0.75

  require(quantreg)
  
  EstLinQuantRegCoefDist <- function(dF, n, numSamples, tau = 0.75){
    
    regCoefs <- as.data.frame(matrix(NA, nrow = n, ncol=2))
    
    for(i in 1:n){
      rowVec <- sample(1:nrow(dF), numSamples, replace=T)
      tmp <- dF[rowVec, ]
      mod <- rq(log10(weight)~log10(length), data=tmp, tau = tau)
      
      regCoefs[i,1] <- mod$coefficients[1]
      regCoefs[i,2] <- mod$coefficients[2]
    }
    names(regCoefs) <- c("intercept", "slope")
    return(regCoefs)
  }

  refDist75 <- EstLinQuantRegCoefDist(wae, 20000, 100, 0.75)
  normDist <- rnorm(20000, mean(refDist75[,2]), sd(refDist75[,2]))  
  ks.test(refDist75[,2], normDist)
  

  hist(refDist75[,2], breaks=100, col="gray")
  MakeColorsTransparent("red")
  hist(normDist, breaks=100, add=T, col=MakeColorsTransparent("red"))
  
  mean(refDist75[,2])
  mean(refDist50[,2])
  
  
  
  
  #indData
  waeInd <- read.table('WAE_independent.txt', header=T)

  #Define GA and SD populations
  
  waeGA <- waeInd[waeInd$State == "GA", ]
 
  ga2Dist <- EstLinQuantRegCoefDist(waeGA[waeGA$lake == "2",], 10000, 100, 0.75)
  ga3Dist <- EstLinQuantRegCoefDist(waeGA[waeGA$lake == "3",], 10000, 100, 0.75)
  ga4Dist <- EstLinQuantRegCoefDist(waeGA[waeGA$lake == "4",], 10000, 100, 0.75)
  
  
  