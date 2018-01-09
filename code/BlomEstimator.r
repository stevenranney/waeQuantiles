
  tmp <- rnorm(50, 2855, 20.34)

  quantile(tmp, .75, type=4)
  quantile(tmp, .75, type=5)
  quantile(tmp, .75, type=6)
  quantile(tmp, .75, type=7)
  quantile(tmp, .75, type=8)
  quantile(tmp, .75, type=9)


  #Blom estimator
  tmp <- rnorm(33597, 1035.698, 940.428)
  
  tmp <- sort(tmp)
  
  P <- c(.25)#, .5, .75)
  
  nDf <- P*(length(tmp)+1)
  
  L <- round(nDf, 0)
  
  s <- abs(nDf-L)
  
  tmp[L]+(s*(tmp[L+1]-tmp[L]))
  
  #Blom estimator in R = quantile(x, .75, type=6)
  quantile(tmp, c(.25), type=6)

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
 
