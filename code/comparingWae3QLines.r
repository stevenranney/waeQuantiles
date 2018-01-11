################################################################################
# This code is annotated below.  It 
################################################################################ 
  
  require(quantreg)
  require(ggplot2)
  library(purrr)
  library(dplyr)
  library(gplots)

  set.seed(256)

  setwd('C:/Users/Steven Harris Ranney/Desktop/Manuscripts/Ws - Quantiles/code')

  #read in and manipulate the WAE dataset from Ranney et al. 2010 and 2011; data
  # filtering protocols have already been applied
  wae <- read.table('WAE Clean.txt', header=T)
  attach(wae)
  names(wae)


  #Assigns fish to a length category (Gabelhouse 1984)
  psd <- with(wae,
              ifelse((length>=250)&(length<380), "S-Q",
              ifelse((length>=380)&(length<510), "Q-P",
              ifelse((length>=510)&(length<630), "P-M",
              ifelse((length>=630)&(length<760), "M-T",
              ifelse(length>=760, ">T",
              "SS"))))))

  #Assigns fish to 10-mm length class
  l.c <- with(wae,
              ifelse((length>=160)&(length<170),"165",
              ifelse((length>=170)&(length<180),"175",
              ifelse((length>=180)&(length<190),"185",
              ifelse((length>=190)&(length<200),"195",
              ifelse((length>=200)&(length<210),"205",
              ifelse((length>=210)&(length<220),"215",
              ifelse((length>=220)&(length<230),"225",
              ifelse((length>=230)&(length<240),"235",
              ifelse((length>=240)&(length<250),"245",
              ifelse((length>=250)&(length<260),"255",
              ifelse((length>=260)&(length<270),"265",
              ifelse((length>=270)&(length<280),"275",
              ifelse((length>=280)&(length<290),"285",
              ifelse((length>=290)&(length<300),"295",
              ifelse((length>=300)&(length<310),"305",
              ifelse((length>=310)&(length<320),"315",
              ifelse((length>=320)&(length<330),"325",
              ifelse((length>=330)&(length<340),"335",
              ifelse((length>=340)&(length<350),"345",
              ifelse((length>=350)&(length<360),"355",
              ifelse((length>=360)&(length<370),"365",
              ifelse((length>=370)&(length<380),"375",
              ifelse((length>=380)&(length<390),"385",
              ifelse((length>=390)&(length<400),"395",
              ifelse((length>=400)&(length<410),"405",
              ifelse((length>=410)&(length<420),"415",
              ifelse((length>=420)&(length<430),"425",
              ifelse((length>=430)&(length<440),"435",
              ifelse((length>=440)&(length<450),"445",
              ifelse((length>=450)&(length<460),"455",
              ifelse((length>=460)&(length<470),"465",
              ifelse((length>=470)&(length<480),"475",
              ifelse((length>=480)&(length<490),"485",
              ifelse((length>=490)&(length<500),"495",
              ifelse((length>=500)&(length<510),"505",
              ifelse((length>=510)&(length<520),"515",
              ifelse((length>=520)&(length<530),"525",
              ifelse((length>=530)&(length<540),"535",
              ifelse((length>=540)&(length<550),"545",
              ifelse((length>=550)&(length<560),"555",
              ifelse((length>=560)&(length<570),"565",
              ifelse((length>=570)&(length<580),"575",
              ifelse((length>=580)&(length<590),"585",
              ifelse((length>=590)&(length<600),"595",
              ifelse((length>=600)&(length<610),"605",
              ifelse((length>=610)&(length<620),"615",
              ifelse((length>=620)&(length<630),"625",
              ifelse((length>=630)&(length<640),"635",
              NA)))))))))))))))))))))))))))))))))))))))))))))))))
  
  l.c.1 <- with(wae,
                ifelse((length>=640)&(length<650),"645",
                ifelse((length>=650)&(length<660),"655",
                ifelse((length>=660)&(length<670),"665",
                ifelse((length>=670)&(length<680),"675",
                ifelse((length>=680)&(length<690),"685",
                ifelse((length>=690)&(length<700),"695",
                ifelse((length>=700)&(length<710),"705",
                ifelse((length>=710)&(length<720),"715",
                ifelse((length>=720)&(length<730),"725",
                ifelse((length>=730)&(length<740),"735",
                l.c)))))))))))

  #combines the PSD values and lengthClass categories to the original dataset
  wae <- cbind(wae, psd, l.c.1)

  #Build nonlinear quantile regression models;
  #USED ONLY FOR CALCULATING A TABLE OF PREDICTED QUANTILE VALUES
  wae.mod.75 <- nlrq(weight~alpha*length^beta, data=wae, tau=0.75,                         #Requires different parm. ests. for alpha and beta
                     start=list(alpha=0.00001, beta=3))
  wae.mod.5 <- nlrq(weight~alpha*length^beta, data=wae, tau=0.5,                           #Requires different parm. ests. for alpha and beta
                    start=list(alpha=0.00001, beta=3))
  wae.mod.25 <- nlrq(weight~alpha*length^beta, data=wae, tau=0.25,                         #Requires different parm. ests. for alpha and beta
                     start=list(alpha=0.00001, beta=3))

  
  #Build linear quantile regression models
  wae.linMod.75 <- rq(log10(weight)~log10(length), data=wae, tau=0.75)
  wae.linMod.5 <- rq(log10(weight)~log10(length), data=wae, tau=0.5)
  wae.linMod.25 <- rq(log10(weight)~log10(length), data=wae, tau=0.25)

 
  #Create a vector of lengths to include in plots
  x <- seq(165,735, by=10)#; x

  #plot data and all lines
#  plot(weight~length, data=wae, xlab="Length (mm)", xlim=c(100, 800), ylim=c(0,6000), 
#       ylab="Weight (g)", pch="", xaxs="i", yaxs="i", bty="n")

  #Create a vector to use in plotting the quantile models from above
#  wae.moda <- seq(150,740)
#    lines(wae.moda, predict(wae.mod.75, list(length=wae.moda)), lwd=2, lty=2)
#    lines(wae.moda, predict(wae.mod.5, list(length=wae.moda)), lwd=2)
#    lines(wae.moda, predict(wae.mod.25, list(length=wae.moda)), lwd=2, lty=2)
  #add midpoints of length categories to plots
#  points(tapply(weight, l.c.1, mean)~x, pch=1)
  #add error bars to points, 95%CI
#  errbar(x, 
#         tapply(weight, l.c.1, mean), 
#         tapply(weight, l.c.1, mean)+(1.96*tapply(weight, l.c.1, sd)), 
#         tapply(weight, l.c.1, mean)-(1.96*tapply(weight, l.c.1, sd)), 
#         pch=" ", add=T)

  #Add text to plot 
#  text(250, 4000, labels=expression(bold(Walleye)), cex=1.5)

  #Create vector of predicted values from each model
  wae.pred.75 <- as.vector(predict(wae.mod.75, list(length=seq(155,745, by=10)), interval="confidence"))
  wae.pred.5 <- as.vector(predict(wae.mod.5, list(length=seq(155,745, by=10)), interval="confidence"))
  wae.pred.25 <- as.vector(predict(wae.mod.25, list(length=seq(155,745, by=10)), interval="confidence"))

  #Create a table of lengths and predicted values for 25th, 50th, and 75th percentiles

  #creates a matrix 60 rows long by 3 columsn wide
  wae.pred.values <- matrix(NA,60,3, dimnames=list(seq(155,745,by=10), c("Q1", "Q2", "Q3")))
  #program to put all predicted values from each model into the matrix "pred.values"
  for (i in 1:60){
    wae.pred.values[i,1]=as.vector(wae.pred.25[i])
    wae.pred.values[i,2]=as.vector(wae.pred.5[i])
    wae.pred.values[i,3]=as.vector(wae.pred.75[i])
  }
  
  wae.pred.values <- as.data.frame(wae.pred.values, header=T)
  names(wae.pred.values) <- c("25th", "50th", "75th")
#  wae.pred.values 


  #create a table of intercept and slope values from each nlrq 
  #THESE ARE NONLINEAR SLOPES FROM THE NLRQ MODELS ABOVE
#  wae.int.slope <- matrix(NA,3,2, dimnames=list(c("25th", "50th", "75th"), c("Intercept", "Slope")))
#    wae.int.slope[1,1]=summary(wae.mod.25)$coef[1,1]
#    wae.int.slope[1,2]=summary(wae.mod.25)$coef[2,1]
#    wae.int.slope[2,1]=summary(wae.mod.5)$coef[1,1]
#    wae.int.slope[2,2]=summary(wae.mod.5)$coef[2,1]
#    wae.int.slope[3,1]=summary(wae.mod.75)$coef[1,1]
#    wae.int.slope[3,2]=summary(wae.mod.75)$coef[2,1]

#  wae.int.slope


#-------------------------------------------------------------------------------
  #Read in independent data
  
  waeInd <- read.table('WAE_independent.txt', header=T)

  #Assigns fish to a length category (Gabelhouse 1984)
  psd <- with(waeInd,
              ifelse((length>=250)&(length<380), "S-Q",
              ifelse((length>=380)&(length<510), "Q-P",
              ifelse((length>=510)&(length<630), "P-M",
              ifelse((length>=630)&(length<760), "M-T",
              ifelse(length>=760, ">T",
              "SS"))))))
  
  waeInd <- cbind(waeInd, psd)
              
  #Define GA and SD populations
  
  waeGA <- waeInd[waeInd$State == "GA", ]
  waeSD <- waeInd[waeInd$State == "SD", ]
  
  #plot 3Q line and GA and SD data
#  for(i in unique(waeGA$lake)){
#    jpeg(paste0("gaPlot", i, ".jpg"), quality=100, res=100, units="in", 
#         width=6.5, height=6.5)
#    plot(weight~length, data=wae, xlab="Length (mm)", xlim=c(100, 800), ylim=c(0,6000), 
#        ylab="Weight (g)", pch="", xaxs="i", yaxs="i", bty="n")
#    wae.moda <- seq(150,740)
#      lines(wae.moda, predict(wae.mod.75, list(length=wae.moda)), lwd=2, lty=2, col="red")
#      lines(wae.moda, predict(wae.mod.5, list(length=wae.moda)), lwd=2, col="red")
#      lines(wae.moda, predict(wae.mod.25, list(length=wae.moda)), lwd=2, lty=2, col="red")
  
#    points(weight~length, data=waeGA[waeGA$lake == i,])
#    dev.off()
#  }
  
#  for(i in unique(waeSD$lake)){
#    jpeg(paste0("sdPlot", i, ".jpg"), quality=100, res=100, units="in", 
#         width=6.5, height=6.5)
#    plot(weight~length, data=wae, xlab="Length (mm)", xlim=c(100, 800), ylim=c(0,6000), 
#        ylab="Weight (g)", pch="", xaxs="i", yaxs="i", bty="n")
#    wae.moda <- seq(150,740)
#      lines(wae.moda, predict(wae.mod.75, list(length=wae.moda)), lwd=2, lty=2, col="red")
#      lines(wae.moda, predict(wae.mod.5, list(length=wae.moda)), lwd=2, col="red")
#      lines(wae.moda, predict(wae.mod.25, list(length=wae.moda)), lwd=2, lty=2, col="red")
  
#    points(weight~length, data=waeSD[waeSD$lake == i,])
#    dev.off()
#  }

  #From the plots above, GA lakes 2, 3, and 4 have pops that may have Q3 higher 
  #than or equal to the ref 
  ga2 <- rq(log10(weight)~log10(length), data=waeGA[waeGA$lake == "2",], tau=0.75)
  ga3 <- rq(log10(weight)~log10(length), data=waeGA[waeGA$lake == "3",], tau=0.75)
  ga4 <- rq(log10(weight)~log10(length), data=waeGA[waeGA$lake == "4",], tau=0.75)

#  plot(log10(weight)~log10(length), data=wae, xlab="log10Length (mm)", ylab="log10Weight (g)", 
#      pch="", xaxs="i", yaxs="i", bty="n", xlim = c(2.5,2.6), ylim=c(2.4,2.8))
  
#  abline(wae.linMod.75)
#  abline(ga2, lty=2)
#  abline(ga3, lty=3)
#  abline(ga4, lty=4)
#  legend("bottomright", legend=c("Reference", "Population 2", "Population 3", 
#         "Population 4"), lty=c(1, 2, 3, 4), bty="n")
  
  
  par(mfrow=c(3,1))
    plot(log10(weight)~log10(length), data=wae, xlab="log10Length (mm)", ylab="log10Weight (g)", 
      pch="", xaxs="i", yaxs="i", bty="n", xlim = c(2.5,2.6), ylim=c(2.0,3.0))
    abline(wae.linMod.75)
    abline(ga2, col="red")
    plot(log10(weight)~log10(length), data=wae, xlab="log10Length (mm)", ylab="log10Weight (g)", 
      pch="", xaxs="i", yaxs="i", bty="n", xlim = c(2.5,2.6), ylim=c(2.0,3.0))
    abline(wae.linMod.75)
    abline(ga3, lty=2, col="red")
    plot(log10(weight)~log10(length), data=wae, xlab="log10Length (mm)", ylab="log10Weight (g)", 
      pch="", xaxs="i", yaxs="i", bty="n", xlim = c(2.5,2.6), ylim=c(2.0,3.0))
    abline(wae.linMod.75)
    abline(ga4, col="red")
  par(mfrow=c(1,1))

  #Plot all lines on the same figure
  width = 1.5
  plot(log10(weight)~log10(length), data=wae, xlab="log10Length (mm)", ylab="log10Weight (g)", 
       pch="", xaxs="i", yaxs="i", bty="n", xlim = c(2.5,2.6), ylim=c(2.4,2.8))
  abline(wae.linMod.75, lty=1, lwd=width)
  abline(ga2, lty=2, lwd=width)
  abline(ga3, lty=4, lwd=width)
  abline(ga4, lty=5, lwd=width)
  legend("bottomright", bty="n", lty=c(1, 2, 4, 5), lwd = width, 
         legend=c("Reference", "GA 1", "GA 2", "GA 3"))
  

  #From the plots above, SD lakes 4, 13, and 25 have pops that may have Q3 higher 
  #than or equal to the ref 
  sd4 <- rq(log10(weight)~log10(length), data=waeSD[waeSD$lake == "4",], tau=0.75)
  sd13 <- rq(log10(weight)~log10(length), data=waeSD[waeSD$lake == "13",], tau=0.75)
  sd25 <- rq(log10(weight)~log10(length), data=waeSD[waeSD$lake == "25",], tau=0.75)

#  par(mfrow=c(3,1))
#  plot(log10(weight)~log10(length), data=wae, xlab="log10Length (mm)", ylab="log10Weight (g)", 
#      pch="", xaxs="i", yaxs="i", bty="n", xlim = c(2.5,2.6), ylim=c(2.0,3.0))
#  abline(wae.linMod.75)
#  abline(sd4, col="red")
#    plot(log10(weight)~log10(length), data=wae, xlab="log10Length (mm)", ylab="log10Weight (g)", 
#      pch="", xaxs="i", yaxs="i", bty="n", xlim = c(2.5,2.6), ylim=c(2.0,3.0))
#  abline(wae.linMod.75)
#  abline(sd13, lty=2, col="red")
#    plot(log10(weight)~log10(length), data=wae, xlab="log10Length (mm)", ylab="log10Weight (g)", 
#      pch="", xaxs="i", yaxs="i", bty="n", xlim = c(2.5,2.6), ylim=c(2.0,3.0))
#  abline(wae.linMod.75)
#  abline(sd25, col="red")
#  par(mfrow=c(1,1))

  #Plot all lines on the same figure
  plot(log10(weight)~log10(length), data=wae, xlab="log10Length (mm)", ylab="log10Weight (g)", 
       pch="", xaxs="i", yaxs="i", bty="n", xlim = c(2.5,2.6), ylim=c(2.4,2.8))
  abline(wae.linMod.75, lty=1, lwd=width)
  abline(sd4, lty=2, lwd=width)
  abline(sd13, lty=4, lwd=width)
  abline(sd25, lty=5, lwd=width)
  legend("bottomright", bty="n", lty=c(1, 2, 4, 5), legend=c("Reference", "SD 1", 
                                                             "SD 2", "SD 3"))
  
  #-----------------------------------------------------------------------------

  #Calculate the Blom estimator for Q3 values for all length categories
    BlomEstimator <- function(vec, percentile){
      #The Blom method here is as described from Gerow 2009
      #Sort the vector from smallest to largest
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
    
  sSRef <- BlomEstimator(wae$weight[wae$psd=="SS"], c(.25, .5,.75))
  sQRef <- BlomEstimator(wae$weight[wae$psd=="S-Q"], c(.25, .5,.75))
  qPRef <- BlomEstimator(wae$weight[wae$psd=="Q-P"], c(.25, .5,.75))
  pMRef <- BlomEstimator(wae$weight[wae$psd=="P-M"], c(.25, .5,.75))
  mTRef <- BlomEstimator(wae$weight[wae$psd=="M-T"], c(.25, .5,.75))
  gTRef <- c(NA, NA, NA)
 
  #GA lake 2
  sSGa2 <- BlomEstimator(waeGA$weight[waeGA$psd=="SS" & waeGA$lake == "2"], c(.25, .5, .75))
  sQGa2 <- BlomEstimator(waeGA$weight[waeGA$psd=="S-Q" & waeGA$lake == "2"], c(.25, .5, .75))
  qPGa2 <- BlomEstimator(waeGA$weight[waeGA$psd=="Q-P" & waeGA$lake == "2"], c(.25, .5, .75))
  pMGa2 <- BlomEstimator(waeGA$weight[waeGA$psd=="P-M" & waeGA$lake == "2"], c(.25, .5, .75))
  mTGa2 <- BlomEstimator(waeGA$weight[waeGA$psd=="M-T" & waeGA$lake == "2"], c(.25, .5, .75))
  gTGa2 <- c(NA, NA, NA)
  
  #GA lake 3
  sSGa3 <- BlomEstimator(waeGA$weight[waeGA$psd=="SS" & waeGA$lake == "3"], c(.25, .5, .75))
  sQGa3 <- BlomEstimator(waeGA$weight[waeGA$psd=="S-Q" & waeGA$lake == "3"], c(.25, .5, .75))
  qPGa3 <- BlomEstimator(waeGA$weight[waeGA$psd=="Q-P" & waeGA$lake == "3"], c(.25, .5, .75))
  pMGa3 <- BlomEstimator(waeGA$weight[waeGA$psd=="P-M" & waeGA$lake == "3"], c(.25, .5, .75))
  mTGa3 <- BlomEstimator(waeGA$weight[waeGA$psd=="M-T" & waeGA$lake == "3"], c(.25, .5, .75))
  gTGa3 <- c(NA, NA, NA)

  #GA lake 4
  sSGa4 <- BlomEstimator(waeGA$weight[waeGA$psd=="SS" & waeGA$lake == "4"], c(.25, .5, .75))
  sQGa4 <- BlomEstimator(waeGA$weight[waeGA$psd=="S-Q" & waeGA$lake == "4"], c(.25, .5, .75))
  qPGa4 <- BlomEstimator(waeGA$weight[waeGA$psd=="Q-P" & waeGA$lake == "4"], c(.25, .5, .75))
  pMGa4 <- BlomEstimator(waeGA$weight[waeGA$psd=="P-M" & waeGA$lake == "4"], c(.25, .5, .75))
  mTGa4 <- BlomEstimator(waeGA$weight[waeGA$psd=="M-T" & waeGA$lake == "4"], c(.25, .5, .75))
  gTGa4 <- c(NA, NA, NA)

  #SD lake 4
  sSSd4 <- BlomEstimator(waeSD$weight[waeSD$psd=="SS" & waeSD$lake == "4"], c(.25, .5, .75))
  sQSd4 <- BlomEstimator(waeSD$weight[waeSD$psd=="S-Q" & waeSD$lake == "4"], c(.25, .5, .75))
  qPSd4 <- BlomEstimator(waeSD$weight[waeSD$psd=="Q-P" & waeSD$lake == "4"], c(.25, .5, .75))
  pMSd4 <- BlomEstimator(waeSD$weight[waeSD$psd=="P-M" & waeSD$lake == "4"], c(.25, .5, .75))
  mTSd4 <- BlomEstimator(waeSD$weight[waeSD$psd=="M-T" & waeSD$lake == "4"], c(.25, .5, .75))
  gTSd4 <- c(NA, NA, NA)
  
  #SD lake 13
  sSSd13 <- BlomEstimator(waeSD$weight[waeSD$psd=="SS" & waeSD$lake == "13"], c(.25, .5, .75))
  sQSd13 <- BlomEstimator(waeSD$weight[waeSD$psd=="S-Q" & waeSD$lake == "13"], c(.25, .5, .75))
  qPSd13 <- BlomEstimator(waeSD$weight[waeSD$psd=="Q-P" & waeSD$lake == "13"], c(.25, .5, .75))
  pMSd13 <- BlomEstimator(waeSD$weight[waeSD$psd=="P-M" & waeSD$lake == "13"], c(.25, .5, .75))
  mTSd13 <- BlomEstimator(waeSD$weight[waeSD$psd=="M-T" & waeSD$lake == "13"], c(.25, .5, .75))
  gTSd13 <- c(NA, NA, NA)
 
  #SD lake 25
  sSSd25 <- BlomEstimator(waeSD$weight[waeSD$psd=="SS" & waeSD$lake == "25"], c(.25, .5, .75))
  sQSd25 <- BlomEstimator(waeSD$weight[waeSD$psd=="S-Q" & waeSD$lake == "25"], c(.25, .5, .75))
  qPSd25 <- BlomEstimator(waeSD$weight[waeSD$psd=="Q-P" & waeSD$lake == "25"], c(.25, .5, .75))
  pMSd25 <- BlomEstimator(waeSD$weight[waeSD$psd=="P-M" & waeSD$lake == "25"], c(.25, .5, .75))
  mTSd25 <- BlomEstimator(waeSD$weight[waeSD$psd=="M-T" & waeSD$lake == "25"], c(.25, .5, .75))
  gTSd25 <- c(NA, NA, NA)
 
  #Create a table of all Blom estimators, .25, .5, and .75 from the Reference 
  #populations and the three GA and SD pops
  quarTable <- matrix(c(sSRef, sQRef, qPRef, pMRef, mTRef, gTRef, 
                        sSGa2, sQGa2, qPGa2, pMGa2, mTGa2, gTGa2, 
                        sSGa3, sQGa3, qPGa3, pMGa3, mTGa3, gTGa3, 
                        sSGa4, sQGa4, qPGa4, pMGa4, mTGa4, gTGa4, 
                        sSSd4, sQSd4, qPSd4, pMSd4, mTSd4, gTSd4, 
                        sSSd13, sQSd13, qPSd13, pMSd13, mTSd13, gTSd13, 
                        sSSd25, sQSd25, qPSd25, pMSd25, mTSd25, gTSd25), 
                        nrow=7, ncol=18, dimnames=list(c("Reference", 
                                                         "GA 1", "GA 2", "GA 3", 
                                                         "SD 1", "SD 2", "SD 3"), 
                                                       c(".25",".5", ".75",
                                                         ".25",".5", ".75",
                                                         ".25",".5", ".75",
                                                         ".25",".5", ".75",
                                                         ".25",".5", ".75",
                                                         ".25",".5", ".75")),
                        byrow=T)

  #Set new output wd() to write the quarTable
#  setwd("C:/Users/Steven Harris Ranney/Desktop/Manuscripts/Ws - Quantiles/output")
  
  #Write the quarTable to the directory
#  write.csv(quarTable, paste0(getwd(), "/quarTable.csv"), row.names=T)
  
  #-----------------------------------------------------------------------------
  #Estimate differences between Reference slope and intercept and other models

  #Bootstrap the sampling distributions of slope and intercept values by 
  #running n regressions with numSamples drawn randomly with replacement from 
  #a given dataset
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

  #Mann-Whitney U test to evaluate for differences between 
  refDist <- EstLinQuantRegCoefDist(wae, 10000, 100, tau=0.75)

  ga2Dist <- EstLinQuantRegCoefDist(waeGA[waeGA$lake == "2",], 10000, 100, 0.75)
#  hist(ga2Dist[,2], breaks=100, col="gray")
#  hist(refDist[,2], breaks=100, add=T, col=MakeColorsTransparent("red"))
#  wilcox.test(ga2Dist[,2], refDist[,2])
#  wilcox.test(ga2Dist[,1], refDist[,1])
 
  ga3Dist <- EstLinQuantRegCoefDist(waeGA[waeGA$lake == "3",], 10000, 100, 0.75)
#  hist(ga3Dist[,2], breaks=100, col="gray")
#  hist(refDist[,2], breaks=100, add=T, col=MakeColorsTransparent("red"))
#  wilcox.test(ga3Dist[,2], refDist[,2])
#  wilcox.test(ga3Dist[,1], refDist[,1])
   
  ga4Dist <- EstLinQuantRegCoefDist(waeGA[waeGA$lake == "4",], 10000, 100, 0.75)
#  hist(ga4Dist[,2], breaks=100, col="gray")
#  hist(refDist[,2], breaks=100, add=T, col=MakeColorsTransparent("red"))
#  wilcox.test(ga4Dist[,2], refDist[,2])
#  wilcox.test(ga4Dist[,1], refDist[,1])
  
  sd4Dist <- EstLinQuantRegCoefDist(waeSD[waeSD$lake == "4",], 10000, 100, 0.75)
#  hist(sd4Dist[,2], breaks=100, col="gray")
#  hist(refDist[,2], breaks=100, add=T, col=MakeColorsTransparent("red"))
#  wilcox.test(sd4Dist[,2], refDist[,2])
#  wilcox.test(sd4Dist[,1], refDist[,1])
 
  sd13Dist <- EstLinQuantRegCoefDist(waeSD[waeSD$lake == "13",], 10000, 100, 0.75)
#  hist(sd13Dist[,2], breaks=100, col="gray")
#  hist(refDist[,2], breaks=100, add=T, col=MakeColorsTransparent("red"))
#  wilcox.test(sd13Dist[,2], refDist[,2])
#  wilcox.test(sd13Dist[,1], refDist[,1])
 
  sd25Dist <- EstLinQuantRegCoefDist(waeSD[waeSD$lake == "25",], 10000, 100, 0.75)
#  hist(sd25Dist[,2], breaks=100, col="gray")
#  hist(refDist[,2], breaks=100, add=T, col=MakeColorsTransparent("red"))
#  wilcox.test(sd25Dist[,2], refDist[,2])
#  wilcox.test(sd25Dist[,1], refDist[,1])
   
  #Create a table of slope coefficients and bootstrapped confidence intervals
   slopeVals <- as.data.frame(matrix(data = NA, ncol=3, nrow=7, byrow=T, 
                             dimnames = list(c("Reference", "GA 1", 
                                               "GA 2", "GA 3", "SD 1",
                                               "SD 2", "SD 3"), 
                                             c("2.5", "mean", "97.5"))))
                                             
  meanVals <- c(mean(refDist[,2]), mean(ga2Dist[,2]), mean(ga3Dist[,2]), 
                mean(ga4Dist[,2]), mean(sd4Dist[,2]), mean(sd13Dist[,2]), 
                mean(sd25Dist[,2]))
  lowerVals <- c(BlomEstimator(refDist[,2], 0.025), 
                 BlomEstimator(ga2Dist[,2], 0.025),
                 BlomEstimator(ga3Dist[,2], 0.025),
                 BlomEstimator(ga4Dist[,2], 0.025),
                 BlomEstimator(sd4Dist[,2], 0.025),
                 BlomEstimator(sd13Dist[,2], 0.025),
                 BlomEstimator(sd25Dist[,2], 0.025))
  upperVals <- c(BlomEstimator(refDist[,2], 0.975), 
                 BlomEstimator(ga2Dist[,2], 0.975),
                 BlomEstimator(ga3Dist[,2], 0.975),
                 BlomEstimator(ga4Dist[,2], 0.975),
                 BlomEstimator(sd4Dist[,2], 0.975),
                 BlomEstimator(sd13Dist[,2], 0.975),
                 BlomEstimator(sd25Dist[,2], 0.975))
                 
  slopeVals[,1] <- lowerVals
  slopeVals[,2] <- meanVals
  slopeVals[,3] <- upperVals
                 
                 
#  par(mfrow=c(1,1))
  #Plot slope and bootstrapped 95% confidence intervals on the same plot for all 
  #population data
  require(Hmisc)
  plot(slopeVals[,2], pch=19, ylim=c(2.9, 3.6), yaxs='i', bty="n", ylab="Slope",
       xlab=NA, xaxt="n", mar=c(0,4,4,2))
  abline(h=2.9)
  errbar(1:7, slopeVals[,2], yplus=slopeVals[,3], yminus=slopeVals[,1], add=T)
  axis(1, 1:7, labels=c("Reference", "GA 1", "GA 2", "GA 3", "SD 1", "SD 2", "SD 3"))
  abline(h=c(slopeVals[1,]), lty=c(2,1,2))               

  #Create a table of intercept coefficients and bootstrapped confidence intervals
  intVals <- as.data.frame(matrix(data = NA, ncol=3, nrow=7, byrow=T, 
                             dimnames = list(c("Reference", "GA 1", 
                                               "GA 2", "GA 3", "SD 1",
                                               "SD 2", "SD 3"), 
                                             c("2.5", "mean", "97.5"))))
                                             
  meanVals <- c(mean(refDist[,1]), mean(ga2Dist[,1]), mean(ga3Dist[,1]), 
                mean(ga4Dist[,1]), mean(sd4Dist[,1]), mean(sd13Dist[,1]), 
                mean(sd25Dist[,1]))
  lowerVals <- c(BlomEstimator(refDist[,1], 0.025), 
                 BlomEstimator(ga2Dist[,1], 0.025),
                 BlomEstimator(ga3Dist[,1], 0.025),
                 BlomEstimator(ga4Dist[,1], 0.025),
                 BlomEstimator(sd4Dist[,1], 0.025),
                 BlomEstimator(sd13Dist[,1], 0.025),
                 BlomEstimator(sd25Dist[,1], 0.025))
  upperVals <- c(BlomEstimator(refDist[,1], 0.975), 
                 BlomEstimator(ga2Dist[,1], 0.975),
                 BlomEstimator(ga3Dist[,1], 0.975),
                 BlomEstimator(ga4Dist[,1], 0.975),
                 BlomEstimator(sd4Dist[,1], 0.975),
                 BlomEstimator(sd13Dist[,1], 0.975),
                 BlomEstimator(sd25Dist[,1], 0.975))
                 
  intVals[,1] <- lowerVals
  intVals[,2] <- meanVals
  intVals[,3] <- upperVals
                 
                 
  #Plot slope and bootstrapped 95% confidence intervals on the same plot for all 
  #population data
#  require(Hmisc)
  plot(intVals[,2], pch=19, ylim=c(-6.5, -4.5), yaxs = "i", bty="n", 
       ylab="Intercept", xlab="Population", xaxt="n", mar=c(5,4,0,2))
  abline(h=-6.5)
  errbar(1:7, intVals[,2], yplus=intVals[,3], yminus=intVals[,1], add=T)
  axis(1, 1:7, labels=c("Reference", "GA 1", "GA 2", "GA 3", "SD 1", "SD 2", "SD 3"))
  abline(h=c(intVals[1,]), lty=c(2,1,2))               
  par(mfrow=c(1,1))
  
  
  #Combine the slope and intercept values into one table
  regVals <- 
    cbind(intVals %>%
            rename(`int2.5` = `2.5`, 
                   intmean = mean, 
                   `int97.5` = `97.5`), 
          slopeVals %>%
            rename(`s2.5` = `2.5`, 
                   smean = mean, 
                   `s97.5` = `97.5`))
#  write.csv(regVals, paste0("C:/Users/Steven Harris Ranney/Desktop/Manuscripts/Ws - Quantiles/", 
#                      "/output/regVals.csv"))
#  write.csv(wae.pred.values, paste0("C:/Users/Steven Harris Ranney/Desktop/Manuscripts/Ws - Quantiles/", 
#                      "/output/wae.pred.values.csv"))
  
  # Do the 95% estimates of slope for each population overlap the slope value for the reference population?
  regVals %>%
    mutate(row_num = row_number(),
           pop = row.names(regVals),
           is_slope_overlap = ifelse(smean > regVals$s2.5[1] &
                                 smean < regVals$s97.5[1], TRUE, FALSE))
    
  

  #Plot all slope and intercept values on the same plot?  With horizontal and 
  #vertical error bars?
  xlab <- bquote(.("Slope" ~ (beta[1])))
  ylab <- bquote(.("Intercept" ~ (beta[0])))

  plot(regVals[1,2]~regVals[1,5], ylim=c(-6.5, -4.5), xlim=c(2.9, 3.6), yaxs="i", 
       xaxs="i", ylab=ylab, xlab=xlab, bty="n", pch=8)
    points(regVals[2,2]~regVals[2,5], pch=0)
    points(regVals[3,2]~regVals[3,5], pch=2)  
    points(regVals[4,2]~regVals[4,5], pch=5)
    points(regVals[5,2]~regVals[5,5], pch=15)
    points(regVals[6,2]~regVals[6,5], pch=17)  
    points(regVals[7,2]~regVals[7,5], pch=18)
  plotCI(regVals[,5], regVals[,2], uiw = 0.5, ui=regVals[,3], li=regVals[,1], err="y", 
         add=T, pch=NA, gap=0.3)
  plotCI(regVals[,5], regVals[,2], uiw = 0.5, ui=regVals[,6], li=regVals[,4], err="x", 
         add=T, pch=NA, gap=0.3)
  legend("topright", pch=c(8,0,2,5,15,17,18), legend=c("Reference", "GA 1", "GA 2", 
                                                       "GA 3", "SD 1", "SD 2", "SD 3"),
         bty="n")
         
  #IF REFERENCE POPULATION SLOPE AND INTERCEPT 95% CI DO NOT OVERLAP WITH POPULATION
  #LEVEL SLOPE AND INTERCEPT, DETERMINE WHETHER THE INTERVAL AROUND THE DIFFERENCE
  #IN SLOPES CONTAINS ZERO; FROM POPE AND KRUSE 2007 AIFFD, page 433
  
  (mean(refDist[,2])-mean(sd25Dist[,2]))+(1.96*sqrt(((sd(refDist[,2])/(sqrt(nrow(refDist))))^2+(sd(refDist[,2])/(sqrt(nrow(refDist))))^2)))
  (mean(refDist[,2])-mean(sd25Dist[,2]))-(1.96*sqrt(((sd(refDist[,2])/(sqrt(nrow(refDist))))^2+(sd(refDist[,2])/(sqrt(nrow(refDist))))^2)))

  
  std_error <- function(x){
    sd(x)/sqrt(length(x))
  }
  

  #Fx to return the interval around the differences in slopes
  interval_around_differences <- function(population, reference){
    
    upper <- (mean(reference$slope) - mean(population$slope)) + (1.96*sqrt(std_error(reference$slope)^2 + std_error(population$slope)^2))
    lower <- (mean(reference$slope) - mean(population$slope)) - (1.96*sqrt(std_error(reference$slope)^2 + std_error(population$slope)^2))
    
    data.frame(upper = upper, 
               lower = lower) %>%
      return()
  }
  
  list(ga2Dist, ga3Dist, ga4Dist, 
    sd4Dist, sd13Dist, sd25Dist) %>%
    map_dfr(interval_around_differences, refDist) %>%
    mutate(pop = c("GA1", "GA2", "GA3", "SD1", "SD2", "SD3"))
  

    #Create table of 



  #BELOW IS NOW SUPERFLUOUS
  #-------------------------------------------------------------------------------
  #-------------------------------------------------------------------------------
  #-------------------------------------------------------------------------------
  #-------------------------------------------------------------------------------
  #-------------------------------------------------------------------------------
  #-------------------------------------------------------------------------------
  #with simulations?
  #WROTE A FUNCTION TO DO THE NEXT ~90 LINES; FUNCTION AND ANALYSES STARTS ON
  #LINE 503
  
  #Perhaps test for equivalence between the slope and intercept values?  
  #require(equivalence)
  #tost()
  
  #Bootstrap 10,000 slope values from the Q3 regressions based upon the Koenker
  #2004 95% confidence intervals and compare the values with a t-test.
#  tmpGa2 <- rnorm(10000, summary(ga2, se="rank", alpha=0.05)$coefficients[2,1], 
#              (summary(ga2, se="rank", alpha=.05)$coefficients[2,1]-summary(ga2, se="rank", alpha=.05)$coefficients[2,2])/1.96)
#  tmpRef <- rnorm(10000, summary(wae.linMod.75, se="rank", alpha=0.05)$coefficients[2,1], 
#              (summary(wae.linMod.75, se="rank", alpha=.05)$coefficients[2,1]-summary(wae.linMod.75, se="rank", alpha=.05)$coefficients[2,2])/1.96)
#  t.test(tmpRef, tmpGa2)
#  tost(tmpGa2, tmpRef, epsilon=.001)

#  tmpGa3 <- rnorm(10000, summary(ga3, se="rank", alpha=0.05)$coefficients[2,1], 
#              (summary(ga3, se="rank", alpha=.05)$coefficients[2,1]-summary(ga3, se="rank", alpha=.05)$coefficients[2,2])/1.96)
#  tmpRef <- rnorm(10000, summary(wae.linMod.75, se="rank", alpha=0.05)$coefficients[2,1], 
#              (summary(wae.linMod.75, se="rank", alpha=.05)$coefficients[2,1]-summary(wae.linMod.75, se="rank", alpha=.05)$coefficients[2,2])/1.96)
#  t.test(tmpRef, tmpGa3)
  
#  tmpGa4 <- rnorm(10000, summary(ga4, se="rank", alpha=0.05)$coefficients[2,1], 
#              (summary(ga4, se="rank", alpha=.05)$coefficients[2,1]-summary(ga4, se="rank", alpha=.05)$coefficients[2,2])/1.96)
#  tmpRef <- rnorm(10000, summary(wae.linMod.75, se="rank", alpha=0.05)$coefficients[2,1], 
#              (summary(wae.linMod.75, se="rank", alpha=.05)$coefficients[2,1]-summary(wae.linMod.75, se="rank", alpha=.05)$coefficients[2,2])/1.96)
#  t.test(tmpRef, tmpGa4)

  #NOW WITH SD POPULATION SLOPES
#  tmpSd4 <- rnorm(10000, summary(sd4, se="rank", alpha=0.05)$coefficients[2,1], 
#              (summary(sd4, se="rank", alpha=.05)$coefficients[2,1]-summary(sd4, se="rank", alpha=.05)$coefficients[2,2])/1.96)
#  tmpRef <- rnorm(10000, summary(wae.linMod.75, se="rank", alpha=0.05)$coefficients[2,1], 
#              (summary(wae.linMod.75, se="rank", alpha=.05)$coefficients[2,1]-summary(wae.linMod.75, se="rank", alpha=.05)$coefficients[2,2])/1.96)
#  t.test(tmpRef, tmpSd4)

#  tmpSd13 <- rnorm(10000, summary(sd13, se="rank", alpha=0.05)$coefficients[2,1], 
#              (summary(sd13, se="rank", alpha=.05)$coefficients[2,1]-summary(sd13, se="rank", alpha=.05)$coefficients[2,2])/1.96)
#  tmpRef <- rnorm(10000, summary(wae.linMod.75, se="rank", alpha=0.05)$coefficients[2,1], 
#              (summary(wae.linMod.75, se="rank", alpha=.05)$coefficients[2,1]-summary(wae.linMod.75, se="rank", alpha=.05)$coefficients[2,2])/1.96)
#  t.test(tmpRef, tmpSd13)
  
#  tmpSd25 <- rnorm(10000, summary(sd25, se="rank", alpha=0.05)$coefficients[2,1], 
#              (summary(sd25, se="rank", alpha=.05)$coefficients[2,1]-summary(sd25, se="rank", alpha=.05)$coefficients[2,2])/1.96)
#  tmpRef <- rnorm(10000, summary(wae.linMod.75, se="rank", alpha=0.05)$coefficients[2,1], 
#              (summary(wae.linMod.75, se="rank", alpha=.05)$coefficients[2,1]-summary(wae.linMod.75, se="rank", alpha=.05)$coefficients[2,2])/1.96)
#  t.test(tmpRef, tmpSd25)

#####################  
  #Bootstrap 10,000 INTERCEPT values from the Q3 regressions based upon the Koenker
  #2004 95% confidence intervals and compare the values with a t-test.
#  tmpGa2 <- rnorm(10000, summary(ga2, se="rank", alpha=0.05)$coefficients[1,1], 
#              (summary(ga2, se="rank", alpha=.05)$coefficients[1,1]-summary(ga2, se="rank", alpha=.05)$coefficients[1,2])/1.96)
#  tmpRef <- rnorm(10000, summary(wae.linMod.75, se="rank", alpha=0.05)$coefficients[1,1], 
#              (summary(wae.linMod.75, se="rank", alpha=.05)$coefficients[1,1]-summary(wae.linMod.75, se="rank", alpha=.05)$coefficients[1,2])/1.96)
#  t.test(tmpRef, tmpGa2)

#  tmpGa3 <- rnorm(10000, summary(ga3, se="rank", alpha=0.05)$coefficients[1,1], 
#              (summary(ga3, se="rank", alpha=.05)$coefficients[1,1]-summary(ga3, se="rank", alpha=.05)$coefficients[1,2])/1.96)
#  tmpRef <- rnorm(10000, summary(wae.linMod.75, se="rank", alpha=0.05)$coefficients[1,1], 
#              (summary(wae.linMod.75, se="rank", alpha=.05)$coefficients[1,1]-summary(wae.linMod.75, se="rank", alpha=.05)$coefficients[1,2])/1.96)
#  t.test(tmpRef, tmpGa3)
  
#  tmpGa4 <- rnorm(10000, summary(ga4, se="rank", alpha=0.05)$coefficients[1,1], 
#              (summary(ga4, se="rank", alpha=.05)$coefficients[1,1]-summary(ga4, se="rank", alpha=.05)$coefficients[1,2])/1.96)
#  tmpRef <- rnorm(10000, summary(wae.linMod.75, se="rank", alpha=0.05)$coefficients[1,1], 
#              (summary(wae.linMod.75, se="rank", alpha=.05)$coefficients[1,1]-summary(wae.linMod.75, se="rank", alpha=.05)$coefficients[1,2])/1.96)
#  t.test(tmpRef, tmpGa4)

  #NOW WITH SD POPULATION INTERCEPTS
#  tmpSd4 <- rnorm(10000, summary(sd4, se="rank", alpha=0.05)$coefficients[1,1], 
#              (summary(sd4, se="rank", alpha=.05)$coefficients[1,1]-summary(sd4, se="rank", alpha=.05)$coefficients[1,2])/1.96)
#  tmpRef <- rnorm(10000, summary(wae.linMod.75, se="rank", alpha=0.05)$coefficients[1,1], 
#              (summary(wae.linMod.75, se="rank", alpha=.05)$coefficients[1,1]-summary(wae.linMod.75, se="rank", alpha=.05)$coefficients[1,2])/1.96)
#  t.test(tmpRef, tmpSd4)

#  tmpSd13 <- rnorm(10000, summary(sd13, se="rank", alpha=0.05)$coefficients[1,1], 
#              (summary(sd13, se="rank", alpha=.05)$coefficients[1,1]-summary(sd13, se="rank", alpha=.05)$coefficients[1,2])/1.96)
#  tmpRef <- rnorm(10000, summary(wae.linMod.75, se="rank", alpha=0.05)$coefficients[1,1], 
#              (summary(wae.linMod.75, se="rank", alpha=.05)$coefficients[1,1]-summary(wae.linMod.75, se="rank", alpha=.05)$coefficients[1,2])/1.96)
#  t.test(tmpRef, tmpSd13)
  
#  tmpSd25 <- rnorm(10000, summary(sd25, se="rank", alpha=0.05)$coefficients[1,1], 
#              (summary(sd25, se="rank", alpha=.05)$coefficients[1,1]-summary(sd25, se="rank", alpha=.05)$coefficients[1,2])/1.96)
#  tmpRef <- rnorm(10000, summary(wae.linMod.75, se="rank", alpha=0.05)$coefficients[1,1], 
#              (summary(wae.linMod.75, se="rank", alpha=.05)$coefficients[1,1]-summary(wae.linMod.75, se="rank", alpha=.05)$coefficients[1,2])/1.96)
#  t.test(tmpRef, tmpSd25)
  
  #Equivalence testing
#  tost(tmpSd25, tmpRef, epsilon=0.00001)
  
  #Function to take parameter estimates from model summaries, calculate the SD
  #from the 95% confidence intervals, simulate a distribution of data, and 
  # compare the vectors of numbers.  This function returns the p-value from a 
  #t-test on the two vectors
#  TestForDifferencesInParameters <- function(mod1, mod2, method="rank", alpha=0.05, n, 
#                                        parameter="slope", R=200, bsmethod="xy"){
#      if(parameter == "slope"){
#      param <- 2
#      } else {
#      param <- 1 
#      }
#      #Create a distribution of values from the normal distribution with n samples, 
#      #and mean from the parameter estimate of interest, and the SD calculated from
#      #the 95% confidence intervals
#      if(method == "rank"){
#      tmpMod1 <- rnorm(n, summary(mod1, se=method, alpha=alpha)$coefficients[param, 1],
#                       (summary(mod1, se=method, alpha=alpha)$coefficients[param, 1]-summary(mod1, se=method, alpha=alpha)$coefficients[param, 2])/1.96)
#      tmpMod2 <- rnorm(n, summary(mod2, se=method, alpha=alpha)$coefficients[param, 1],
#                       (summary(mod2, se=method, alpha=alpha)$coefficients[param, 1]-summary(mod2, se=method, alpha=alpha)$coefficients[param, 2])/1.96)
#      }
#      if(method != "rank"){
#      tmpMod1 <- rnorm(n, summary(mod1, se=method, alpha=alpha, R, bsmethod)$coefficients[param, 1],
#                       (summary(mod1, se=method, alpha=alpha)$coefficients[param, 2]*sqrt(nrow(mod1$model))))
#      tmpMod2 <- rnorm(n, summary(mod2, se=method, alpha=alpha)$coefficients[param, 1],
#                       (summary(mod2, se=method, alpha=alpha)$coefficients[param, 2]*sqrt(nrow(mod2$model))))
#      }
#      return(t.test(tmpMod1, tmpMod2)$p.value)
#  }

   #Test of system time to run these summaries with bootstrapped standard errors
#  system.time(tmp1 <- summary(ga2, se="boot", alpha=0.05, R=200, bsmethod="xy"))
#  system.time(tmp2 <- summary(wae.linMod.75, se="boot", alpha=0.05, R=200, bsmethod="xy"))
  
  
#  ga1S <- TestForDifferencesInParameters(ga2, wae.linMod.75, method="boot", n=10000)
#  ga2S <- TestForDifferencesInParameters(ga3, wae.linMod.75, method="boot", n=10000)
#  ga3S <- TestForDifferencesInParameters(ga4, wae.linMod.75, method="boot", n=10000)
#  ga1I <- TestForDifferencesInParameters(ga2, wae.linMod.75, method="boot", n=10000, parameter="intercept")
#  ga2I <- TestForDifferencesInParameters(ga3, wae.linMod.75, method="boot", n=10000, parameter="intercept")
#  ga3I <- TestForDifferencesInParameters(ga4, wae.linMod.75, method="boot", n=10000, parameter="intercept")

#  sd1S <- TestForDifferencesInParameters(sd4, wae.linMod.75, method="boot", n=10000)
#  sd2S <- TestForDifferencesInParameters(sd13, wae.linMod.75, method="boot", n=10000)
#  sd3S <- TestForDifferencesInParameters(sd25, wae.linMod.75, method="boot", n=10000)
#  sd1I <- TestForDifferencesInParameters(sd4, wae.linMod.75, method="boot", n=10000, parameter="intercept")
#  sd2I <- TestForDifferencesInParameters(sd13, wae.linMod.75, method="boot", n=10000, parameter="intercept")
#  sd3I <- TestForDifferencesInParameters(sd25, wae.linMod.75, method="boot", n=10000, parameter="intercept")
 
#  slopeIntPvalues <- as.data.frame(matrix(NA, nrow=7, ncol=2, dimnames=list(c("Reference", "Georgia 1", 
#                                                         "Georgia 2", "Georgia 3", 
#                                                         "South Dakota 1", "South Dakota 2", 
#                                                         "South Dakota 3"), 
#                                                     c("Beta0Pvalue", "Beta1Pvalue"))))
#  
#  slopeIntPvalues[2,1] <- ga1I
#  slopeIntPvalues[3,1] <- ga2I
#  slopeIntPvalues[4,1] <- ga3I
#  slopeIntPvalues[5,1] <- sd1I
#  slopeIntPvalues[6,1] <- sd2I
#  slopeIntPvalues[7,1] <- sd3I
  
#  slopeIntPvalues[2,2] <- ga1S
#  slopeIntPvalues[3,2] <- ga2S
#  slopeIntPvalues[4,2] <- ga3S
#  slopeIntPvalues[5,2] <- sd1S
#  slopeIntPvalues[6,2] <- sd2S
#  slopeIntPvalues[7,2] <- sd3S

#  slopeInt <- cbind(slopeInt, slopeIntPvalues)
 
#  slopeInt <- slopeInt[, c("Beta0", "Beta0Pvalue", "Beta1", "Beta1Pvalue")]
 
  #Wrapper to test for equivalence
  #Null hypothesis is that the parameters are different
  #Alternative hypotheis is that paramaters are similar
#  TestForEquivalence <- function(mod1, mod2, parameter="slope", n, alpha = .05, epsilon){
#    if(parameter == "slope"){
#    param <- 2
#    } else {
#    param <- 1 
#    }
    #Create a distribution of values from the normal distribution with n samples, 
    #and mean from the parameter estimate of interest, and the SD calculated from
    #the 95% confidence intervals
#    tmpMod1 <- rnorm(n, summary(mod1, se="rank", alpha=alpha)$coefficients[param, 1],
#                     (summary(mod1, se="rank", alpha=alpha)$coefficients[param, 1]-summary(mod1, se="rank", alpha=alpha)$coefficients[param, 2])/1.96)
#    tmpMod2 <- rnorm(n, summary(mod2, se="rank", alpha=alpha)$coefficients[param, 1],
#                     (summary(mod2, se="rank", alpha=alpha)$coefficients[param, 1]-summary(mod2, se="rank", alpha=alpha)$coefficients[param, 2])/1.96)
#    return(tost(tmpMod1, tmpMod2, epsilon = epsilon)$p.value)
#  } 

#  ga1S <- TestForEquivalence(ga2, wae.linMod.75, n=10000)
#  ga2S <- TestForEquivalence(ga3, wae.linMod.75, n=10000)
#  ga3S <- TestForEquivalence(ga4, wae.linMod.75, n=10000)
#  ga1I <- TestForEquivalence(ga2, wae.linMod.75, n=10000, parameter="intercept")
#  ga2I <- TestForEquivalence(ga3, wae.linMod.75, n=10000, parameter="intercept")
#  ga3I <- TestForEquivalence(ga4, wae.linMod.75, n=10000, parameter="intercept")

#  sd1S <- TestForEquivalence(sd4, wae.linMod.75, n=10000)
#  sd2S <- TestForEquivalence(sd13, wae.linMod.75, n=10000)
#  sd3S <- TestForEquivalence(sd25, wae.linMod.75, n=10000)
#  sd1I <- TestForEquivalence(sd4, wae.linMod.75, n=10000, parameter="intercept")
#  sd2I <- TestForEquivalence(sd13, wae.linMod.75, n=10000, parameter="intercept")
#  sd3I <- TestForEquivalence(sd25, wae.linMod.75, n=10000, parameter="intercept")
    
  
  
  #-----------------------------------------------------------------------------
  ###BELOW IS FROM BRIAN CADE AT USGS IN FORT COLLINS
  ###First 3 commands creates vectors of rho by 0.05 to 0.95 quantiles 
  ###for models with separate slopes and intercepts by group, 
  ###common slopes by groups, and with just separate intercepts (null). 
  ###Next 2 commands perform the R1 computations, then those are combined 
  ###into common object for matrix plotting of results. 

#  rqR1.separate <- (rq(log10(weight) ~ log10(length),data=wae,tau=5:95/100))$rho 

#  rqR1.common <- (rq(log10(weight) ~ log10(length),data=wae,tau=5:95/100))$rho 

#  rqR1.intercept <- (rq(log10(weight) ~ 1,data=wae,tau=5:95/100))$rho 

#  rqR1.separate.value <- 1 - rqR1.separate/rqR1.intercept 
#  rqR1.common.value <- 1 - rqR1.common/rqR1.intercept 
#  rqR1.both.value <- cbind(rqR1.separate.value,rqR1.common.value) 

#  rqR1.tau <- c(5:95/100) 
#  matplot(rqR1.tau,rqR1.both.value,type="o",cex=0.75,ylim=c(0.85,0.88))   
 
#  wae.mod.1 <- nlrq(weight~alpha*length^beta, start=list(alpha=c(0.000001, 0.000001, 0.000001), 
#                    beta=c(3.33, 3.41, 3.45)), tau=c(0.25,0.5,0.75), data=wae)
#  plot(summary(wae.mod.1, se="rank", iid=F, alpha=0.05), pch=19)
#  wae.mod <- nlrq(log10(weight)~log10(length), data=wae, tau=c(25,50,75)/100,                         #Requires different parm. ests. for alpha and beta
#                 start=list(alpha=0.00001, beta=3), trace=T)
#  Plot regression summaries at once
#  WARNING: THIS WILL TAKE SOME TIME, 10 minutes for the model here
#  system.time(plot(summary(wae.mod,se="rank",iid=F,alpha=0.05)))
#  plot(summary(wae.mod,se="rank",iid=F,alpha=0.05))
 