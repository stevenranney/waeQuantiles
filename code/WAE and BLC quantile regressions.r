################################################################################
#This program creates non-linear quantile regression for tau = 0.75, 0.5, and  #
#0.25 for both WAE and BLC from the CLEAN datasets used in Ranney et al.       #
#(2010).  Weights (by 10mm length class) are predicted from each quantile      #
#regression and tabulated.  A TABLE is then created that provides intercept    #
#(a) and slope (b) for each regression.  Ultimately, this information will be  #
#presented in a manuscript for publication in NAJFM that will offer a novel    #
#method with which to compare fish condition, based on some previous work done #
#by Brian Cade USGS in Fort Collins.                                           #                 
################################################################################ 

setwd('C:/Users/Steven Harris Ranney/Desktop/Manuscripts/Ws - Quantiles/code')

wae <- read.table('WAE Clean.txt', header=T)
attach(wae)
names(wae)

require(quantreg)
require(Hmisc)

#Assigns fish to a length category (Gabelhouse 1984)
psd <- with(wae,
            ifelse((length>=250)&(length<380), 250,
            ifelse((length>=380)&(length<510), 380,
            ifelse((length>=510)&(length<630), 510,
            ifelse((length>=630)&(length<760), 630,
            ifelse(length>=760, 760,
            0))))))

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



#Build a matrix of cutoff values for each lake
#cutoff <- matrix(NA,101,1, dimnames=list(c(1:101), c("abs(dffits) cutoff value")))
#mod <- c()

#for (i in 1:length(unique(lake))){
#     mod=lm(log10(weight[lake==unique(lake)[i]])~log10(length[lake==unique(lake)[i]]))
#     cutoff[i,1]=as.vector(sqrt(2/length(length[lake==unique(lake)[i]])))
#     int.slope[i,2]=as.vector(mod$coef[2])
#     int.slope[i,3]=as.vector(summary(mod)$coef[1,4])
#}


#wae.lm <- lm(log10(Weight)~log10(Length))

#cutoff <- sqrt(2/length(Length))

#length.1 <- ifelse((abs(dffits(wae.lm))>cutoff), NA, 
#            ifelse((Length<160), NA,
#            ifelse((Length>740), NA, Length)))

wae.mod.75 <- nlrq(weight~alpha*length^beta, data=wae, tau=0.75,                         #Requires different parm. ests. for alpha and beta
      start=list(alpha=0.00001, beta=3), trace=T)
wae.mod.5 <- nlrq(weight~alpha*length^beta, data=wae, tau=0.5,                           #Requires different parm. ests. for alpha and beta
      start=list(alpha=0.00001, beta=3), trace=T)
wae.mod.25 <- nlrq(weight~alpha*length^beta, data=wae, tau=0.25,                         #Requires different parm. ests. for alpha and beta
      start=list(alpha=0.00001, beta=3), trace=T)
#wae.mod.0 <- nlrq(weight~alpha*length^beta, data=wae, tau=0.01,                         #Requires different parm. ests. for alpha and beta
#      start=list(alpha=0.00001, beta=3), trace=T)
#Calculate many regressions at once
#wae.mod <- rq(log10(weight)~log10(length), data=wae, tau=5:90/100)#,                         #Requires different parm. ests. for alpha and beta
#      start=list(alpha=0.00001, beta=3), trace=T)
#Plot regression summaries at once
#WARNING: THIS WILL TAKE SOME TIME, 10 minutes for the model here
#system.time(plot(summary(wae.mod,se="rank",iid=F,alpha=0.05)))
#plot(summary(wae.mod,se="rank",iid=F,alpha=0.05))


#midpoint of each legnth category; double check that max length includes the entirety of the length category 
#x <- as.vector(((tapply(length, psd, max)-tapply(length, psd, min))/2)+tapply(length, psd, min)); x
#x <- c(204.5,314.5,444.5,569.5,695)

x <- seq(165,735, by=10)#; x

#plot data and all lines
plot(weight~length, data=wae, xlab="Length (mm)", xlim=c(100, 800), ylim=c(0,6000), 
      ylab="Weight (g)", pch="", xaxs="i", yaxs="i", bty="n")
#  abline(v=c(250,380,510,630,760), lty=2)
  wae.moda <- seq(150,740)
    lines(wae.moda, predict(wae.mod.75, list(length=wae.moda)), lwd=2, lty=2)
    lines(wae.moda, predict(wae.mod.5, list(length=wae.moda)), lwd=2)
    lines(wae.moda, predict(wae.mod.25, list(length=wae.moda)), lwd=2, lty=2)
#add midpoints of length categories to plots
  points(tapply(weight, l.c.1, mean)~x, pch=1)
#add error bars to points, 95%CI
  errbar(x, 
         tapply(weight, l.c.1, mean), 
         tapply(weight, l.c.1, mean)+(1.96*tapply(weight, l.c.1, sd)), 
         tapply(weight, l.c.1, mean)-(1.96*tapply(weight, l.c.1, sd)), 
         pch=" ", add=T)
#  errbar(x, 
#         tapply(weight, psd, median), 
#         c(quantile(weight[psd==0], 0.75, na.rm=T), quantile(weight[psd==250], 0.75, na.rm=T), quantile(weight[psd==380], 0.75, na.rm=T), quantile(weight[psd==510], 0.75, na.rm=T), quantile(weight[psd==630], 0.75, na.rm=T)), 
#         c(quantile(weight[psd==0], 0.25, na.rm=T), quantile(weight[psd==250], 0.25, na.rm=T), quantile(weight[psd==380], 0.25, na.rm=T), quantile(weight[psd==510], 0.25, na.rm=T), quantile(weight[psd==630], 0.25, na.rm=T)), 
#         pch=" ", add=T)
  text(250, 4000, labels=expression(bold(Walleye)), cex=1.5)
#  text(x, c(5250,5250,5250,5250,5250), labels=c("SS", "S-Q", "Q-P", "P-M", "M-T"))

wae.pred.75 <- as.vector(predict(wae.mod.75, list(length=seq(155,745, by=10)), interval="confidence"))
wae.pred.5 <- as.vector(predict(wae.mod.5, list(length=seq(155,745, by=10)), interval="confidence"))
wae.pred.25 <- as.vector(predict(wae.mod.25, list(length=seq(155,745, by=10)), interval="confidence"))

#creates a matrix 60 rows long by 3 columsn wide
wae.pred.values <- matrix(NA,60,3, dimnames=list(seq(155,745,by=10), c("Q1", "Q2", "Q3")))
#for length categories
#wae.lencat.values <- matrix(NA,5,3, dimnames=list(c(250,380,510,630,760), c("Q1","Q2","Q3")))

#program to put all predicted values from each model into the matrix "pred.values"
for (i in 1:60){
  wae.pred.values[i,1]=as.vector(wae.pred.25[i])
  wae.pred.values[i,2]=as.vector(wae.pred.5[i])
  wae.pred.values[i,3]=as.vector(wae.pred.75[i])
  }
  
wae.pred.values <- as.data.frame(wae.pred.values, header=T)
wae.pred.values 

#for (i in 1:5){
#  wae.lencat.values[i,1]=as.vector(wae.pred.25[i])
#  wae.lencat.values[i,2]=as.vector(wae.pred.5[i])
#  wae.lencat.values[i,3]=as.vector(wae.pred.75[i])
#  }
 
#wae.lencat.values <- as.data.frame(wae.lencat.values, header=T)
#wae.lencat.values 

#with(pred.values, 
#  plot(Q1~seq(205,795,by=10)))

#create a table of intercept and slope values from each nlrq 
wae.int.slope <- matrix(NA,3,2, dimnames=list(c("25th", "50th", "75th"), c("Intercept", "Slope")))
     wae.int.slope[1,1]=summary(wae.mod.25)$coef[1,1]
     wae.int.slope[1,2]=summary(wae.mod.25)$coef[2,1]
     wae.int.slope[2,1]=summary(wae.mod.5)$coef[1,1]
     wae.int.slope[2,2]=summary(wae.mod.5)$coef[2,1]
     wae.int.slope[3,1]=summary(wae.mod.75)$coef[1,1]
     wae.int.slope[3,2]=summary(wae.mod.75)$coef[2,1]

wae.int.slope

#wae.mod.1 <- nlrq(weight~alpha*length^beta, start=list(alpha=c(0.000001, 0.000001, 0.000001), 
#  beta=c(3.33, 3.41, 3.45)), tau=c(0.25,0.5,0.75), data=wae)
#plot(summary(wae.mod.1, se="rank", iid=F, alpha=0.05), pch=19)
#plot(summary(wae.mod.75,se="rank",iid=F,alpha=0.05))

#Plot each length category individually
#par(mfrow=c(5,1))
# plot(weight~length, data=wae, xlim=c(100, 250), ylim=c(0, 200), pch="", xaxs="i", yaxs="i", bty="n", ylab="Weight (g)", xlab="Length (mm)")
#   wae.moda <- seq(160,250)
#     lines(wae.moda, predict(wae.mod.75, list(length=wae.moda)), lwd=2, lty=2)
#     lines(wae.moda, predict(wae.mod.5, list(length=wae.moda)), lwd=2)
#     lines(wae.moda, predict(wae.mod.25, list(length=wae.moda)), lwd=2, lty=2)
#  abline(h=0)
#  abline(v=250, lty=2)
#    points(tapply(weight, l.c.1, mean)~x, pch=1)
#     errbar(x,
#          tapply(weight, l.c.1, mean),
#          tapply(weight, l.c.1, mean)+tapply(weight, l.c.1, sd),
#          tapply(weight, l.c.1, mean)-tapply(weight, l.c.1, sd),
#          pch=" ", add=T)
#  text(140,150, label=expression(bold("Substock")), cex=2)
# plot(weight~length, data=wae, xlim=c(250,380), ylim=c(0,1000), pch="", xaxs="i", yaxs="i", bty="n", ylab="Weight (g)", xlab="Length (mm)")
#   wae.moda <- seq(160,740)
#     lines(wae.moda, predict(wae.mod.75, list(length=wae.moda)), lwd=2, lty=2)
#     lines(wae.moda, predict(wae.mod.5, list(length=wae.moda)), lwd=2)
#     lines(wae.moda, predict(wae.mod.25, list(length=wae.moda)), lwd=2, lty=2)
#    abline(h=0)
#    abline(v=c(250,380), lty=2)
#    points(tapply(weight, l.c.1, mean)~x, pch=1)
#     errbar(x,
#          tapply(weight, l.c.1, mean),
#          tapply(weight, l.c.1, mean)+tapply(weight, l.c.1, sd),
#          tapply(weight, l.c.1, mean)-tapply(weight, l.c.1, sd),
#          pch=" ", add=T)
#     text(300,800, labels=expression(bold("S-Q")), cex=2)
# plot(weight~length, data=wae, xlim=c(380,510), ylim=c(0,2000), pch="", xaxs="i", yaxs="i", bty="n", xpd=F, ylab="Weight (g)", xlab="Length (mm)")
#   wae.moda <- seq(160,545)
#     lines(wae.moda, predict(wae.mod.75, list(length=wae.moda)), lwd=2, lty=2)
#     lines(wae.moda, predict(wae.mod.5, list(length=wae.moda)), lwd=2)
#     lines(wae.moda, predict(wae.mod.25, list(length=wae.moda)), lwd=2, lty=2)
#    abline(v=c(380,510), lty=2)
#    abline(h=0)
#    points(tapply(weight, l.c.1, mean)~x, pch=1)
#     errbar(x,
#          tapply(weight, l.c.1, mean),
#          tapply(weight, l.c.1, mean)+tapply(weight, l.c.1, sd),
#          tapply(weight, l.c.1, mean)-tapply(weight, l.c.1, sd),
#          pch=" ", add=T)
#     text(420,1500, label=expression(bold("Q-P")), cex=2)
# plot(weight~length, data=wae, xlim=c(500,630), ylim=c(0,3500), pch="", xaxs="i", yaxs="i", bty="n", ylab="Weight (g)", xlab="Length (mm)")
#   wae.moda <- seq(160,640)
#     lines(wae.moda, predict(wae.mod.75, list(length=wae.moda)), lwd=2, lty=2)
#     lines(wae.moda, predict(wae.mod.5, list(length=wae.moda)), lwd=2)
#     lines(wae.moda, predict(wae.mod.25, list(length=wae.moda)), lwd=2, lty=2)
#    abline(v=c(510,630), lty=2)
#    abline(h=0)
#    points(tapply(weight, l.c.1, mean)~x, pch=1)
#     errbar(x,
#          tapply(weight, l.c.1, mean),
#          tapply(weight, l.c.1, mean)+tapply(weight, l.c.1, sd),
#          tapply(weight, l.c.1, mean)-tapply(weight, l.c.1, sd),
#          pch=" ", add=T)
#     text(540,3000, label=expression(bold("P-M")), cex=2)
# plot(weight~length, data=wae, xlim=c(620,760), ylim=c(2000,5500), pch="", xaxs="i", yaxs="i", bty="n", ylab="Weight (g)", xlab="Length (mm)")
#   wae.moda <- seq(160,740)
#     lines(wae.moda, predict(wae.mod.75, list(length=wae.moda)), lwd=2, lty=2)
#     lines(wae.moda, predict(wae.mod.5, list(length=wae.moda)), lwd=2)
#     lines(wae.moda, predict(wae.mod.25, list(length=wae.moda)), lwd=2, lty=2)
#    abline(v=c(630,760), lty=2)
#    points(tapply(weight, l.c.1, mean)~x, pch=1)
#     errbar(x,
#          tapply(weight, l.c.1, mean),
#          tapply(weight, l.c.1, mean)+tapply(weight, l.c.1, sd),
#          tapply(weight, l.c.1, mean)-tapply(weight, l.c.1, sd),
#          pch=" ", add=T)
#     text(660,5000, label=expression(bold("P-M")), cex=2)
#


#####################################
#Plot 25th, 50th, and 75th P'tiles on the same page as EmP and RLP lines
x.1 <- seq(160, 740, 1)

#EmP and RLP equations are from Ranney et al. (2010)
emp.y <- ((-4.866)+(2.661*log10(x.1))+((0.111)*(log10(x.1))^2))
rlp.y <- (-5.422)+((3.165)*log10(x.1))

plot(weight~length, data=wae, xlab="Length (mm)", xlim=c(100, 800), ylim=c(0,6000), 
      ylab="Weight (g)", pch="", xaxs="i", yaxs="i", bty="n")
  abline(v=c(250,380,510,630,760), lty=2)
  wae.moda <- seq(150,740)
    lines(wae.moda, predict(wae.mod.75, list(length=wae.moda)), lwd=2, lty=2)
    lines(wae.moda, predict(wae.mod.5, list(length=wae.moda)), lwd=2)
    lines(wae.moda, predict(wae.mod.25, list(length=wae.moda)), lwd=2, lty=2)
  legend(150, 5750, bg="white", legend=c("75th Percentile", "50th Percentile", "25th Percentile", "EmP", "RLP"), lty=c(2,1,2,1,1), lwd=2, col=c("black","black","black","red","blue"))
  text(250, 3500, label=expression(bold("Walleye")), cex=2)

lines((10^emp.y)~x.1, col="red", lwd=2)
lines((10^rlp.y)~x.1, col="blue", lwd=2)

#################################################################################
#compare wae.mod.75 line with independent data
wae.ind <- read.table('WAE_independent.txt', header=T)
attach(wae.ind)
names(wae.ind)

l.c <- with(wae.ind,
  ifelse((length>=80)&(length<90),"085",
  ifelse((length>=90)&(length<100),"095",
  ifelse((length>=100)&(length<110),"105",
  ifelse((length>=110)&(length<120),"115",
  ifelse((length>=120)&(length<130),"125",
  ifelse((length>=130)&(length<140),"135",
  ifelse((length>=140)&(length<150),"145",
  ifelse((length>=150)&(length<160),"155",
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
  NA)))))))))))))))))))))))))))))))))))))))))))))))))
  
l.c.1 <- with(wae.ind,
  ifelse((length>=560)&(length<570),"565",
  ifelse((length>=570)&(length<580),"575",
  ifelse((length>=580)&(length<590),"585",
  ifelse((length>=590)&(length<600),"595",
  ifelse((length>=600)&(length<610),"605",
  ifelse((length>=610)&(length<620),"615",
  ifelse((length>=620)&(length<630),"625",
  ifelse((length>=630)&(length<640),"635",
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
  ifelse((length>=740)&(length<750),"745",
  ifelse((length>=750)&(length<760),"755",
  ifelse((length>=760)&(length<770),"765",
  ifelse((length>=770)&(length<780),"775",
  ifelse((length>=780)&(length<790),"785",
  ifelse((length>=790)&(length<800),"795",
  ifelse((length>=800)&(length<810),"805",
  ifelse((length>=810)&(length<820),"815",
  ifelse((length>=820)&(length<830),"825",
  ifelse((length>=830)&(length<840),"835",
  l.c)))))))))))))))))))))))))))))

x <- seq(85,795,10)

plot(weight~length, data=wae.ind, xlab="Length (mm)", xlim=c(100, 800), ylim=c(0,6000), 
      ylab="Weight (g)", pch="", xaxs="i", yaxs="i", bty="n")
  abline(v=c(250,380,510,630,760), lty=2)
points(tapply(weight, l.c.1, mean)~x, data=wae.ind, pch=19)
  wae.moda <- seq(85,795)
    lines(wae.moda, predict(wae.mod.75, list(length=wae.moda)), lwd=2, lty=2)
  legend(150,5000, bg="white", legend=c("Development data 75th Percentile", "Validation data mean"), lty=c(2,0), lwd=c(2,0), pch=c(NA,19))
  
q3 <- matrix(NA,72,1, dimnames=list(seq(85,795,10), "Q3"), byrow=T)

for (i in 1:length(unique(l.c.1))) {
     q3[i] <- quantile(weight[l.c.1==unique(l.c.1)[i]], 0.75, data=wae.ind, na.rm=T)
    }
q3

z <- c(quantile(weight[l.c.1=="085"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="095"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="105"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="115"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="125"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="135"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="145"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="155"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="165"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="175"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="185"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="195"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="205"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="215"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="225"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="235"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="245"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="255"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="265"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="275"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="285"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="295"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="305"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="315"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="325"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="335"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="345"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="355"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="365"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="375"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="385"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="395"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="405"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="415"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="425"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="435"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="445"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="455"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="465"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="475"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="485"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="495"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="505"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="515"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="525"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="535"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="545"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="555"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="565"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="575"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="585"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="595"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="605"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="615"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="625"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="635"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="645"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="655"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="665"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="675"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="685"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="695"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="705"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="715"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="725"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="735"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="745"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="755"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="765"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="775"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="785"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="795"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="805"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="815"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="825"], 0.75, na.rm=T, data=wae.ind),
       quantile(weight[l.c.1=="835"], 0.75, na.rm=T, data=wae.ind))

x <- seq(85,835,10)

wae.ind.75 <- nlrq(weight~alpha*length^beta, data=wae.ind, tau=0.75,                         #Requires different parm. ests. for alpha and beta
      start=list(alpha=0.00001, beta=3), trace=T)
wae.ind.5 <- nlrq(weight~alpha*length^beta, data=wae.ind, tau=0.5,                         #Requires different parm. ests. for alpha and beta
      start=list(alpha=0.00001, beta=3), trace=T)

wae.ind.pred.75 <- coef(wae.ind.75)[1]*x.1^coef(wae.ind.75)[2]      
wae.mod.pred.75 <- 1.266202e-06*x.1^3.349679

plot(weight~length, data=wae.ind, xlab="Length (mm)", xlim=c(100, 800), ylim=c(0,5000), 
      ylab="Weight (g)", pch="", xaxs="i", yaxs="i", bty="n", col="blue")
  abline(v=c(250,380,510,630,760), lty=2)
points(z~x, pch=19, col="blue") #independent data; 75th percentile as a function of 10mm length class
#points(weight~length, data=wae, pch=".", col="red")
  wae.moda <- seq(85,795)
    lines(wae.moda, predict(wae.mod.75, list(length=wae.moda)), lwd=2, col="red")
    lines(wae.moda, predict(wae.ind.75, list(length=wae.moda)), lwd=2, col="blue")
    lines(wae.moda, predict(wae.ind.5, list(length=wae.moda)), lwd=2, lty=2, col="blue")    

{
plot(weight~length, data=wae.ind, xlab="Length (mm)", xlim=c(100, 800), ylim=c(0,5000), 
      ylab="Weight (g)", pch="", xaxs="i", yaxs="i", bty="n", col="blue")
lines((10^emp.y)~x.1, lwd=2, col="red")
lines((10^rlp.y)~x.1, lwd=2, lty=2, col="red")
lines(wae.ind.pred.75~x.1, lwd=2)
lines(wae.mod.pred.75~x.1, lty=2, lwd=2, col="blue")
  legend(150,4750, bg="white", legend=c("Development 75th", "Independent 75th", "EmP", "RLP"), lty=c(2,1,1,2), lwd=2, col=c("blue", "black", "red", "red"))
}

with(wae.ind, 
  hist(length, breaks=100, xaxs="i"))
with(wae, 
  hist(length, breaks=100, add=T, col="red"))
  legend(550,3750, legend=c("Validation data", "Development data"), fill=c("white","red"), cex=1.1, bty="n")  

##################################################################################
mod <- rq(log10(weight)~log10(length), tau=1:99/100,  data=wae, model=T)
mod.1 <- rq(log10(weight)~log10(length), tau=c(0.25, 0.5, 0.75), data=wae, model=T)

system.time(plot(summary(mod, se="rank", iid=F, alpha=0.05), pch=19))
system.time(plot(summary(mod.1, se="rank", iid=F, alpha=0.05), pch=19, xlab=c("","Quantile"), 
  ylab=c("Estimate", "Estimate"), xlim=c(c(0,1),c(0,1)), main=c(expression(paste(log[10] ,beta[0])), expression(paste(log[10] ,beta[1]))), mar=c(5,4,4,2)))





################################################################################
#Now with Black crappie
rm(wae)

blc <- read.table('BLC Clean.txt', header=TRUE)
attach(blc)
names(blc)

psd <- with(blc,
            ifelse((length>=130)&(length<200), 130,
            ifelse((length>=200)&(length<250), 200,
            ifelse((length>=250)&(length<300), 250,
            ifelse((length>=300)&(length<380), 300,
            ifelse(length>=380, 380,
            0))))))

l.c <- with(blc,
  ifelse((length>=130)&(length<140),"135",
  ifelse((length>=140)&(length<150),"145",
  ifelse((length>=150)&(length<160),"155",
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
  NA))))))))))))))))))))))))))




blc.mod.75 <- nlrq(weight~alpha*length^beta, data=blc, tau=0.75,                         #Requires different parm. ests. for alpha and beta
      start=list(alpha=0.00001, beta=3), trace=TRUE)
blc.mod.5 <- nlrq(weight~alpha*length^beta, data=blc, tau=0.5,                           #Requires different parm. ests. for alpha and beta
      start=list(alpha=0.00001, beta=3), trace=TRUE)
blc.mod.25 <- nlrq(weight~alpha*length^beta, data=blc, tau=0.25,                         #Requires different parm. ests. for alpha and beta
      start=list(alpha=0.00001, beta=3), trace=TRUE)

#midpoint of each legnth category; double check that max length includes the entirety of the length category 
x <- as.vector(((c(199,249,299,380)-tapply(length, psd, min))/2)+tapply(length, psd, min)); x

x<- seq(135,375, by=10)

plot(weight~length, data=blc, xlab="Length (mm)", xlim=c(100, 400), ylim=c(0,1200), xaxs="i", yaxs="i", bty="n", 
      ylab="Weight (g)", pch="")
#  abline(v=c(130,200,250,300,380), lty=2)
  blc.moda <- seq(130,375)
    lines(blc.moda, predict(blc.mod.75, list(length=blc.moda)), lwd=2, lty=2)
    lines(blc.moda, predict(blc.mod.5, list(length=blc.moda)), lwd=2)
    lines(blc.moda, predict(blc.mod.25, list(length=blc.moda)), lwd=2, lty=2)
  points(tapply(weight, l.c, mean)~x, pch=1)
#add error bars to points
  errbar(x, 
         tapply(weight, l.c, mean), 
         tapply(weight, l.c, mean)+1.96*tapply(weight, l.c, sd), 
         tapply(weight, l.c, mean)-1.96*tapply(weight, l.c, sd), 
         pch=" ", add=T)
  text(175, 800, labels=expression(bold("Black\nCrappie")), cex=1.5)
#  text(x, c(300, 300, 175, 175), labels=c("S-Q", "Q-P", "P-M", "M-T")) 

#wae.pred.75 <- as.vector(predict(wae.mod.75, list(length=seq(205,745, by=10)), interval="confidence"))
blc.pred.75 <- as.vector(predict(blc.mod.75, list(length=seq(130,420, by=10)), interval="confidence"))
blc.pred.5 <- as.vector(predict(blc.mod.5, list(length=seq(130,420, by=10)), interval="confidence"))
blc.pred.25 <- as.vector(predict(blc.mod.25, list(length=seq(130,420, by=10)), interval="confidence"))

#creates a matrix 60 rows long by 3 columsn wide
blc.pred.values <- matrix(NA,30,3, dimnames=list(seq(130,420,by=10), c("Q1", "Q2", "Q3")))

#program to put all predicted values from each model into the matrix "pred.values"
for (i in 1:30){
  blc.pred.values[i,1]=as.vector(blc.pred.25[i])
  blc.pred.values[i,2]=as.vector(blc.pred.5[i])
  blc.pred.values[i,3]=as.vector(blc.pred.75[i])
  }

blc.pred.values <- as.data.frame(blc.pred.values, header=T)
blc.pred.values

#with(pred.values,
#  plot(Q1~seq(205,795,by=10)))

#create a table of intercept and slope values from each nlrq
blc.int.slope <- matrix(NA,3,2, dimnames=list(c("25th", "50th", "75th"), c("Intercept", "Slope")))
     blc.int.slope[1,1]=summary(blc.mod.25)$coef[1,1]
     blc.int.slope[1,2]=summary(blc.mod.25)$coef[2,1]
     blc.int.slope[2,1]=summary(blc.mod.5)$coef[1,1]
     blc.int.slope[2,2]=summary(blc.mod.5)$coef[2,1]
     blc.int.slope[3,1]=summary(blc.mod.75)$coef[1,1]
     blc.int.slope[3,2]=summary(blc.mod.75)$coef[2,1]

blc.int.slope

#####################################
#Plot 25th, 50th, and 75th P'tiles on the same page as EmP and RLP lines
x.1 <- seq(130, 370, 1)

#EmP and RLP equations are from Ranney et al. (2010)
emp.y <- ((-4.9346)+(2.734*log10(x.1))+((0.135)*(log10(x.1))^2))
rlp.y <- (-5.523)+((3.304)*log10(x.1))


plot(weight~length, data=blc, xlab="Length (mm)", xlim=c(100, 400), ylim=c(0,1200), xaxs="i", yaxs="i", bty="n", 
      ylab="Weight (g)", pch="")
  abline(v=c(130,200,250,300,380), lty=2)
  blc.moda <- seq(130,375)
    lines(blc.moda, predict(blc.mod.75, list(length=blc.moda)), lwd=2, lty=2)
    lines(blc.moda, predict(blc.mod.5, list(length=blc.moda)), lwd=2)
    lines(blc.moda, predict(blc.mod.25, list(length=blc.moda)), lwd=2, lty=2)
  legend(125, 1100, bg="white", legend=c("75th Percentile", "50th Percentile", "25th Percentile", "EmP", "RLP"), lty=c(2,1,2,1,1), lwd=2, col=c("black","black","black","red","blue"))
  text(175, 650, label=expression(bold("Black\nCrappie")), cex=2)

lines((10^emp.y)~x.1, col="red", lwd=2)#, cex=2)
lines((10^rlp.y)~x.1, col="blue", lwd=2)#, cex=2)
