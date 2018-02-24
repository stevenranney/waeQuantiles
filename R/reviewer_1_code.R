###S. Ranney walleye weight length data Feb 2018
###Code by B. S. Cade, 8 Feb 2018

library(dplyr)
library(ggplot2)
library(scales)

walleye.ref <-
  read.table("data/wae_clean.txt",
             header=T)

walleye.state <-
  read.table("data/wae_independent.txt",
             header=T)

walleye.state6 <- 
  walleye.state %>% 
  filter((walleye.state$State=="SD" & walleye.state$lake==4)|
          (walleye.state$State=="SD" & walleye.state$lake==13)|
          (walleye.state$State=="SD" & walleye.state$lake==25)|
          (walleye.state$State=="GA" & walleye.state$lake==2)|
          (walleye.state$State=="GA" & walleye.state$lake==3)|
          (walleye.state$State=="GA" & walleye.state$lake==4))

walleye.comb <- 
  walleye.ref %>% 
  select(length, weight, lake)

walleye.comb <- 
  cbind(State="ref",walleye.comb)

walleye.comb <- 
  walleye.comb %>% 
  bind_rows(walleye.state6)

walleye.comb <-
  walleye.comb %>%
  mutate(State = State %>% as.character())

walleye.comb <-
  walleye.comb %>% 
  mutate(State = ifelse(State == "ref", "ref", 
                        ifelse(State == "GA" & lake == 2, "GA2", 
                               ifelse(State == "GA" & lake == 3, "GA3", 
                                      ifelse(State == "GA" & lake == 4, "GA4", 
                                             ifelse(State == "SD" & lake == 4, "SD4", 
                                                    ifelse(State == "SD" & lake == 13, "SD13", "SD25")))))), 
         State = State %>% as.factor())

walleye.comb %>%
  ggplot(aes(x = log10(length), y = log10(weight))) +
  geom_point(alpha = 0.25) + 
  facet_wrap(~State)

################################################################################
# Creating quantile regressions

library(quantreg)

#Quantile regressions from the reference data
w.ref.75 <- rq(log10(weight)~log10(length), data = walleye.ref, tau = 0.75)
w.ref.50 <- rq(log10(weight)~log10(length), data = walleye.ref, tau = 0.50)
w.ref.25 <- rq(log10(weight)~log10(length), data = walleye.ref, tau = 0.25)

###five quantiles at once with their predictions by 10mm increments
w.ref.10to90 <-
  rq(log10(weight)~log10(length), data = walleye.ref, tau = c(0.10, 0.25, 0.50, 0.75, 0.90))

by10mm <- seq(155, 745, by = 10)
by10mm <- data.frame(length = by10mm)

predict.by.10mm <- predict(w.ref.10to90, newdata = by10mm, confidence = none)
predict.by.10mm <- 10^(predict.by.10mm) %>% round(1) %>% comma()
predict.by.10mm <- 
  cbind(by10mm, predict.by.10mm) %>%
  rename(`Total length (mm)` = length, 
         `0.10` = "tau= 0.10", 
         `0.25` = "tau= 0.25", 
         `0.50` = "tau= 0.50", 
         `0.75` = "tau= 0.75", 
         `0.90` = "tau= 0.90") 

predict.by.10mm %>%
  write.csv(paste0("output/", Sys.Date(), "_predicted_values.csv"), 
            row.names = FALSE)
         

###Three different approaches for obtaining standard errors of esitmates
###se="nid" is default approach for large sample sizes and based on the 
###sandwich form of the variance/covariance matrix for heterogeneous 
###models.
###se="wxy" is a weighted exponential tilting bootstrap using entire 
###sample.
###se="xy",R=1000, mofn=5000 is bootstrap of xy-pairs 5000 of n samples 
###made 1000 times.
summary(w.ref.75,se="nid")
summary(w.ref.75,se="boot",bsmethod="wxy",R=100)
summary(w.ref.75,se="boot",bsmethod="xy",R=1000,mofn=5000)

###Combined 6 states and reference data
###Make the ref data the base level in this estimate of 0.75 quantile.
walleye.comb <-
  walleye.comb %>%
  mutate(State = State %>% relevel(ref = "ref"))

w.comb.75 <- 
  walleye.comb %>% 
  rq(log10(weight)~log10(length) + State + log10(length):State, data = ., 
     contrasts = list(State="contr.treatment"), tau = 0.75)

w.comb.75.estimates <- summary(w.comb.75, se = "boot", bsmethod = "xy", R = 1000, mofn = 5000)

w.comb.75.estimates <- w.comb.75.estimates$coef
w.comb.75.estimates <- data.frame(w.comb.75.estimates)

###Calculate 95% confidence intervals using bootstrap estimates of SE.
resid.df <- nrow(w.comb.75$x) - ncol(w.comb.75$x)

w.comb.75.estimates <- 
  w.comb.75.estimates %>%
  mutate(Lwr95CI = Value + Std..Error * qt(0.025,resid.df), 
         Upr95CI = Value + Std..Error * qt(0.975,resid.df), 
         name = row.names(.)) %>%
  select(name, Lwr95CI, Value, Upr95CI, t.value, Pr...t..)

###Now estimating alternative form of same model so that estimates are 
###obtained for each population by removing the intercept term from
###model.
w.comb.75.alt <- 
  walleye.comb %>% 
  rq(log10(weight)~State + log10(length):State - 1, data = ., 
     contrasts = list(State = "contr.treatment"), tau = 0.75)

w.comb.75.alt.estimates <- summary(w.comb.75.alt, se = "boot", bsmethod = "xy", R = 1000, mofn = 5000)
w.comb.75.alt.estimates <- w.comb.75.alt.estimates$coef
w.comb.75.alt.estimates <- data.frame(w.comb.75.alt.estimates)

###Calculate 95% confidence intervals using bootstrap estimates of SE.
resid.df <- nrow(w.comb.75.alt$x) - ncol(w.comb.75.alt$x)
w.comb.75.alt.estimates <- 
  w.comb.75.alt.estimates %>%
  mutate(Lwr95CI = Value + Std..Error * qt(0.025, resid.df), 
         Upr95CI = Value + Std..Error * qt(0.975,resid.df), 
         names = row.names(.)) %>%
  select(names, Lwr95CI, Value, Upr95CI)

###To estimate predictions at length
###Create new data set for predictions using midpoint lengths of length
###intervals in Ranney 2018
walleye.comb.new <-
  data.frame(State = rep(c("ref", "GA2", "GA3", "GA4", "SD4", "SD13", "SD25"), 5), 
             length = rep(c(125, 315,445,570,695), each = 7), 
             weight = rep(NA, 35))
# walleye.comb.new <- walleye.comb[1:35,]
# walleye.comb.new$weight <- "na"
# walleye.comb.new$length[1:7] <- 125
# walleye.comb.new$length[8:14] <- 315
# walleye.comb.new$length[15:21] <- 445
# walleye.comb.new$length[22:28] <- 570
# walleye.comb.new$length[29:35] <- 695
# walleye.comb.new$State[c(1,8,15,22,29)] <- "ref"
# walleye.comb.new$State[c(2,9,16,23,30)] <- "GA2"
# walleye.comb.new$State[c(3,10,17,24,31)] <- "GA3"
# walleye.comb.new$State[c(4,11,18,25,32)] <- "GA4"
# walleye.comb.new$State[c(5,12,19,26,33)] <- "SD4"
# walleye.comb.new$State[c(6,13,20,27,34)] <- "SD13"
# walleye.comb.new$State[c(7,14,21,28,35)] <- "SD25"
# walleye.comb.new <- walleye.comb.new[,-5]
walleye.comb.pred.75 <- predict(w.comb.75,newdata=walleye.comb.new,
                                type="percentile",se="boot",bsmethod="xy", R
                                =1000,mofn=5000,interval="confidence",level=0.95)

###Exponentiate to get weights into g.
walleye.comb.pred.cat.midpoints.75 <- 10^walleye.comb.pred.75
walleye.comb.pred.cat.midpoints.75 <-
  data.frame(walleye.comb.pred.cat.midpoints.75)
walleye.comb.pred.cat.midpoints.75 <-
  cbind(walleye.comb.new$State,walleye.comb.new$length,walleye.comb.pred.cat.midpoints.75)

###To create categorical factor for length intervals similar to what 
###Ranney used in his Table 4.
walleye.comb$lengthCat <- "SS"
walleye.comb$lengthCat[walleye.comb$length >=250 & walleye.comb$length
                       <380] <- "S-Q"
walleye.comb$lengthCat[walleye.comb$length >=380 & walleye.comb$length
                       <510] <- "Q-P"
walleye.comb$lengthCat[walleye.comb$length >=510 & walleye.comb$length
                       <630] <- "P-M"
walleye.comb$lengthCat[walleye.comb$length >=630 & walleye.comb$length
                       <760] <- "M-T"
walleye.comb$lengthCat[walleye.comb$length >=760] <- "T"
walleye.comb$lengthCat <- as.factor(walleye.comb$lengthCat)
walleye.comb.cat.75 <- rq(log10(weight)~lengthCat:State -1,
                          data=walleye.comb,tau=0.75,contrasts=list(lengthCat="contr.treatment",
                                                                    State="contr.treatment"))
10^(walleye.comb.cat.75$coef)
walleye.comb.cat.50 <- rq(log10(weight)~lengthCat:State -1,
                          data=walleye.comb,tau=0.50,contrasts=list(lengthCat="contr.treatment",
                                                                    State="contr.treatment"))
10^(walleye.comb.cat.50$coef)
walleye.comb.cat.25 <- rq(log10(weight)~lengthCat:State -1,
                          data=walleye.comb,tau=0.25,contrasts=list(lengthCat="contr.treatment",
                                                                    State="contr.treatment"))
10^(walleye.comb.cat.25$coef)

###Nonlinear QR estimates. Do not use these as they don't correspond to 
###log10 transformed estimates, are unstable, lack flexibility for 
###including additional parameters, and have poorly developed infernce 
###procedures compared to the linear QR model.
w.ref.nlrq75 <-
  nlrq(weight~b0*length^b1,data=walleye.ref,tau=0.75,start=list(b0=0.0001,b1=3)
  )
w.ref.nlrq50 <-
  nlrq(weight~b0*length^b1,data=walleye.ref,tau=0.50,start=list(b0=0.0001,b1=3)
  )
w.ref.nlrq25 <-
  nlrq(weight~b0*length^b1,data=walleye.ref,tau=0.25,start=list(b0=0.0001,b1=3)
  )
predict(w.ref.nlrq75,newdata=walleye.comb.new)