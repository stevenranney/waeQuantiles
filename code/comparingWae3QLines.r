################################################################################
# This code is annotated below.  It 
################################################################################ 
  
require(quantreg)
require(ggplot2)
library(purrr)
library(dplyr)
library(gplots)

source("code/helper_functions.R")

# For repeatability
set.seed(256)

# Read in and manipulate the WAE dataset from Ranney et al. 2010 and 2011; data
# filtering protocols have already been applied
wae <- read.table("data/wae_clean.txt", header = T)

# Assigns fish to a length category (Gabelhouse 1984)
wae <- 
  wae %>%
  mutate(psd = ifelse((length >= 250) & (length < 380), "S-Q",
                      ifelse((length >= 380) & (length < 510), "Q-P",
                             ifelse((length >= 510) & (length < 630), "P-M",
                                    ifelse((length >= 630) & (length < 760), "M-T",
                                           ifelse(length >= 760, ">T", "SS"))))))

# Assign fish to 10mm length class at midpoint
# round_down() fx is in helper_functions.R
wae <- 
  wae %>%
  mutate(l.c = round_down(length) +5)


# Build nonlinear quantile regression models;
# USED ONLY FOR CALCULATING A TABLE OF PREDICTED QUANTILE VALUES
wae.mod.75 <- nlrq(weight~alpha*length^beta, data=wae, tau=0.75,                         
                   start=list(alpha=0.00001, beta=3))
wae.mod.5 <- nlrq(weight~alpha*length^beta, data=wae, tau=0.5,                           
                  start=list(alpha=0.00001, beta=3))
wae.mod.25 <- nlrq(weight~alpha*length^beta, data=wae, tau=0.25,                         
                   start=list(alpha=0.00001, beta=3))


# Build linear quantile regression models
wae.linMod.75 <- rq(log10(weight)~log10(length), data=wae, tau=0.75)
wae.linMod.5 <- rq(log10(weight)~log10(length), data=wae, tau=0.5)
wae.linMod.25 <- rq(log10(weight)~log10(length), data=wae, tau=0.25)

# Create a data.frame of predicted values from the non-linear quantile regression models above
wae.pred.values <- 
  data.frame("25th" = as.vector(predict(wae.mod.25, list(length=seq(155,745, by=10)), interval="confidence")), 
             "50th" = as.vector(predict(wae.mod.5, list(length=seq(155,745, by=10)), interval="confidence")),
             "75th" = as.vector(predict(wae.mod.75, list(length=seq(155,745, by=10)), interval="confidence"))) %>%
  rename("25th" = `X25th`, 
         "50th" = `X50th`, 
         "75th" = `X75th`)

wae.pred.values %>%
  write.csv(paste0("output/", Sys.Date(), "_wae.pred.values.csv"))

#-------------------------------------------------------------------------------
# Read in independent data

waeInd <- read.table("data/wae_independent.txt", header=T)

# Assigns fish to a length category (Gabelhouse 1984)
waeInd <- 
  waeInd %>%
  mutate(psd = ifelse((length>=250)&(length<380), "S-Q",
                      ifelse((length>=380)&(length<510), "Q-P",
                             ifelse((length>=510)&(length<630), "P-M",
                                    ifelse((length>=630)&(length<760), "M-T",
                                           ifelse(length>=760, ">T", "SS"))))))

# Define GA and SD populations
waeGA <- 
  waeInd %>% 
  filter(State == "GA")

waeSD <- 
  waeInd %>%
  filter(State == "SD")

# Random GA populations
ga2 <- rq(log10(weight)~log10(length), data=waeGA %>% filter(lake == "2"), tau=0.75)
ga3 <- rq(log10(weight)~log10(length), data=waeGA %>% filter(lake == "3"), tau=0.75)
ga4 <- rq(log10(weight)~log10(length), data=waeGA %>% filter(lake == "4"), tau=0.75)


# Random SD populations
sd4 <- rq(log10(weight)~log10(length), data=waeSD %>% filter(lake == "4"), tau=0.75)
sd13 <- rq(log10(weight)~log10(length), data=waeSD %>% filter(lake == "13"), tau=0.75)
sd25 <- rq(log10(weight)~log10(length), data=waeSD %>% filter(lake == "25"), tau=0.75)


#-----------------------------------------------------------------------------

# Calculate the Blom estimator for Q3 values for all length categories
# BlomEstimator fx comes from helper_functions.R

# Reference population
sSRef <- BlomEstimator(wae$weight[wae$psd=="SS"], c(.25, .5,.75))
sQRef <- BlomEstimator(wae$weight[wae$psd=="S-Q"], c(.25, .5,.75))
qPRef <- BlomEstimator(wae$weight[wae$psd=="Q-P"], c(.25, .5,.75))
pMRef <- BlomEstimator(wae$weight[wae$psd=="P-M"], c(.25, .5,.75))
mTRef <- BlomEstimator(wae$weight[wae$psd=="M-T"], c(.25, .5,.75))
gTRef <- c(NA, NA, NA)

# GA 1
sSGa2 <- BlomEstimator(waeGA$weight[waeGA$psd=="SS" & waeGA$lake == "2"], c(.25, .5, .75))
sQGa2 <- BlomEstimator(waeGA$weight[waeGA$psd=="S-Q" & waeGA$lake == "2"], c(.25, .5, .75))
qPGa2 <- BlomEstimator(waeGA$weight[waeGA$psd=="Q-P" & waeGA$lake == "2"], c(.25, .5, .75))
pMGa2 <- BlomEstimator(waeGA$weight[waeGA$psd=="P-M" & waeGA$lake == "2"], c(.25, .5, .75))
mTGa2 <- BlomEstimator(waeGA$weight[waeGA$psd=="M-T" & waeGA$lake == "2"], c(.25, .5, .75))
gTGa2 <- c(NA, NA, NA)

# GA 2
sSGa3 <- BlomEstimator(waeGA$weight[waeGA$psd=="SS" & waeGA$lake == "3"], c(.25, .5, .75))
sQGa3 <- BlomEstimator(waeGA$weight[waeGA$psd=="S-Q" & waeGA$lake == "3"], c(.25, .5, .75))
qPGa3 <- BlomEstimator(waeGA$weight[waeGA$psd=="Q-P" & waeGA$lake == "3"], c(.25, .5, .75))
pMGa3 <- BlomEstimator(waeGA$weight[waeGA$psd=="P-M" & waeGA$lake == "3"], c(.25, .5, .75))
mTGa3 <- BlomEstimator(waeGA$weight[waeGA$psd=="M-T" & waeGA$lake == "3"], c(.25, .5, .75))
gTGa3 <- c(NA, NA, NA)

# GA 3
sSGa4 <- BlomEstimator(waeGA$weight[waeGA$psd=="SS" & waeGA$lake == "4"], c(.25, .5, .75))
sQGa4 <- BlomEstimator(waeGA$weight[waeGA$psd=="S-Q" & waeGA$lake == "4"], c(.25, .5, .75))
qPGa4 <- BlomEstimator(waeGA$weight[waeGA$psd=="Q-P" & waeGA$lake == "4"], c(.25, .5, .75))
pMGa4 <- BlomEstimator(waeGA$weight[waeGA$psd=="P-M" & waeGA$lake == "4"], c(.25, .5, .75))
mTGa4 <- BlomEstimator(waeGA$weight[waeGA$psd=="M-T" & waeGA$lake == "4"], c(.25, .5, .75))
gTGa4 <- c(NA, NA, NA)

# SD 1
sSSd4 <- BlomEstimator(waeSD$weight[waeSD$psd=="SS" & waeSD$lake == "4"], c(.25, .5, .75))
sQSd4 <- BlomEstimator(waeSD$weight[waeSD$psd=="S-Q" & waeSD$lake == "4"], c(.25, .5, .75))
qPSd4 <- BlomEstimator(waeSD$weight[waeSD$psd=="Q-P" & waeSD$lake == "4"], c(.25, .5, .75))
pMSd4 <- BlomEstimator(waeSD$weight[waeSD$psd=="P-M" & waeSD$lake == "4"], c(.25, .5, .75))
mTSd4 <- BlomEstimator(waeSD$weight[waeSD$psd=="M-T" & waeSD$lake == "4"], c(.25, .5, .75))
gTSd4 <- c(NA, NA, NA)

# SD 2
sSSd13 <- BlomEstimator(waeSD$weight[waeSD$psd=="SS" & waeSD$lake == "13"], c(.25, .5, .75))
sQSd13 <- BlomEstimator(waeSD$weight[waeSD$psd=="S-Q" & waeSD$lake == "13"], c(.25, .5, .75))
qPSd13 <- BlomEstimator(waeSD$weight[waeSD$psd=="Q-P" & waeSD$lake == "13"], c(.25, .5, .75))
pMSd13 <- BlomEstimator(waeSD$weight[waeSD$psd=="P-M" & waeSD$lake == "13"], c(.25, .5, .75))
mTSd13 <- BlomEstimator(waeSD$weight[waeSD$psd=="M-T" & waeSD$lake == "13"], c(.25, .5, .75))
gTSd13 <- c(NA, NA, NA)

# SD 3
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

# Write the quarTable to the output directory
quarTable %>%
  write.csv(paste0("output/", Sys.Date(), "_quarTable.csv"), row.names=T)

#-----------------------------------------------------------------------------
# Estimate differences between Reference slope and intercept and other models

# Bootstrap the sampling distributions of slope and intercept values by 
# running n regressions with numSamples drawn randomly with replacement from 
# a given dataset
# EstLinQuantRegCoefDist fx comes from the helper_functions.R file

refDist <- EstLinQuantRegCoefDist(wae, 10000, 100, tau=0.75)

ga2Dist <- EstLinQuantRegCoefDist(waeGA %>% filter(lake == "2"), 10000, 100, 0.75)

ga3Dist <- EstLinQuantRegCoefDist(waeGA %>% filter(lake == "3"), 10000, 100, 0.75)

ga4Dist <- EstLinQuantRegCoefDist(waeGA %>% filter(lake == "4"), 10000, 100, 0.75)

sd4Dist <- EstLinQuantRegCoefDist(waeSD %>% filter(lake == "4"), 10000, 100, 0.75)

sd13Dist <- EstLinQuantRegCoefDist(waeSD %>% filter(lake == "13"), 10000, 100, 0.75)

sd25Dist <- EstLinQuantRegCoefDist(waeSD %>% filter(lake == "25"), 10000, 100, 0.75)

# Create a table of slope coefficients and bootstrapped confidence intervals
slopeVals <- as.data.frame(matrix(data = NA, ncol=3, nrow=7, byrow=T, 
                           dimnames = list(c("Reference", "GA 1", 
                                             "GA 2", "GA 3", "SD 1",
                                             "SD 2", "SD 3"), 
                                           c("2.5", "mean", "97.5"))))

slopeVals <- 
  data.frame(`2.5` = c(BlomEstimator(refDist[,2], 0.025), 
                       BlomEstimator(ga2Dist[,2], 0.025),
                       BlomEstimator(ga3Dist[,2], 0.025),
                       BlomEstimator(ga4Dist[,2], 0.025),
                       BlomEstimator(sd4Dist[,2], 0.025),
                       BlomEstimator(sd13Dist[,2], 0.025),
                       BlomEstimator(sd25Dist[,2], 0.025)), 
             mean = c(mean(refDist[,2]), mean(ga2Dist[,2]), mean(ga3Dist[,2]), 
                      mean(ga4Dist[,2]), mean(sd4Dist[,2]), mean(sd13Dist[,2]), 
                      mean(sd25Dist[,2])), 
             `97.5` = c(BlomEstimator(refDist[,2], 0.975), 
                        BlomEstimator(ga2Dist[,2], 0.975),
                        BlomEstimator(ga3Dist[,2], 0.975),
                        BlomEstimator(ga4Dist[,2], 0.975),
                        BlomEstimator(sd4Dist[,2], 0.975),
                        BlomEstimator(sd13Dist[,2], 0.975),
                        BlomEstimator(sd25Dist[,2], 0.975))) %>%
  mutate(pop = c("Reference", "GA 1", "GA 2", "GA 3", "SD 1", "SD 2", "SD 3")) %>%
  rename("2.5" = "X2.5", 
         "97.5" = "X97.5")

slopeVals %>%
  ggplot(aes(x = pop, y = mean)) +
  geom_point() +
  labs(x = "Population", y = "Slope") +
  geom_errorbar(aes(x = pop, ymin = `2.5`, ymax = `97.5`)) +
  geom_hline(yintercept = slopeVals %>% filter(pop == "Reference") %>% .$`2.5`, 
             linetype = 2, size = 0.1) +
  geom_hline(yintercept = slopeVals %>% filter(pop == "Reference") %>% .$`97.5`, 
             linetype = 2, size = 0.1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Create a table of intercept coefficients and bootstrapped confidence intervals

intVals <- 
  data.frame(`2.5` = c(BlomEstimator(refDist[,1], 0.025), 
                       BlomEstimator(ga2Dist[,1], 0.025),
                       BlomEstimator(ga3Dist[,1], 0.025),
                       BlomEstimator(ga4Dist[,1], 0.025),
                       BlomEstimator(sd4Dist[,1], 0.025),
                       BlomEstimator(sd13Dist[,1], 0.025),
                       BlomEstimator(sd25Dist[,1], 0.025)), 
             mean = c(mean(refDist[,1]), mean(ga2Dist[,1]), mean(ga3Dist[,1]), 
                      mean(ga4Dist[,1]), mean(sd4Dist[,1]), mean(sd13Dist[,1]), 
                      mean(sd25Dist[,1])), 
             `97.5` = c(BlomEstimator(refDist[,1], 0.975), 
                        BlomEstimator(ga2Dist[,1], 0.975),
                        BlomEstimator(ga3Dist[,1], 0.975),
                        BlomEstimator(ga4Dist[,1], 0.975),
                        BlomEstimator(sd4Dist[,1], 0.975),
                        BlomEstimator(sd13Dist[,1], 0.975),
                        BlomEstimator(sd25Dist[,1], 0.975))) %>%
  mutate(pop = c("Reference", "GA 1", "GA 2", "GA 3", "SD 1", "SD 2", "SD 3")) %>%
  rename("2.5" = "X2.5", 
         "97.5" = "X97.5")

intVals %>%
  ggplot(aes(x = pop, y = mean)) +
  geom_point() +
  labs(x = "Population", y = "Intercept") +
  geom_errorbar(aes(x = pop, ymin = `2.5`, ymax = `97.5`)) +
  geom_hline(yintercept = intVals %>% filter(pop == "Reference") %>% .$`2.5`, 
             linetype = 2, size = 0.1) +
  geom_hline(yintercept = intVals %>% filter(pop == "Reference") %>% .$`97.5`, 
             linetype = 2, size = 0.1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


# Combine the slope and intercept values into one table
regVals <- 
  bind_cols(intVals %>%
          rename(`int2.5` = `2.5`, 
                 intmean = mean, 
                 `int97.5` = `97.5`) %>%
            select(int2.5, intmean, int97.5), 
        slopeVals %>%
          rename(`s2.5` = `2.5`, 
                 smean = mean, 
                 `s97.5` = `97.5`) %>%
          select(s2.5, smean, s97.5, pop))

#Write tables to file
regVals %>%
  write.csv(paste0("output/", Sys.Date(), "_regVals.csv"))

# Do the 95% estimates of slope for each population overlap the slope value for the reference population?

regVals <- 
  regVals %>%
  mutate(row_num = row_number(),
         is_slope_overlap = ifelse(smean > (regVals %>% filter(pop == "Reference") %>% .$s2.5) &
                               smean < (regVals %>% filter(pop == "Reference") %>% .$s97.5), TRUE, FALSE)) %>%
  mutate(pop = pop %>% factor(levels = c("Reference", "GA 1", "GA 2", "GA 3", 
                                            "SD 1", "SD 2", "SD 3")))


#-----------------------------------------------------------  
# CONVERT BELOW TO A GGPLOT WITH H AND V ERRBARS

xlab <- bquote(.("Slope" ~ (beta[1])))
ylab <- bquote(.("Intercept" ~ (beta[0])))

regVals %>%
  ggplot(aes(x = smean, y = intmean, shape = pop)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = int2.5, ymax = int97.5), width = 0.01) +
  geom_errorbarh(aes(xmin = s2.5, xmax = s97.5), height = 0.07) + 
  scale_x_continuous(breaks = seq(3.0, 3.5, 0.05)) +
  scale_y_continuous(breaks = seq(-6.4, -5, 0.2)) +
  scale_size(guide = "none") +
  xlab(bquote(.("Slope" ~ (beta[1])))) +
  ylab(bquote(.("Intercept" ~ (beta[0])))) +
  scale_shape_manual(name = "Population", 
                     values = c("Reference" = 13, 
                                "GA 1" = 0, 
                                "GA 2" = 1,
                                "GA 3" = 2, 
                                "SD 1" = 15,
                                "SD 2" = 16, 
                                "SD 3" = 17)) +
  theme(legend.position = "bottom", 
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))


#-----------------------------------------------------------  


       
# IF REFERENCE POPULATION SLOPE AND INTERCEPT 95% CI DO NOT OVERLAP WITH POPULATION
# LEVEL SLOPE AND INTERCEPT, DETERMINE WHETHER THE INTERVAL AROUND THE DIFFERENCE
# IN SLOPES CONTAINS ZERO; FROM POPE AND KRUSE 2007 AIFFD, page 433 (i.e., rows 
# that are FALSE in regVals$is_slope_overlap)


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
  mutate(pop = c("GA 1", "GA 2", "GA 3", "SD 1", "SD 2", "SD 3"))


