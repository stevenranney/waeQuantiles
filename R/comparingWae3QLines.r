#-----------------------------------------------------------------------------
# Code to accompany the manuscript 
# Quantile Regression Estimates of Body Weight for Walleye
# S. H. Ranney, 2018
#-----------------------------------------------------------------------------

library(quantreg)
library(ggplot2)
library(purrr)
library(dplyr)
library(gplots)
library(scales)

source("R/helper_functions.R")

# For repeatability
set.seed(256)

# Data handling. Read in the reference and state "independent" datasets, manipulate, 
# and combine
wae <- read.table("data/wae_clean.txt", header = T)

wae_ref <- 
  wae %>% 
  mutate(State = "ref", 
         State = State %>% as.factor) %>%
  select(State, length, weight, lake)

waeInd <- 
  read.table("data/wae_independent.txt", header=T) %>%
  filter((State=="SD" & lake==4)|
           (State=="SD" & lake==13)|
           (State=="SD" & lake==25)|
           (State=="GA" & lake==2)|
           (State=="GA" & lake==3)|
           (State=="GA" & lake==4))

wae <- 
  wae_ref %>%
  bind_rows(waeInd) %>%
  mutate(State = ifelse(State == "ref", "ref", 
                        ifelse(State == "GA" & lake == 2, "GA2", 
                               ifelse(State == "GA" & lake == 3, "GA3", 
                                      ifelse(State == "GA" & lake == 4, "GA4", 
                                             ifelse(State == "SD" & lake == 4, "SD4", 
                                                    ifelse(State == "SD" & lake == 13, "SD13", "SD25")))))), 
         State = State %>% as.factor(), 
         psd = assign_wae_psd(length), 
         l.c = round_down(length) +5)

wae %>%
  ggplot(aes(x = log10(length), y = log10(weight))) +
  geom_point(alpha = 0.25) + 
  facet_wrap(~State)

#-----------------------------------------------------------------------------
# Table of quantiles and predictions of weight-at-length

# Five quantiles at once with their predictions by 10mm increments
wae_ref_10to90 <-
  wae_ref %>%
  rq(log10(weight)~log10(length), data = ., tau = c(0.10, 0.25, 0.50, 0.75, 0.90))

by10mm <- 
  data.frame(length = seq(155, 745, by = 10))

predict_by_10mm <- 
  predict(wae_ref_10to90, newdata = by10mm, confidence = none)

# Create, as much as possible, the prediction table in R so not as much 
# Excel or Word formatting needs to be done.
predict_by_10mm <- 
  10^(predict_by_10mm) %>% #Exponentiate all values
  round(1) %>% #Round to 1 decimal place
  comma() %>% #Add comma
  cbind(by10mm, .) %>%
  rename(`Total length (mm)` = length, 
         `0.10` = "tau= 0.10", 
         `0.25` = "tau= 0.25", 
         `0.50` = "tau= 0.50", 
         `0.75` = "tau= 0.75", 
         `0.90` = "tau= 0.90")

predict_by_10mm %>%
  write.csv(paste0("output/", Sys.Date(), "_predicted_values.csv"), 
            row.names = FALSE)


#-----------------------------------------------------------------------------
# Approach for obtaining estimates of the differences among populations
# se="xy",R=1000, mofn=5000 is bootstrap of xy-pairs 5000 of n samples 
# made 1000 times.

# Make the ref data the base level in this estimate of 0.75 quantile.
wae <-
  wae %>%
  mutate(State = State %>% relevel(ref = "ref"))

wae_75 <- 
  wae %>% 
  rq(log10(weight)~log10(length) + State + log10(length):State, data = ., 
     contrasts = list(State="contr.treatment"), tau = 0.75)

wae_75_diff <- summary(wae_75, se = "boot", bsmethod = "xy", R = 1000, mofn = 5000)

wae_75_diff <- 
  data.frame(wae_75_diff$coef) %>%
  mutate(name = row.names(.)) %>%
  select(name, Value, Std..Error, t.value, Pr...t..) %>%
  rename(Estimate = Value, 
         SE = `Std..Error`,
         `t value` = t.value, 
         `p value` = `Pr...t..`)

# Calculate 95% confidence intervals around the estimate of the differences in 
# slope/int among populationsusing bootstrap estimates of SE.
resid_df <- nrow(wae_75$x) - ncol(wae_75$x)

wae_75_diff <- 
  wae_75_diff %>%
  mutate(Lwr95CI = Estimate + SE * qt(0.025,resid_df), 
         Upr95CI = Estimate + SE * qt(0.975,resid_df)) %>%
  select(name, Lwr95CI, Estimate, Upr95CI, `t value`, `p value`)

wae_75_diff %>%
  write.csv(paste0("output/", Sys.Date(), "_differences_in_slope_int.csv"), 
            row.names = FALSE)

# Retrieve slope and intercept for each population
# Same model as above but removing the intercept term so that I can find slope/int
# estimates for each population, including ref

wae_75_slope_int <- 
  wae %>% 
  rq(log10(weight)~State + log10(length):State - 1, data = ., 
     contrasts = list(State = "contr.treatment"), tau = 0.75)

wae_75_slope_int_est <- summary(wae_75_slope_int, se = "boot", bsmethod = "xy", R = 1000, mofn = 5000)

wae_75_slope_int_est <- 
  data.frame(wae_75_slope_int_est$coefficients) %>%
  mutate(name = row.names(.)) %>%
  select(name, Value, Std..Error, t.value, Pr...t..) %>%
  rename(`Point estimate` = Value, 
         SE = `Std..Error`,
         `t value` = t.value, 
         `p value` = `Pr...t..`)


###Calculate 95% confidence intervals using bootstrap estimates of SE.
resid_df <- nrow(wae_75_slope_int$x) - ncol(wae_75_slope_int$x)

wae_75_slope_int_est <- 
  wae_75_slope_int_est %>%
  mutate(Lwr95CI = `Point estimate` + SE * qt(0.025, resid_df), 
         Upr95CI = `Point estimate` + SE * qt(0.975, resid_df)) %>%
  select(name, Lwr95CI, `Point estimate`, Upr95CI)

wae_75_slope_int_est %>%
  write.csv(paste0("output/", Sys.Date(), "_slope_int_estimates.csv"), 
            row.names = FALSE)


#-----------------------------------------------------------------------------
# For every population from every quantile seq(0.05, 0.95, by = 0.05), predict 
# weight at lenght at the midpoints of the Gabelhouse length categories

wae_new <-
  data.frame(State = rep(c("ref", "GA2", "GA3", "GA4", "SD4", "SD13", "SD25"), 5), 
             length = rep(c(125, 315,445,570,695), each = 7), 
             weight = rep(NA, 35))

#Empty list to store values
predicted_output <- list()

taus <- seq(0.05, 0.95, by = 0.05)

start <- Sys.time()

for(i in 1:length(taus)){
  
  #Create model for each tau
  wae_mod <- 
    wae %>% 
    rq(log10(weight)~log10(length) + State + log10(length):State, data = ., 
       contrasts = list(State="contr.treatment"), tau = taus[i])
  
  #predict weights at length in wae_new for each tau 
  wae_pred <- predict(wae_mod, newdata = wae_new,
                      type = "percentile", se = "boot", bsmethod = "xy", R = 1000,
                      mofn = 5000, interval = "confidence", level = 0.95)
  
  #exponentiate weights into g
  wae_pred_midpoints <- 10^wae_pred
  wae_pred_midpoints <- data.frame(wae_pred_midpoints)
  wae_pred_midpoints <-
    cbind(wae_new$State,wae_new$length,wae_pred_midpoints) %>%
    mutate(tau = taus[i])
  
  #Store output in the empty list
  predicted_output[[i]] <- wae_pred_midpoints
  
  }

end <- Sys.time()

end - start

# Convert list of output into a single dataframe
predicted_output <- 
  do.call("rbind", predicted_output ) %>%
  as.data.frame() %>%
  rename(state = `wae_new$State`, 
         length = `wae_new$length`)

# Save to RDS because this takes a while, 8.3 minutes
predicted_output %>%
  saveRDS(paste0(Sys.Date(), "_predicted_weight_at_length.rds"))


# Figure of predicted weights at specific lengths at all quantiles
predicted_output %>%
  rename(weight = fit) %>%
  mutate(length = paste0("TL = ", length)) %>%
  filter(state != "ref") %>%
  ggplot(aes(x = tau, y = weight, colour = state)) +
#  geom_point() +
  geom_line() +
  geom_ribbon(aes(x = tau, ymin = lower, ymax = higher, fill = state, alpha = 0.05)) +
  facet_wrap(~length, scales = "free_y")
  









#-----------------------------------------------------------------------------
# BELOW IS LIKELY SUPERFLUOUS NOW

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

ggsave(paste0("output/", Sys.Date(), "_slopes.png"))

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

ggsave(paste0("output/", Sys.Date(), "_intercepts.png"))


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
  write.csv(paste0("output/", Sys.Date(), "_regVals.csv"), row.names = F)

# Do the 95% estimates of slope for each population overlap the slope value for the reference population?

regVals <- 
  regVals %>%
  mutate(row_num = row_number(),
         is_slope_overlap = ifelse(smean > (regVals %>% filter(pop == "Reference") %>% .$s2.5) &
                               smean < (regVals %>% filter(pop == "Reference") %>% .$s97.5), TRUE, FALSE), 
         is_int_overlap = ifelse(intmean > (regVals %>% filter(pop == "Reference") %>% .$int2.5) &
                                     intmean < (regVals %>% filter(pop == "Reference") %>% .$int97.5), TRUE, FALSE)) %>%
  mutate(pop = pop %>% factor(levels = c("Reference", "GA 1", "GA 2", "GA 3", 
                                            "SD 1", "SD 2", "SD 3")))

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

slope_overlap <- 
  list(ga2Dist, ga3Dist, ga4Dist, 
     sd4Dist, sd13Dist, sd25Dist) %>%
  map_dfr(interval_around_differences, refDist) %>%
  mutate(pop = c("GA 1", "GA 2", "GA 3", "SD 1", "SD 2", "SD 3")) %>%
  full_join(regVals)

slope_overlap %>%
  write.csv(paste0("output/", Sys.Date(), "_slope_overlap.csv"), row.names = F)
  


#-----------------------------------------------------------  

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
  guides(shape = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom", 
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

ggsave(paste0("output/", Sys.Date(), "_main_fig.png"))
ggsave(paste0("output/", Sys.Date(), "_main_fig.tiff"))



#-----------------------------------------------------------  


