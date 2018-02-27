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
# weight at length at the midpoints of the Gabelhouse length categories

wae_new <-
  data.frame(State = rep(c("ref", "GA2", "GA3", "GA4", "SD4", "SD13", "SD25"), 5), 
             length = rep(c(125, 315,445,570,695), each = 7), 
             weight = rep(NA, 35))

#Empty list to store values
predicted_output <- list()

taus <- seq(0.05, 0.95, by = 0.05)

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

# Convert list of output into a single dataframe
predicted_output <- 
  do.call("rbind", predicted_output) %>%
  as.data.frame() %>%
  rename(state = `wae_new$State`, 
         length = `wae_new$length`) %>%
  mutate(state = ifelse(state == "GA2", "GA1", 
                        ifelse(state == "GA3", "GA2", 
                               ifelse(state == "GA4", "GA3", 
                                      ifelse(state == "SD4", "SD1", 
                                             ifelse(state == "SD13", "SD2", 
                                                    ifelse(state == "SD25", "SD3", "Reference")))))), 
         state = state %>% as.factor(), 
         state = state %>% relevel(ref = "Reference"))

# Save to RDS because this takes a while, 8.3 minutes
predicted_output %>%
  saveRDS(paste0(Sys.Date(), "_predicted_weight_at_length.rds"))

predicted_output <- 
  readRDS(paste0(Sys.Date(), "_predicted_weight_at_length.rds")) %>%
  mutate(state = ifelse(state == "GA2", "GA1", 
                        ifelse(state == "GA3", "GA2", 
                               ifelse(state == "GA4", "GA3", 
                                      ifelse(state == "SD4", "SD1", 
                                             ifelse(state == "SD13", "SD2", 
                                                    ifelse(state == "SD25", "SD3", "Reference")))))), 
         state = state %>% as.factor(),
         state = state %>% relevel(ref = "Reference"))


# SD LAKES
predicted_output %>%
  rename(weight = fit) %>%
  mutate(length = paste0("TL = ", length)) %>%
  filter(state %in% c("Reference", "SD1", "SD2", "SD3")) %>%
  ggplot(aes(x = tau, y = weight, fill = state)) +
  geom_line(aes(linetype = state), lwd = 0.65) +
  geom_ribbon(aes(x = tau, ymin = lower, ymax = higher, fill = state, alpha = 0.05)) +
  facet_wrap(~length, scales = "free_y") +
  labs(x= "Quantile", y = "Weight (g)") +
  scale_fill_manual(name = "State", 
                    labels = c("Reference", "SD1", "SD2", "SD3"), 
                    values = gray.colors(4, start = 0.05, end = 0.8, gamma = 2.2, alpha = 0.5)) +
  scale_linetype_manual(name = "State", 
                        labels = c("Reference", "SD1", "SD2", "SD3"),
                        values = c(1,3,4,6)) +
  scale_alpha(guide = "none") +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(breaks = seq(0.05, 0.95, by = .15)) +
  theme_bw() +
  theme(legend.position = c(.9, .1), 
        legend.justification = c(1, 0), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(paste0("output/", Sys.Date(), "_sd_plots.png"))
ggsave(paste0("output/", Sys.Date(), "_sd_plots.tiff"))

#Color, .tiff
predicted_output %>%
  rename(weight = fit) %>%
  mutate(length = paste0("TL = ", length)) %>%
  filter(state %in% c("Reference", "SD1", "SD2", "SD3")) %>%
  ggplot(aes(x = tau, y = weight, fill = state)) +
  geom_line(aes(linetype = state), lwd = 0.65) +
  geom_ribbon(aes(x = tau, ymin = lower, ymax = higher, fill = state, alpha = 0.05)) +
  facet_wrap(~length, scales = "free_y") +
  labs(x= "Quantile", y = "Weight (g)") +
  scale_fill_manual(name = "State", 
                    labels = c("Reference", "SD1", "SD2", "SD3"), 
                    values = alpha(hue_pal()(4), alpha = 0.5)) +
  scale_linetype_manual(name = "State", 
                        labels = c("Reference", "SD1", "SD2", "SD3"), 
                        values = c(1,3,4,6)) +
  scale_alpha(guide = "none") +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(breaks = seq(0.05, 0.95, by = .15)) +
  theme_bw() +
  theme(legend.position = c(.9, .1), 
        legend.justification = c(1, 0), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(paste0("output/", Sys.Date(), "_sd_plots_color.png"))
ggsave(paste0("output/", Sys.Date(), "_sd_plots_color.tiff"))

  
# GA LAKES
predicted_output %>%
  rename(weight = fit) %>%
  mutate(length = paste0("TL = ", length)) %>%
  filter(state %in% c("Reference", "GA1", "GA2", "GA3")) %>%
  ggplot(aes(x = tau, y = weight, fill = state)) +
  geom_line(aes(linetype = state), lwd = 0.65) +
  geom_ribbon(aes(x = tau, ymin = lower, ymax = higher, fill = state, alpha = 0.05)) +
  facet_wrap(~length, scales = "free_y") +
  labs(x= "Quantile", y = "Weight (g)") +
  scale_fill_manual(name = "State", 
                    labels = c("Reference", "GA1", "GA2", "GA3"), 
                    values = gray.colors(4, start = 0.05, end = 0.8, gamma = 2.2, alpha = 0.5)) +
  scale_linetype_manual(name = "State", 
                        labels = c("Reference", "GA1", "GA2", "GA3"),
                        values = c(1,3,4,6)) +
  scale_alpha(guide = "none") +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(breaks = seq(0.05, 0.95, by = .15)) +
  theme_bw() +
  theme(legend.position = c(.9, .1), 
        legend.justification = c(1, 0), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(paste0("output/", Sys.Date(), "_ga_plots.png"))
ggsave(paste0("output/", Sys.Date(), "_ga_plots.tiff"))

#Color, .tiff
predicted_output %>%
  rename(weight = fit) %>%
  mutate(length = paste0("TL = ", length)) %>%
  filter(state %in% c("Reference", "GA1", "GA2", "GA3")) %>%
  ggplot(aes(x = tau, y = weight, fill = state)) +
  geom_line(aes(linetype = state)) +
  geom_ribbon(aes(x = tau, ymin = lower, ymax = higher, fill = state, alpha = 0.05)) +
  facet_wrap(~length, scales = "free_y") +
  labs(x= "Quantile", y = "Weight (g)") +
  scale_fill_manual(name = "State", 
                    labels = c("Reference", "GA1", "GA2", "GA3"), 
                    values = alpha(hue_pal()(4), alpha = 0.5)) +
  scale_linetype_manual(name = "State", 
                        labels = c("Reference", "GA1", "GA2", "GA3"), 
                        values = c(1,3,4,6)) +
  scale_alpha(guide = "none") +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(breaks = seq(0.05, 0.95, by = .15)) +
  theme_bw() +
  theme(legend.position = c(.9, .1), 
        legend.justification = c(1, 0), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  

ggsave(paste0("output/", Sys.Date(), "_ga_plots_color.png"))
ggsave(paste0("output/", Sys.Date(), "_ga_plots_color.tiff"))



