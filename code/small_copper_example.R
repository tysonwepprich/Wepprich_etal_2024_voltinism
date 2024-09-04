# Note: this code has been modified from   
# Fagard-Jenkin C (2022). _rGAI: Generalised Abundance Index for Seasonal Invertibrates_. R package version 1.0.
# https://github.com/calliste-fagard-jenkin/rGAI/tree/master/examples
# Modification is in lines 126-324 to demonstrate model selection 

# Author E. B. Dennis

# Example 2 - modelling multiple broods
# This example demonstrates use of the rGAI package to a species with
# three broods per year by applications to Small Copper Lycaena phlaeas

# Note: Modified to try out 4 and 5 broods to see if methods choose correct number of broods with AIC

# Load packages
library(rGAI)
library(ggplot2)
library(dplyr)

# Read in example data for Small Copper
data("smallcopper_data")
head(smallcopper_data)

# Summarise and plot the observed count data
count_summary <- smallcopper_data %>% group_by(occasion) %>%
                      summarise(mean = mean(count, na.rm=TRUE),
                                lower = quantile(count, 0.05, na.rm = TRUE),
                                upper = quantile(count, 0.95, na.rm = TRUE))
   
# Plot average counts across sites per week
ggplot(count_summary, aes(occasion, mean))+
                geom_point(size = 2)+
                geom_errorbar(aes(ymin = lower, ymax = upper))+
                geom_line()+
                ylab("Count")+
                xlab("Occasion (week)")+
                theme_classic()+
                theme(text = element_text(size = 18))+
                scale_x_continuous(breaks = seq(0,25,5)) 

# Fit the GAI model with several different options as shown in the paper

# 1. Fit for B = 1 brood, Poisson distribution
starting_guesses1 <- list(mu = 13, sigma = 2)
transformed_starting_guesses1 <- transform_starting_values(starting_values = starting_guesses1,
                                                          a_choice = "mixture", 
                                                          dist_choice = "P",
                                                          options = list(B = 1),
                                                          DF = smallcopper_data)
GAI_B1_P <- fit_GAI(start = transformed_starting_guesses1,
                          DF = smallcopper_data, 
                          a_choice = "mixture",
                          dist_choice = "P", 
                          options = list(B = 1),
                          verbose = T, hessian = T)

# 2. Fit for B = 2 broods with shared sigma, Poisson distribution
starting_guesses2 <- list(mu = c(7, 20), sigma = 2, w = c(.2, .8), shared_sigma = TRUE)
transformed_starting_guesses2 <- transform_starting_values(starting_values = starting_guesses2,
                                                          a_choice = "mixture", 
                                                          dist_choice = "P",
                                                          options = list(B = 2, shared_sigma = TRUE),
                                                          DF = smallcopper_data)
GAI_B2_P <- fit_GAI(start = transformed_starting_guesses2,
                          DF = smallcopper_data, 
                          a_choice = "mixture",
                          dist_choice = "P", 
                          options = list(B = 2, shared_sigma = TRUE),
                          verbose = T, hessian = T)


# 3. Fit for B = 3 broods with shared sigma, Poisson distribution
starting_guesses3 <- list(mu = c(8, 17, 25), sigma = 2, w = c(.3, .3, .4), shared_sigma = TRUE)
transformed_starting_guesses3 <- transform_starting_values(starting_values = starting_guesses3,
                                                          a_choice = "mixture", 
                                                          dist_choice = "P",
                                                          options = list(B = 3, shared_sigma = TRUE),
                                                          DF = smallcopper_data)
GAI_B3_P <- fit_GAI(start = transformed_starting_guesses3,
                           DF = smallcopper_data, 
                           a_choice = "mixture",
                           dist_choice = "P", options = list(B = 3, shared_sigma = TRUE),
                           verbose = T, hessian = T)



# 4. Fit for B = 3 broods with shared sigma, Zero-inflated Poisson distribution
starting_guesses4 <- list(mu = c(8, 17, 25), sigma = 2, w = c(.3, .3, .4), shared_sigma = TRUE)
transformed_starting_guesses4 <- transform_starting_values(starting_values = starting_guesses4,
                                                          a_choice = "mixture", 
                                                          dist_choice = "ZIP",
                                                          options = list(B = 3, shared_sigma = TRUE),
                                                          DF = smallcopper_data)
GAI_B3_ZIP <- fit_GAI(start = transformed_starting_guesses4,
                           DF = smallcopper_data, 
                           a_choice = "mixture",
                           dist_choice = "ZIP", options = list(B = 3, shared_sigma = TRUE),
                           verbose = T, hessian = T)


# 5. Fit for B = 3 broods with shared sigma, negative-binomial distribution
starting_guesses5 <- list(mu = c(8, 17, 25), sigma = 2, w = c(.3, .3, .4), shared_sigma = TRUE)
transformed_starting_guesses5 <- transform_starting_values(starting_values = starting_guesses5,
                                                          a_choice = "mixture", 
                                                          dist_choice = "NB",
                                                          options = list(B = 3, shared_sigma = TRUE),
                                                          DF = smallcopper_data)
GAI_B3_NB <- fit_GAI(start = transformed_starting_guesses5,
                               DF = smallcopper_data, 
                               a_choice = "mixture",
                               dist_choice = "NB", 
                               options = list(B = 3, shared_sigma = TRUE),
                               verbose = T, hessian = T)


# 6. Fit for B = 3 broods with separate sigma, negative-binomial distribution
starting_guesses6 <- list(mu = c(8, 17, 25), sigma = c(2, 2, 2), w= c(.3, .3, .4), shared_sigma = FALSE)
transformed_starting_guesses6 <- transform_starting_values(starting_values = starting_guesses6,
                                                          a_choice = "mixture", dist_choice = "NB",
                                                          options = list(B = 3, shared_sigma = FALSE),
                                                          DF = smallcopper_data)
GAI_B3_NB_sig <- fit_GAI(start = transformed_starting_guesses6,
                              DF = smallcopper_data, a_choice = "mixture",
                              dist_choice = "NB", options = list(B = 3, shared_sigma = FALSE),
                              verbose = T, hessian = T)

# Start modification for demonstration ----
# Tyson adds 4 broods
# 6. Fit for B = 4 broods with separate sigma, negative-binomial distribution
starting_guesses7 <- list(mu = c(8, 17, 20, 25), sigma = c(2, 2, 2, 2), w= c(.25, .25, .25, .25), shared_sigma = FALSE)
transformed_starting_guesses7 <- transform_starting_values(starting_values = starting_guesses7,
                                                           a_choice = "mixture", dist_choice = "NB",
                                                           options = list(B = 4, shared_sigma = FALSE),
                                                           DF = smallcopper_data)
GAI_B4_NB_sig <- fit_GAI(start = transformed_starting_guesses7,
                         DF = smallcopper_data, a_choice = "mixture",
                         dist_choice = "NB", options = list(B = 4, shared_sigma = FALSE),
                         verbose = T, hessian = T)

# Tyson adds 4 broods
# 6. Fit for B = 4 broods with same sigma, negative-binomial distribution
starting_guesses7 <- list(mu = c(8, 17, 20, 25), sigma = 2, w= c(.25, .25, .25, .25), shared_sigma = TRUE)
transformed_starting_guesses7 <- transform_starting_values(starting_values = starting_guesses7,
                                                           a_choice = "mixture", dist_choice = "NB",
                                                           options = list(B = 4, shared_sigma = TRUE),
                                                           DF = smallcopper_data)
GAI_B4_NB <- fit_GAI(start = transformed_starting_guesses7,
                         DF = smallcopper_data, a_choice = "mixture",
                         dist_choice = "NB", options = list(B = 4, shared_sigma = TRUE),
                         verbose = T, hessian = T)

# Tyson adds 5 broods
# Additional broods have lower AIC but are biologically meaningless
# 8. Fit for B = 5 broods with separate sigma, negative-binomial distribution
starting_guesses8 <- list(mu = c(5, 8, 17, 20, 25), sigma = c(2, 2, 2, 2, 2), w= c(.2, .2, .2, .2, .2), shared_sigma = FALSE)
transformed_starting_guesses8 <- transform_starting_values(starting_values = starting_guesses8,
                                                           a_choice = "mixture", dist_choice = "NB",
                                                           options = list(B = 5, shared_sigma = FALSE),
                                                           DF = smallcopper_data)
GAI_B5_NB_sig <- fit_GAI(start = transformed_starting_guesses8,
                         DF = smallcopper_data, a_choice = "mixture",
                         dist_choice = "NB", options = list(B = 5, shared_sigma = FALSE),
                         verbose = T, hessian = T)

modelcomparison <- data.frame(mod = 1:9,
                              npar = c(length(GAI_B1_P$par),
                                       length(GAI_B2_P$par),
                                       length(GAI_B3_P$par),
                                       length(GAI_B3_ZIP$par),
                                       length(GAI_B3_NB$par),
                                       length(GAI_B3_NB_sig$par),
                                       length(GAI_B4_NB$par),
                                       length(GAI_B4_NB_sig$par),
                                       length(GAI_B5_NB_sig$par)),
                              AIC = c(AIC(GAI_B1_P),
                                      AIC(GAI_B2_P),
                                      AIC(GAI_B3_P),
                                      AIC(GAI_B3_ZIP),
                                      AIC(GAI_B3_NB),
                                      AIC(GAI_B3_NB_sig),
                                      AIC(GAI_B4_NB),
                                      AIC(GAI_B4_NB_sig),
                                      AIC(GAI_B5_NB_sig)))

modelcomparison$deltaAIC <- modelcomparison$AIC - min(modelcomparison$AIC)



# Tyson plot distributions ----
# library(RColorBrewer)
# brewer.pal(n=5,"Dark2")
# # [1] "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E"


# best model from paper
u <- seq(from = 0,
         to = 30,
         length.out = 1e+5)
v1 <- dnorm(x = u,
            mean = 9,
            sd = 1.9)
v2 <- dnorm(x = u,
            mean = 17.4,
            sd = 1.71)
v3 <- dnorm(x = u,
            mean = 25.4,
            sd = 2.6)

w1 <- 0.064
w2 <- 0.23
w3 <- 1-w1-w2

matplot(x = u,
        y = cbind(v1*w1, v2*w2, v3*w3),
        type = "l",
        lty = 1,
        col = c("#1B9E77", "#D95F02", "#7570B3" ),
        xlab = "Weeks",
        ylab = "Densities",
        main = expression(paste("NB, B = 3, ", sigma, )["1,2,3"]))


# plot with 4 generations 
# separate sigma
u <- seq(from = 0,
         to = 30,
         length.out = 1e+5)
v1 <- dnorm(x = u,
            mean = 8.95,
            sd = 1.89)
v2 <- dnorm(x = u,
            mean = 17.57,
            sd = 1.78)
v3 <- dnorm(x = u,
            mean = 24.3,
            sd = 1.92)
v4 <- dnorm(x = u,
            mean = 27.0,
            sd = .55)

w1 <- 0.055
w2 <- 0.21
w3 <- 0.40
w4 <- 1-w1-w2-w3

matplot(x = u,
        y = cbind(v1*w1, v2*w2, v3*w3, v4*w4),
        type = "l",
        lty = 1,
        col = c("#1B9E77" ,"#D95F02" ,"#7570B3", "#E7298A"),
        xlab = "Weeks",
        ylab = "Densities",
        main = expression(paste("NB, B = 4, ", sigma, )["1,2,3,4"]))

# plot with 4 generations
# shared sigma

u <- seq(from = 0,
         to = 30,
         length.out = 1e+5)
v1 <- dnorm(x = u,
            mean = 8.91,
            sd = 1.83)
v2 <- dnorm(x = u,
            mean = 17.60,
            sd = 1.83)
v3 <- dnorm(x = u,
            mean = 24.03,
            sd = 1.83)
v4 <- dnorm(x = u,
            mean = 28.5,
            sd = 1.83)

w1 <- 0.045
w2 <- 0.17
w3 <- 0.27
w4 <- 1-w1-w2-w3

matplot(x = u,
        y = cbind(v1*w1, v2*w2, v3*w3, v4*w4),
        type = "l",
        lty = 1,
        col = c("#1B9E77", "#D95F02" ,"#7570B3" ,"#E7298A"),
        xlab = "Weeks",
        ylab = "Densities",
        main = expression(paste("NB, B = 4, ", sigma, )["shared"]))

# plot with 5 generations
u <- seq(from = 0,
         to = 30,
         length.out = 1e+5)
v1 <- dnorm(x = u,
            mean = 8.76,
            sd = 1.39)
v2 <- dnorm(x = u,
            mean = 11.8,
            sd = 4.14)
v3 <- dnorm(x = u,
            mean = 17.4,
            sd = 1.57)
v4 <- dnorm(x = u,
            mean = 25.5,
            sd = 2.74)
v5 <- dnorm(x = u,
            mean = 36.4,
            sd = 8.45)

w1 <- 0.042
w2 <- 0.02
w3 <- 0.18
w4 <- 0.64
w5 <- 1-w1-w2-w3-w4

matplot(x = u,
        y = cbind(v1*w1, v2*w2, v3*w3, v4*w4, v5*w5),
        type = "l",
        lty = 1,
        col = c("#1B9E77", "#D95F02" ,"#7570B3", "#E7298A" ,"#66A61E"),
        xlab = "Weeks",
        ylab = "Densities",
        main = expression(paste("NB, B = 5, ", sigma, )["1,2,3,4,5"]))



# End modification for demonstration ----

# Collate models for comparison by AIC (Table 2 in the paper)
modelcomparison <- data.frame(mod = 1:6,
                   npar = c(length(GAI_B1_P$par),
                            length(GAI_B2_P$par),
                            length(GAI_B3_P$par),
                            length(GAI_B3_ZIP$par),
                            length(GAI_B3_NB$par),
                            length(GAI_B3_NB_sig$par)),
                   AIC = c(AIC(GAI_B1_P),
                           AIC(GAI_B2_P),
                           AIC(GAI_B3_P),
                           AIC(GAI_B3_ZIP),
                           AIC(GAI_B3_NB),
                           AIC(GAI_B3_NB_sig)))
                   
modelcomparison$deltaAIC <- modelcomparison$AIC - min(modelcomparison$AIC)

# Now we focus on the "best" model
# Model summary 
summary(GAI_B4_NB_sig)

# Plot output using the in-built function
plot(GAI_B3_NB_sig)

# Get estimates on the parameter scale 
transform_output(GAI_B4_NB_sig)

# Parametric bootstrap to get confidence intervals for all parameters, including seasonal pattern
best_GAI_boot <- bootstrap(GAI_B3_NB_sig, R = 500, refit = FALSE,
                                   alpha = 0.05)
best_GAI_boot$par

# Produce Figure 3
# Get the predicted average count per occasion (week)
best_GAI_wcount <- data.frame(Week <- 1:26,
                     avcount = transform_output(GAI_B3_NB_sig, DF = smallcopper_data[1,], provide_A=TRUE)$A[1,]*mean(GAI_B3_NB_sig$N))

# Add 5% and 95% quantiles for the predicted counts per occasion (week)
best_GAI_wcount$predcount_upper <- apply(GAI_B3_NB_sig$N*GAI_B3_NB_sig$A, 2, quantile, 0.95)
best_GAI_wcount$predcount_lower <- apply(GAI_B3_NB_sig$N*GAI_B3_NB_sig$A, 2, quantile, 0.05)

ggplot(count_summary, aes(occasion, mean))+
  geom_point(size = 2)+
  geom_errorbar(aes(ymin = lower, ymax = upper))+
  geom_line(aes(Week, avcount), data = best_GAI_wcount, col = "blue", size = 1)+
  geom_ribbon(aes(Week, avcount, ymin = predcount_lower, ymax = predcount_upper), 
              alpha = .3, fill = "blue", data = best_GAI_wcount)+
  ylab("Count")+
  xlab("Occasion (week)")+
  theme_classic()+
  theme(text = element_text(size = 18))+
  scale_x_continuous(breaks = seq(0,25,5)) 


