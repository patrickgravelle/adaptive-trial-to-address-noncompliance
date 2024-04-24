library(dplyr)
library(ggplot2)
library(xtable)
library(broom)


################################################################################
################################################################################

# SIMULATION A - Multiple Configurations

################################################################################
################################################################################


bracket_mat <- matrix(NA, nrow = 6, ncol = 3)

bracket_mat[,1] <- " ("

bracket_mat[,2] <- ")"

bracket_mat[,3] <- ","

bracket_mat <- data.frame(bracket_mat)

names(bracket_mat) <- c("left_brack", "right_brack", "comma")





load("ITT_vals_model_new_Sim_A_C1_estimates.Rdata")
load("ITT_vals_model_new_Sim_A_C1_upper.Rdata")
load("ITT_vals_model_new_Sim_A_C1_lower.Rdata")
load("ITT_vals_model_new_Sim_A_C1_variance.Rdata")

load("ITT_vals_model_new_Sim_A_C2_estimates.Rdata")
load("ITT_vals_model_new_Sim_A_C2_upper.Rdata")
load("ITT_vals_model_new_Sim_A_C2_lower.Rdata")
load("ITT_vals_model_new_Sim_A_C2_variance.Rdata")

load("ITT_vals_model_new_Sim_A_C3_estimates.Rdata")
load("ITT_vals_model_new_Sim_A_C3_upper.Rdata")
load("ITT_vals_model_new_Sim_A_C3_lower.Rdata")
load("ITT_vals_model_new_Sim_A_C3_variance.Rdata")

load("ITT_vals_model_new_Sim_A_C4_estimates.Rdata")
load("ITT_vals_model_new_Sim_A_C4_upper.Rdata")
load("ITT_vals_model_new_Sim_A_C4_lower.Rdata")
load("ITT_vals_model_new_Sim_A_C4_variance.Rdata")

load("true_data_new_simulation_Sim_A.Rdata")

load("ITT_vals_model_new_Sim_A_X1_estimates.Rdata")
load("ITT_vals_model_new_Sim_A_X1_variance.Rdata")
load("ITT_vals_model_new_Sim_A_X2_estimates.Rdata")
load("ITT_vals_model_new_Sim_A_X2_variance.Rdata")








# Only got first half of simulations to run because of server shutdown


ITT_vals_model_new_Sim_A_C1_estimates <- ITT_vals_model_new_Sim_A_C1_estimates[,1:48]
ITT_vals_model_new_Sim_A_C1_upper <- ITT_vals_model_new_Sim_A_C1_upper[,1:48]
ITT_vals_model_new_Sim_A_C1_lower <- ITT_vals_model_new_Sim_A_C1_lower[,1:48]
ITT_vals_model_new_Sim_A_C1_variance <- ITT_vals_model_new_Sim_A_C1_variance[,1:48]

ITT_vals_model_new_Sim_A_C2_estimates <- ITT_vals_model_new_Sim_A_C2_estimates[,1:48]
ITT_vals_model_new_Sim_A_C2_upper <- ITT_vals_model_new_Sim_A_C2_upper[,1:48]
ITT_vals_model_new_Sim_A_C2_lower <- ITT_vals_model_new_Sim_A_C2_lower[,1:48]
ITT_vals_model_new_Sim_A_C2_variance <- ITT_vals_model_new_Sim_A_C2_variance[,1:48]

ITT_vals_model_new_Sim_A_C3_estimates <- ITT_vals_model_new_Sim_A_C3_estimates[,1:48]
ITT_vals_model_new_Sim_A_C3_upper <- ITT_vals_model_new_Sim_A_C3_upper[,1:48]
ITT_vals_model_new_Sim_A_C3_lower <- ITT_vals_model_new_Sim_A_C3_lower[,1:48]
ITT_vals_model_new_Sim_A_C3_variance <- ITT_vals_model_new_Sim_A_C3_variance[,1:48]

ITT_vals_model_new_Sim_A_C4_estimates <- ITT_vals_model_new_Sim_A_C4_estimates[,1:48]
ITT_vals_model_new_Sim_A_C4_upper <- ITT_vals_model_new_Sim_A_C4_upper[,1:48]
ITT_vals_model_new_Sim_A_C4_lower <- ITT_vals_model_new_Sim_A_C4_lower[,1:48]
ITT_vals_model_new_Sim_A_C4_variance <- ITT_vals_model_new_Sim_A_C4_variance[,1:48]


true_data_new_simulation_Sim_A <- true_data_new_simulation_Sim_A[1,1:48]

ITT_vals_model_new_Sim_A_X1_estimates <- ITT_vals_model_new_Sim_A_X1_estimates[,1:48]
ITT_vals_model_new_Sim_A_X1_variance <- ITT_vals_model_new_Sim_A_X1_variance[,1:48]
ITT_vals_model_new_Sim_A_X2_estimates <- ITT_vals_model_new_Sim_A_X2_estimates[,1:48]
ITT_vals_model_new_Sim_A_X2_variance <- ITT_vals_model_new_Sim_A_X2_variance[,1:48]




# determine number of simulations and replications

nsims <- length(true_data_new_simulation_Sim_A)
nreps <- nrow(ITT_vals_model_new_Sim_A_C1_estimates)

Model_C1_coverage <- matrix(NA, nrow = 1, ncol = nsims)
Model_C2_coverage <- matrix(NA, nrow = 1, ncol = nsims)
Model_C3_coverage <- matrix(NA, nrow = 1, ncol = nsims)
Model_C4_coverage <- matrix(NA, nrow = 1, ncol = nsims)

Model_C1_variance <- matrix(NA, nrow = 1, ncol = nsims)
Model_C2_variance <- matrix(NA, nrow = 1, ncol = nsims)
Model_C3_variance <- matrix(NA, nrow = 1, ncol = nsims)
Model_C4_variance <- matrix(NA, nrow = 1, ncol = nsims)

Model_C1_bias <- matrix(NA, nrow = 1, ncol = nsims)
Model_C2_bias <- matrix(NA, nrow = 1, ncol = nsims)
Model_C3_bias <- matrix(NA, nrow = 1, ncol = nsims)
Model_C4_bias <- matrix(NA, nrow = 1, ncol = nsims)


Model_C1_mabias <- matrix(NA, nrow = 1, ncol = nsims)
Model_C2_mabias <- matrix(NA, nrow = 1, ncol = nsims)
Model_C3_mabias <- matrix(NA, nrow = 1, ncol = nsims)
Model_C4_mabias <- matrix(NA, nrow = 1, ncol = nsims)


Model_X1_coverage <- matrix(NA, nrow = 1, ncol = nsims)
Model_X1_variance <- matrix(NA, nrow = 1, ncol = nsims)
Model_X1_bias <- matrix(NA, nrow = 1, ncol = nsims)
Model_X1_mabias <- matrix(NA, nrow = 1, ncol = nsims)

Model_X2_coverage <- matrix(NA, nrow = 1, ncol = nsims)
Model_X2_variance <- matrix(NA, nrow = 1, ncol = nsims)
Model_X2_bias <- matrix(NA, nrow = 1, ncol = nsims)
Model_X2_mabias <- matrix(NA, nrow = 1, ncol = nsims)


for (n in 1:nsims){
  
  # Coverage in each configuration
  Model_C1_coverage[1,n] <- mean(ifelse((rep(true_data_new_simulation_Sim_A[n],100) <= ITT_vals_model_new_Sim_A_C1_upper[,n] & 
                                           rep(true_data_new_simulation_Sim_A[n],100) >= ITT_vals_model_new_Sim_A_C1_lower[,n]) , 1, 0))
  
  Model_C2_coverage[1,n] <- mean(ifelse((rep(true_data_new_simulation_Sim_A[n],100) <= ITT_vals_model_new_Sim_A_C2_upper[,n] & 
                                           rep(true_data_new_simulation_Sim_A[n],100) >= ITT_vals_model_new_Sim_A_C2_lower[,n]) , 1, 0))
  
  Model_C3_coverage[1,n] <- mean(ifelse((rep(true_data_new_simulation_Sim_A[n],100) <= ITT_vals_model_new_Sim_A_C3_upper[,n] & 
                                           rep(true_data_new_simulation_Sim_A[n],100) >= ITT_vals_model_new_Sim_A_C3_lower[,n]) , 1, 0))
  
  Model_C4_coverage[1,n] <- mean(ifelse((rep(true_data_new_simulation_Sim_A[n],100) <= ITT_vals_model_new_Sim_A_C4_upper[,n] & 
                                           rep(true_data_new_simulation_Sim_A[n],100) >= ITT_vals_model_new_Sim_A_C4_lower[,n]) , 1, 0))
  
  # Mean bias in each configuration
  Model_C1_bias[1,n] <- mean((rep(true_data_new_simulation_Sim_A[n],100) - ITT_vals_model_new_Sim_A_C1_estimates[,n]))
  
  Model_C2_bias[1,n] <- mean((rep(true_data_new_simulation_Sim_A[n],100) - ITT_vals_model_new_Sim_A_C2_estimates[,n]))
  
  Model_C3_bias[1,n] <- mean((rep(true_data_new_simulation_Sim_A[n],100) - ITT_vals_model_new_Sim_A_C3_estimates[,n]))
  
  Model_C4_bias[1,n] <- mean((rep(true_data_new_simulation_Sim_A[n],100) - ITT_vals_model_new_Sim_A_C4_estimates[,n]))
  
  
  # Mean absolute bias in each configuration
  Model_C1_mabias[1,n] <- mean(abs(rep(true_data_new_simulation_Sim_A[n],100) - ITT_vals_model_new_Sim_A_C1_estimates[,n]))
  
  Model_C2_mabias[1,n] <- mean(abs(rep(true_data_new_simulation_Sim_A[n],100) - ITT_vals_model_new_Sim_A_C2_estimates[,n]))
  
  Model_C3_mabias[1,n] <- mean(abs(rep(true_data_new_simulation_Sim_A[n],100) - ITT_vals_model_new_Sim_A_C3_estimates[,n]))
  
  Model_C4_mabias[1,n] <- mean(abs(rep(true_data_new_simulation_Sim_A[n],100) - ITT_vals_model_new_Sim_A_C4_estimates[,n]))
  
  
  # Mean Variance in each configuration
  Model_C1_variance[1,n] <- mean(ITT_vals_model_new_Sim_A_C1_variance[,n])
  
  Model_C2_variance[1,n] <- mean(ITT_vals_model_new_Sim_A_C2_variance[,n])
  
  Model_C3_variance[1,n] <- mean(ITT_vals_model_new_Sim_A_C3_variance[,n])
  
  Model_C4_variance[1,n] <- mean(ITT_vals_model_new_Sim_A_C4_variance[,n])
  
  
  
  # Model X1 Results
  ITT_vals_model_new_Sim_A_X1_upper <- ITT_vals_model_new_Sim_A_X1_estimates[,n] + 1.96*sqrt(ITT_vals_model_new_Sim_A_X1_variance[,n])
  
  ITT_vals_model_new_Sim_A_X1_lower <- ITT_vals_model_new_Sim_A_X1_estimates[,n] - 1.96*sqrt(ITT_vals_model_new_Sim_A_X1_variance[,n])
  
  Model_X1_coverage[1,n] <- mean(ifelse((rep(true_data_new_simulation_Sim_A[n],100) <= ITT_vals_model_new_Sim_A_X1_upper & 
                                           rep(true_data_new_simulation_Sim_A[n],100) >= ITT_vals_model_new_Sim_A_X1_lower) , 1, 0))
  
  Model_X1_bias[1,n] <- mean((rep(true_data_new_simulation_Sim_A[n],100) - ITT_vals_model_new_Sim_A_X1_estimates[,n]))
  
  Model_X1_mabias[1,n] <- mean(abs(rep(true_data_new_simulation_Sim_A[n],100) - ITT_vals_model_new_Sim_A_X1_estimates[,n]))
  
  Model_X1_variance[1,n] <- mean(ITT_vals_model_new_Sim_A_X1_variance[,n])
  
  # Model X2 Results
  ITT_vals_model_new_Sim_A_X2_upper <- ITT_vals_model_new_Sim_A_X2_estimates[,n] + 1.96*sqrt(ITT_vals_model_new_Sim_A_X2_variance[,n])
  
  ITT_vals_model_new_Sim_A_X2_lower <- ITT_vals_model_new_Sim_A_X2_estimates[,n] - 1.96*sqrt(ITT_vals_model_new_Sim_A_X2_variance[,n])
  
  Model_X2_coverage[1,n] <- mean(ifelse((rep(true_data_new_simulation_Sim_A[n],100) <= ITT_vals_model_new_Sim_A_X2_upper & 
                                           rep(true_data_new_simulation_Sim_A[n],100) >= ITT_vals_model_new_Sim_A_X2_lower) , 1, 0))
  
  Model_X2_bias[1,n] <- mean((rep(true_data_new_simulation_Sim_A[n],100) - ITT_vals_model_new_Sim_A_X2_estimates[,n]))
  
  Model_X2_mabias[1,n] <- mean(abs(rep(true_data_new_simulation_Sim_A[n],100) - ITT_vals_model_new_Sim_A_X2_estimates[,n]))
  
  Model_X2_variance[1,n] <- mean(ITT_vals_model_new_Sim_A_X2_variance[,n])
  
}





coverage_Sim_A <- c(mean(Model_C1_coverage),
                    mean(Model_C2_coverage),
                    mean(Model_C3_coverage),
                    mean(Model_C4_coverage),
                    mean(Model_X1_coverage),
                    mean(Model_X2_coverage))

coverage_Sim_A_sd <- c(sd(Model_C1_coverage),
                    sd(Model_C2_coverage),
                    sd(Model_C3_coverage),
                    sd(Model_C4_coverage),
                    sd(Model_X1_coverage),
                    sd(Model_X2_coverage))


bias_Sim_A <- c(mean(Model_C1_bias),
                    mean(Model_C2_bias),
                    mean(Model_C3_bias),
                    mean(Model_C4_bias),
                    mean(Model_X1_bias),
                    mean(Model_X2_bias))

bias_Sim_A_sd <- c(sd(Model_C1_bias),
                sd(Model_C2_bias),
                sd(Model_C3_bias),
                sd(Model_C4_bias),
                sd(Model_X1_bias),
                sd(Model_X2_bias))


mab_bias_Sim_A <- c(mean(Model_C1_mabias),
                    mean(Model_C2_mabias),
                    mean(Model_C3_mabias),
                    mean(Model_C4_mabias),
                    mean(Model_X1_mabias),
                    mean(Model_X2_mabias))

mab_bias_Sim_A_sd <- c(sd(Model_C1_mabias),
                    sd(Model_C2_mabias),
                    sd(Model_C3_mabias),
                    sd(Model_C4_mabias),
                    sd(Model_X1_mabias),
                    sd(Model_X2_mabias))


variance_Sim_A <- c(mean(Model_C1_variance),
                    mean(Model_C2_variance),
                    mean(Model_C3_variance),
                    mean(Model_C4_variance),
                    mean(Model_X1_variance),
                    mean(Model_X2_variance))

variance_Sim_A_sd <- c(sd(Model_C1_variance),
                    sd(Model_C2_variance),
                    sd(Model_C3_variance),
                    sd(Model_C4_variance),
                    sd(Model_X1_variance),
                    sd(Model_X2_variance))

models_names <- c("Indep", "EqVar", "DifVar", "DifCov", "BothTRT", "AugTRT")


Sim_A_results <- data.frame(models_names, 
                            round(coverage_Sim_A,3), round(coverage_Sim_A_sd,2), 
                            round(variance_Sim_A, 4), round(variance_Sim_A_sd, 4), 
                            round(bias_Sim_A, 3), round(bias_Sim_A_sd, 3), 
                            round(mab_bias_Sim_A, 3), round(mab_bias_Sim_A_sd, 3))

names(Sim_A_results) <- c("Model", "Coverage", "Cov_sd", "Variance", "Var_sd", "Bias", "Bias_sd", "MAB", "MAB_sd")


Sim_A_results$Cov_full <- paste0(Sim_A_results$Coverage, bracket_mat$left_brack, Sim_A_results$Cov_sd, bracket_mat$right_brack)
Sim_A_results$Var_full <- paste0(Sim_A_results$Variance, bracket_mat$left_brack, Sim_A_results$Var_sd, bracket_mat$right_brack)
Sim_A_results$Bias_full <- paste0(Sim_A_results$Bias, bracket_mat$left_brack, Sim_A_results$Bias_sd, bracket_mat$right_brack)
Sim_A_results$MAB_full <- paste0(Sim_A_results$MAB, bracket_mat$left_brack, Sim_A_results$MAB_sd, bracket_mat$right_brack)


Sim_A_results_full <- Sim_A_results %>% select(Model, Cov_full, Var_full, Bias_full, MAB_full)

names(Sim_A_results_full) <- c("Model", "Coverage", "Variance", "Bias", "MAB")


# Sim_A_results_xtable <- xtable(t(Sim_A_results_full))
# 
# print(Sim_A_results_xtable, include.rownames = T, include.colnames = F)


Sim_A_results_xtable <- xtable((Sim_A_results_full))

print(Sim_A_results_xtable, include.rownames = F, include.colnames = T)


bias_Sim_A[6]









# full_coverage_Sim_A <- data.frame((Model_C1_coverage),
#                                   (Model_C2_coverage),
#                                   (Model_C3_coverage),
#                                   (Model_C4_coverage),
#                                   (Model_X1_coverage),
#                                   (Model_X2_coverage))



full_coverage_Sim_A <- c((Model_C1_coverage),
                         (Model_C2_coverage),
                         (Model_C3_coverage),
                         (Model_C4_coverage),
                         (Model_X1_coverage),
                         (Model_X2_coverage))

# the variables set to vary
p_assumed <- c(0.35, 0.5, 0.65)

delta_mean <- c(0.5, 1)

ndat_true <- c(1000, 2000)

prop_true <- c(0.1, 0.2)

lambda_pe <- c(0.5, 0.75)


# determine all combinations
varcombos_matrix <- expand.grid(p_assumed, delta_mean, ndat_true, prop_true, lambda_pe)

# set variables to vector from matrix above
p_assumed <- rep(varcombos_matrix$Var1, 6)

delta_mean <- rep(varcombos_matrix$Var2, 6)

ndat_true <- rep(varcombos_matrix$Var3, 6)

prop_true <- rep(varcombos_matrix$Var4, 6)

lambda_pe <- rep(varcombos_matrix$Var5, 6)


total_config <- data.frame(full_coverage_Sim_A, p_assumed, delta_mean, ndat_true, prop_true, lambda_pe)


names(total_config) <- c("coverage", "p_assumed", "trt_effect", "ndat", "prop_Z10", "lambda_pe")


important_factors_full <- aov(coverage ~ p_assumed*trt_effect*ndat*prop_Z10*lambda_pe, 
                                    data = total_config)

important_factors_full <- aov(coverage ~ p_assumed + trt_effect + ndat + prop_Z10 + lambda_pe, 
                              data = total_config)

summary(important_factors_full)



# makes output clean
ancova_check <- tidy(important_factors_full)
# just keeps variables and their MSE
ancova_dat <- data.frame(ancova_check$term,ancova_check$meansq)
# sorts the variables by decreasing MSE
ancova_dat <- ancova_dat[order(ancova_dat$ancova_check.meansq, decreasing = T), ]

ancova_dat

# # Just keeps the top 5 MSE
# ancova_dat <- ancova_dat[1:5,]

















################################################################################
################################################################################

# SIMULATION B - Multiple Configurations

################################################################################
################################################################################


load("ITT_vals_model_new_Sim_B_C1_estimates.Rdata")
load("ITT_vals_model_new_Sim_B_C1_upper.Rdata")
load("ITT_vals_model_new_Sim_B_C1_lower.Rdata")
load("ITT_vals_model_new_Sim_B_C1_variance.Rdata")

load("ITT_vals_model_new_Sim_B_C2_estimates.Rdata")
load("ITT_vals_model_new_Sim_B_C2_upper.Rdata")
load("ITT_vals_model_new_Sim_B_C2_lower.Rdata")
load("ITT_vals_model_new_Sim_B_C2_variance.Rdata")

load("ITT_vals_model_new_Sim_B_C3_estimates.Rdata")
load("ITT_vals_model_new_Sim_B_C3_upper.Rdata")
load("ITT_vals_model_new_Sim_B_C3_lower.Rdata")
load("ITT_vals_model_new_Sim_B_C3_variance.Rdata")

load("ITT_vals_model_new_Sim_B_C4_estimates.Rdata")
load("ITT_vals_model_new_Sim_B_C4_upper.Rdata")
load("ITT_vals_model_new_Sim_B_C4_lower.Rdata")
load("ITT_vals_model_new_Sim_B_C4_variance.Rdata")

load("true_data_new_simulation_Sim_B.Rdata")

load("ITT_vals_model_new_Sim_B_X1_estimates.Rdata")
load("ITT_vals_model_new_Sim_B_X1_variance.Rdata")
load("ITT_vals_model_new_Sim_B_X2_estimates.Rdata")
load("ITT_vals_model_new_Sim_B_X2_variance.Rdata")


nsims <- ncol(true_data_new_simulation_Sim_B)
nreps <- nrow(ITT_vals_model_new_Sim_B_C1_estimates)

Model_C1_coverage <- matrix(NA, nrow = 1, ncol = nsims)
Model_C2_coverage <- matrix(NA, nrow = 1, ncol = nsims)
Model_C3_coverage <- matrix(NA, nrow = 1, ncol = nsims)
Model_C4_coverage <- matrix(NA, nrow = 1, ncol = nsims)

Model_C1_variance <- matrix(NA, nrow = 1, ncol = nsims)
Model_C2_variance <- matrix(NA, nrow = 1, ncol = nsims)
Model_C3_variance <- matrix(NA, nrow = 1, ncol = nsims)
Model_C4_variance <- matrix(NA, nrow = 1, ncol = nsims)

Model_C1_bias <- matrix(NA, nrow = 1, ncol = nsims)
Model_C2_bias <- matrix(NA, nrow = 1, ncol = nsims)
Model_C3_bias <- matrix(NA, nrow = 1, ncol = nsims)
Model_C4_bias <- matrix(NA, nrow = 1, ncol = nsims)


Model_C1_mabias <- matrix(NA, nrow = 1, ncol = nsims)
Model_C2_mabias <- matrix(NA, nrow = 1, ncol = nsims)
Model_C3_mabias <- matrix(NA, nrow = 1, ncol = nsims)
Model_C4_mabias <- matrix(NA, nrow = 1, ncol = nsims)


Model_X1_coverage <- matrix(NA, nrow = 1, ncol = nsims)
Model_X1_variance <- matrix(NA, nrow = 1, ncol = nsims)
Model_X1_bias <- matrix(NA, nrow = 1, ncol = nsims)
Model_X1_mabias <- matrix(NA, nrow = 1, ncol = nsims)

Model_X2_coverage <- matrix(NA, nrow = 1, ncol = nsims)
Model_X2_variance <- matrix(NA, nrow = 1, ncol = nsims)
Model_X2_bias <- matrix(NA, nrow = 1, ncol = nsims)
Model_X2_mabias <- matrix(NA, nrow = 1, ncol = nsims)

for (n in 1:nsims){
  
  # Coverage in each configuration
  Model_C1_coverage[1,n] <- mean(ifelse((rep(true_data_new_simulation_Sim_B[1,n],100) <= ITT_vals_model_new_Sim_B_C1_upper[,n] & 
                                           rep(true_data_new_simulation_Sim_B[1,n],100) >= ITT_vals_model_new_Sim_B_C1_lower[,n]) , 1, 0))
  
  Model_C2_coverage[1,n] <- mean(ifelse((rep(true_data_new_simulation_Sim_B[1,n],100) <= ITT_vals_model_new_Sim_B_C2_upper[,n] & 
                                           rep(true_data_new_simulation_Sim_B[1,n],100) >= ITT_vals_model_new_Sim_B_C2_lower[,n]) , 1, 0))
  
  Model_C3_coverage[1,n] <- mean(ifelse((rep(true_data_new_simulation_Sim_B[1,n],100) <= ITT_vals_model_new_Sim_B_C3_upper[,n] & 
                                           rep(true_data_new_simulation_Sim_B[1,n],100) >= ITT_vals_model_new_Sim_B_C3_lower[,n]) , 1, 0))
  
  Model_C4_coverage[1,n] <- mean(ifelse((rep(true_data_new_simulation_Sim_B[1,n],100) <= ITT_vals_model_new_Sim_B_C4_upper[,n] & 
                                           rep(true_data_new_simulation_Sim_B[1,n],100) >= ITT_vals_model_new_Sim_B_C4_lower[,n]) , 1, 0))
  
  # Mean absolute bias in each configuration
  Model_C1_mabias[1,n] <- mean(abs(rep(true_data_new_simulation_Sim_B[1,n],100) - ITT_vals_model_new_Sim_B_C1_estimates[,n]))
  
  Model_C2_mabias[1,n] <- mean(abs(rep(true_data_new_simulation_Sim_B[1,n],100) - ITT_vals_model_new_Sim_B_C2_estimates[,n]))
  
  Model_C3_mabias[1,n] <- mean(abs(rep(true_data_new_simulation_Sim_B[1,n],100) - ITT_vals_model_new_Sim_B_C3_estimates[,n]))
  
  Model_C4_mabias[1,n] <- mean(abs(rep(true_data_new_simulation_Sim_B[1,n],100) - ITT_vals_model_new_Sim_B_C4_estimates[,n]))
  
  # Mean bias
  Model_C1_bias[1,n] <- mean((rep(true_data_new_simulation_Sim_B[1,n],100) - ITT_vals_model_new_Sim_B_C1_estimates[,n]))
  
  Model_C2_bias[1,n] <- mean((rep(true_data_new_simulation_Sim_B[1,n],100) - ITT_vals_model_new_Sim_B_C2_estimates[,n]))
  
  Model_C3_bias[1,n] <- mean((rep(true_data_new_simulation_Sim_B[1,n],100) - ITT_vals_model_new_Sim_B_C3_estimates[,n]))
  
  Model_C4_bias[1,n] <- mean((rep(true_data_new_simulation_Sim_B[1,n],100) - ITT_vals_model_new_Sim_B_C4_estimates[,n]))
  
  
  # Mean Variance in each configuration
  Model_C1_variance[1,n] <- mean(ITT_vals_model_new_Sim_B_C1_variance[,n])
  
  Model_C2_variance[1,n] <- mean(ITT_vals_model_new_Sim_B_C2_variance[,n])
  
  Model_C3_variance[1,n] <- mean(ITT_vals_model_new_Sim_B_C3_variance[,n])
  
  Model_C4_variance[1,n] <- mean(ITT_vals_model_new_Sim_B_C4_variance[,n])
  
  
  
  # Model X1 Results
  ITT_vals_model_new_Sim_B_X1_upper <- ITT_vals_model_new_Sim_B_X1_estimates[,n] + 1.96*sqrt(ITT_vals_model_new_Sim_B_X1_variance[,n])
  
  ITT_vals_model_new_Sim_B_X1_lower <- ITT_vals_model_new_Sim_B_X1_estimates[,n] - 1.96*sqrt(ITT_vals_model_new_Sim_B_X1_variance[,n])
  
  Model_X1_coverage[1,n] <- mean(ifelse((rep(true_data_new_simulation_Sim_B[1,n],100) <= ITT_vals_model_new_Sim_B_X1_upper & 
                                           rep(true_data_new_simulation_Sim_B[1,n],100) >= ITT_vals_model_new_Sim_B_X1_lower) , 1, 0))
  
  Model_X1_mabias[1,n] <- mean(abs(rep(true_data_new_simulation_Sim_B[1,n],100) - ITT_vals_model_new_Sim_B_X1_estimates[,n]))
  
  Model_X1_bias[1,n] <- mean((rep(true_data_new_simulation_Sim_B[1,n],100) - ITT_vals_model_new_Sim_B_X1_estimates[,n]))
  
  Model_X1_variance[1,n] <- mean(ITT_vals_model_new_Sim_B_X1_variance[,n])
  
  # Model X2 Results
  ITT_vals_model_new_Sim_B_X2_upper <- ITT_vals_model_new_Sim_B_X2_estimates[,n] + 1.96*sqrt(ITT_vals_model_new_Sim_B_X2_variance[,n])
  
  ITT_vals_model_new_Sim_B_X2_lower <- ITT_vals_model_new_Sim_B_X2_estimates[,n] - 1.96*sqrt(ITT_vals_model_new_Sim_B_X2_variance[,n])
  
  Model_X2_coverage[1,n] <- mean(ifelse((rep(true_data_new_simulation_Sim_B[1,n],100) <= ITT_vals_model_new_Sim_B_X2_upper & 
                                           rep(true_data_new_simulation_Sim_B[1,n],100) >= ITT_vals_model_new_Sim_B_X2_lower) , 1, 0))
  
  Model_X2_mabias[1,n] <- mean(abs(rep(true_data_new_simulation_Sim_B[1,n],100) - ITT_vals_model_new_Sim_B_X2_estimates[,n]))
  
  Model_X2_bias[1,n] <- mean((rep(true_data_new_simulation_Sim_B[1,n],100) - ITT_vals_model_new_Sim_B_X2_estimates[,n]))
  
  Model_X2_variance[1,n] <- mean(ITT_vals_model_new_Sim_B_X2_variance[,n])
  
}





coverage_Sim_B <- c(mean(Model_C1_coverage),
                    mean(Model_C2_coverage),
                    mean(Model_C3_coverage),
                    mean(Model_C4_coverage),
                    mean(Model_X1_coverage),
                    mean(Model_X2_coverage))

coverage_Sim_B_sd <- c(sd(Model_C1_coverage),
                       sd(Model_C2_coverage),
                       sd(Model_C3_coverage),
                       sd(Model_C4_coverage),
                       sd(Model_X1_coverage),
                       sd(Model_X2_coverage))


bias_Sim_B <- c(mean(Model_C1_bias),
                mean(Model_C2_bias),
                mean(Model_C3_bias),
                mean(Model_C4_bias),
                mean(Model_X1_bias),
                mean(Model_X2_bias))

bias_Sim_B_sd <- c(sd(Model_C1_bias),
                   sd(Model_C2_bias),
                   sd(Model_C3_bias),
                   sd(Model_C4_bias),
                   sd(Model_X1_bias),
                   sd(Model_X2_bias))


mab_bias_Sim_B <- c(mean(Model_C1_mabias),
                    mean(Model_C2_mabias),
                    mean(Model_C3_mabias),
                    mean(Model_C4_mabias),
                    mean(Model_X1_mabias),
                    mean(Model_X2_mabias))

mab_bias_Sim_B_sd <- c(sd(Model_C1_mabias),
                       sd(Model_C2_mabias),
                       sd(Model_C3_mabias),
                       sd(Model_C4_mabias),
                       sd(Model_X1_mabias),
                       sd(Model_X2_mabias))


variance_Sim_B <- c(mean(Model_C1_variance),
                    mean(Model_C2_variance),
                    mean(Model_C3_variance),
                    mean(Model_C4_variance),
                    mean(Model_X1_variance),
                    mean(Model_X2_variance))

variance_Sim_B_sd <- c(sd(Model_C1_variance),
                       sd(Model_C2_variance),
                       sd(Model_C3_variance),
                       sd(Model_C4_variance),
                       sd(Model_X1_variance),
                       sd(Model_X2_variance))

models_names <- c("Indep", "EqVar", "DifVar", "DifCov", "BothTRT", "AugTRT")


Sim_B_results <- data.frame(models_names, 
                            round(coverage_Sim_B,3), round(coverage_Sim_B_sd,2), 
                            round(variance_Sim_B, 4), round(variance_Sim_B_sd, 4), 
                            round(bias_Sim_B, 3), round(bias_Sim_B_sd, 3), 
                            round(mab_bias_Sim_B, 3), round(mab_bias_Sim_B_sd, 3))

names(Sim_B_results) <- c("Model", "Coverage", "Cov_sd", "Variance", "Var_sd", "Bias", "Bias_sd", "MAB", "MAB_sd")


Sim_B_results$Cov_full <- paste0(Sim_B_results$Coverage, bracket_mat$left_brack, Sim_B_results$Cov_sd, bracket_mat$right_brack)
Sim_B_results$Var_full <- paste0(Sim_B_results$Variance, bracket_mat$left_brack, Sim_B_results$Var_sd, bracket_mat$right_brack)
Sim_B_results$Bias_full <- paste0(Sim_B_results$Bias, bracket_mat$left_brack, Sim_B_results$Bias_sd, bracket_mat$right_brack)
Sim_B_results$MAB_full <- paste0(Sim_B_results$MAB, bracket_mat$left_brack, Sim_B_results$MAB_sd, bracket_mat$right_brack)


Sim_B_results_full <- Sim_B_results %>% select(Model, Cov_full, Var_full, Bias_full, MAB_full)

names(Sim_B_results_full) <- c("Model", "Coverage", "Variance", "Bias", "MAB")


# Sim_B_results_xtable <- xtable(t(Sim_B_results_full))
# 
# print(Sim_B_results_xtable, include.rownames = T, include.colnames = F)

Sim_B_results_xtable <- xtable((Sim_B_results_full))

print(Sim_B_results_xtable, include.rownames = F, include.colnames = T)








################################################################################
################################################################################

# SIMULATION C - Multiple Configurations

################################################################################
################################################################################


load("ITT_vals_model_new_Sim_C_C1_estimates.Rdata")
load("ITT_vals_model_new_Sim_C_C1_upper.Rdata")
load("ITT_vals_model_new_Sim_C_C1_lower.Rdata")
load("ITT_vals_model_new_Sim_C_C1_variance.Rdata")

load("ITT_vals_model_new_Sim_C_C2_estimates.Rdata")
load("ITT_vals_model_new_Sim_C_C2_upper.Rdata")
load("ITT_vals_model_new_Sim_C_C2_lower.Rdata")
load("ITT_vals_model_new_Sim_C_C2_variance.Rdata")

load("ITT_vals_model_new_Sim_C_C3_estimates.Rdata")
load("ITT_vals_model_new_Sim_C_C3_upper.Rdata")
load("ITT_vals_model_new_Sim_C_C3_lower.Rdata")
load("ITT_vals_model_new_Sim_C_C3_variance.Rdata")

load("ITT_vals_model_new_Sim_C_C4_estimates.Rdata")
load("ITT_vals_model_new_Sim_C_C4_upper.Rdata")
load("ITT_vals_model_new_Sim_C_C4_lower.Rdata")
load("ITT_vals_model_new_Sim_C_C4_variance.Rdata")

load("true_data_new_simulation_Sim_C.Rdata")

load("ITT_vals_model_new_Sim_C_X1_estimates.Rdata")
load("ITT_vals_model_new_Sim_C_X1_variance.Rdata")
load("ITT_vals_model_new_Sim_C_X2_estimates.Rdata")
load("ITT_vals_model_new_Sim_C_X2_variance.Rdata")


nsims <- ncol(true_data_new_simulation_Sim_C)
nreps <- nrow(ITT_vals_model_new_Sim_C_C1_estimates)

Model_C1_coverage <- matrix(NA, nrow = 1, ncol = nsims)
Model_C2_coverage <- matrix(NA, nrow = 1, ncol = nsims)
Model_C3_coverage <- matrix(NA, nrow = 1, ncol = nsims)
Model_C4_coverage <- matrix(NA, nrow = 1, ncol = nsims)

Model_C1_variance <- matrix(NA, nrow = 1, ncol = nsims)
Model_C2_variance <- matrix(NA, nrow = 1, ncol = nsims)
Model_C3_variance <- matrix(NA, nrow = 1, ncol = nsims)
Model_C4_variance <- matrix(NA, nrow = 1, ncol = nsims)

Model_C1_bias <- matrix(NA, nrow = 1, ncol = nsims)
Model_C2_bias <- matrix(NA, nrow = 1, ncol = nsims)
Model_C3_bias <- matrix(NA, nrow = 1, ncol = nsims)
Model_C4_bias <- matrix(NA, nrow = 1, ncol = nsims)


Model_C1_mabias <- matrix(NA, nrow = 1, ncol = nsims)
Model_C2_mabias <- matrix(NA, nrow = 1, ncol = nsims)
Model_C3_mabias <- matrix(NA, nrow = 1, ncol = nsims)
Model_C4_mabias <- matrix(NA, nrow = 1, ncol = nsims)


Model_X1_coverage <- matrix(NA, nrow = 1, ncol = nsims)
Model_X1_variance <- matrix(NA, nrow = 1, ncol = nsims)
Model_X1_bias <- matrix(NA, nrow = 1, ncol = nsims)
Model_X1_mabias <- matrix(NA, nrow = 1, ncol = nsims)

Model_X2_coverage <- matrix(NA, nrow = 1, ncol = nsims)
Model_X2_variance <- matrix(NA, nrow = 1, ncol = nsims)
Model_X2_bias <- matrix(NA, nrow = 1, ncol = nsims)
Model_X2_mabias <- matrix(NA, nrow = 1, ncol = nsims)

for (n in 1:nsims){
  
  # Coverage in each configuration
  Model_C1_coverage[1,n] <- mean(ifelse((rep(true_data_new_simulation_Sim_C[1,n],100) <= ITT_vals_model_new_Sim_C_C1_upper[,n] & 
                                           rep(true_data_new_simulation_Sim_C[1,n],100) >= ITT_vals_model_new_Sim_C_C1_lower[,n]) , 1, 0))
  
  Model_C2_coverage[1,n] <- mean(ifelse((rep(true_data_new_simulation_Sim_C[1,n],100) <= ITT_vals_model_new_Sim_C_C2_upper[,n] & 
                                           rep(true_data_new_simulation_Sim_C[1,n],100) >= ITT_vals_model_new_Sim_C_C2_lower[,n]) , 1, 0))
  
  Model_C3_coverage[1,n] <- mean(ifelse((rep(true_data_new_simulation_Sim_C[1,n],100) <= ITT_vals_model_new_Sim_C_C3_upper[,n] & 
                                           rep(true_data_new_simulation_Sim_C[1,n],100) >= ITT_vals_model_new_Sim_C_C3_lower[,n]) , 1, 0))
  
  Model_C4_coverage[1,n] <- mean(ifelse((rep(true_data_new_simulation_Sim_C[1,n],100) <= ITT_vals_model_new_Sim_C_C4_upper[,n] & 
                                           rep(true_data_new_simulation_Sim_C[1,n],100) >= ITT_vals_model_new_Sim_C_C4_lower[,n]) , 1, 0))
  
  # Mean absolute bias in each configuration
  Model_C1_mabias[1,n] <- mean(abs(rep(true_data_new_simulation_Sim_C[1,n],100) - ITT_vals_model_new_Sim_C_C1_estimates[,n]))
  
  Model_C2_mabias[1,n] <- mean(abs(rep(true_data_new_simulation_Sim_C[1,n],100) - ITT_vals_model_new_Sim_C_C2_estimates[,n]))
  
  Model_C3_mabias[1,n] <- mean(abs(rep(true_data_new_simulation_Sim_C[1,n],100) - ITT_vals_model_new_Sim_C_C3_estimates[,n]))
  
  Model_C4_mabias[1,n] <- mean(abs(rep(true_data_new_simulation_Sim_C[1,n],100) - ITT_vals_model_new_Sim_C_C4_estimates[,n]))
  
  # Mean bias
  Model_C1_bias[1,n] <- mean((rep(true_data_new_simulation_Sim_C[1,n],100) - ITT_vals_model_new_Sim_C_C1_estimates[,n]))
  
  Model_C2_bias[1,n] <- mean((rep(true_data_new_simulation_Sim_C[1,n],100) - ITT_vals_model_new_Sim_C_C2_estimates[,n]))
  
  Model_C3_bias[1,n] <- mean((rep(true_data_new_simulation_Sim_C[1,n],100) - ITT_vals_model_new_Sim_C_C3_estimates[,n]))
  
  Model_C4_bias[1,n] <- mean((rep(true_data_new_simulation_Sim_C[1,n],100) - ITT_vals_model_new_Sim_C_C4_estimates[,n]))
  
  
  # Mean Variance in each configuration
  Model_C1_variance[1,n] <- mean(ITT_vals_model_new_Sim_C_C1_variance[,n])
  
  Model_C2_variance[1,n] <- mean(ITT_vals_model_new_Sim_C_C2_variance[,n])
  
  Model_C3_variance[1,n] <- mean(ITT_vals_model_new_Sim_C_C3_variance[,n])
  
  Model_C4_variance[1,n] <- mean(ITT_vals_model_new_Sim_C_C4_variance[,n])
  
  
  
  # Model X1 Results
  ITT_vals_model_new_Sim_C_X1_upper <- ITT_vals_model_new_Sim_C_X1_estimates[,n] + 1.96*sqrt(ITT_vals_model_new_Sim_C_X1_variance[,n])
  
  ITT_vals_model_new_Sim_C_X1_lower <- ITT_vals_model_new_Sim_C_X1_estimates[,n] - 1.96*sqrt(ITT_vals_model_new_Sim_C_X1_variance[,n])
  
  Model_X1_coverage[1,n] <- mean(ifelse((rep(true_data_new_simulation_Sim_C[1,n],100) <= ITT_vals_model_new_Sim_C_X1_upper & 
                                           rep(true_data_new_simulation_Sim_C[1,n],100) >= ITT_vals_model_new_Sim_C_X1_lower) , 1, 0))
  
  Model_X1_mabias[1,n] <- mean(abs(rep(true_data_new_simulation_Sim_C[1,n],100) - ITT_vals_model_new_Sim_C_X1_estimates[,n]))
  
  Model_X1_bias[1,n] <- mean((rep(true_data_new_simulation_Sim_C[1,n],100) - ITT_vals_model_new_Sim_C_X1_estimates[,n]))
  
  Model_X1_variance[1,n] <- mean(ITT_vals_model_new_Sim_C_X1_variance[,n])
  
  # Model X2 Results
  ITT_vals_model_new_Sim_C_X2_upper <- ITT_vals_model_new_Sim_C_X2_estimates[,n] + 1.96*sqrt(ITT_vals_model_new_Sim_C_X2_variance[,n])
  
  ITT_vals_model_new_Sim_C_X2_lower <- ITT_vals_model_new_Sim_C_X2_estimates[,n] - 1.96*sqrt(ITT_vals_model_new_Sim_C_X2_variance[,n])
  
  Model_X2_coverage[1,n] <- mean(ifelse((rep(true_data_new_simulation_Sim_C[1,n],100) <= ITT_vals_model_new_Sim_C_X2_upper & 
                                           rep(true_data_new_simulation_Sim_C[1,n],100) >= ITT_vals_model_new_Sim_C_X2_lower) , 1, 0))
  
  Model_X2_mabias[1,n] <- mean(abs(rep(true_data_new_simulation_Sim_C[1,n],100) - ITT_vals_model_new_Sim_C_X2_estimates[,n]))
  
  Model_X2_bias[1,n] <- mean((rep(true_data_new_simulation_Sim_C[1,n],100) - ITT_vals_model_new_Sim_C_X2_estimates[,n]))
  
  Model_X2_variance[1,n] <- mean(ITT_vals_model_new_Sim_C_X2_variance[,n])
  
}





coverage_Sim_C <- c(mean(Model_C1_coverage),
                    mean(Model_C2_coverage),
                    mean(Model_C3_coverage),
                    mean(Model_C4_coverage),
                    mean(Model_X1_coverage),
                    mean(Model_X2_coverage))

coverage_Sim_C_sd <- c(sd(Model_C1_coverage),
                       sd(Model_C2_coverage),
                       sd(Model_C3_coverage),
                       sd(Model_C4_coverage),
                       sd(Model_X1_coverage),
                       sd(Model_X2_coverage))


bias_Sim_C <- c(mean(Model_C1_bias),
                mean(Model_C2_bias),
                mean(Model_C3_bias),
                mean(Model_C4_bias),
                mean(Model_X1_bias),
                mean(Model_X2_bias))

bias_Sim_C_sd <- c(sd(Model_C1_bias),
                   sd(Model_C2_bias),
                   sd(Model_C3_bias),
                   sd(Model_C4_bias),
                   sd(Model_X1_bias),
                   sd(Model_X2_bias))


mab_bias_Sim_C <- c(mean(Model_C1_mabias),
                    mean(Model_C2_mabias),
                    mean(Model_C3_mabias),
                    mean(Model_C4_mabias),
                    mean(Model_X1_mabias),
                    mean(Model_X2_mabias))

mab_bias_Sim_C_sd <- c(sd(Model_C1_mabias),
                       sd(Model_C2_mabias),
                       sd(Model_C3_mabias),
                       sd(Model_C4_mabias),
                       sd(Model_X1_mabias),
                       sd(Model_X2_mabias))


variance_Sim_C <- c(mean(Model_C1_variance),
                    mean(Model_C2_variance),
                    mean(Model_C3_variance),
                    mean(Model_C4_variance),
                    mean(Model_X1_variance),
                    mean(Model_X2_variance))

variance_Sim_C_sd <- c(sd(Model_C1_variance),
                       sd(Model_C2_variance),
                       sd(Model_C3_variance),
                       sd(Model_C4_variance),
                       sd(Model_X1_variance),
                       sd(Model_X2_variance))

models_names <- c("Indep", "EqVar", "DifVar", "DifCov", "BothTRT", "AugTRT")


Sim_C_results <- data.frame(models_names, 
                            round(coverage_Sim_C,3), round(coverage_Sim_C_sd,2), 
                            round(variance_Sim_C, 4), round(variance_Sim_C_sd, 4), 
                            round(bias_Sim_C, 3), round(bias_Sim_C_sd, 3), 
                            round(mab_bias_Sim_C, 3), round(mab_bias_Sim_C_sd, 3))

names(Sim_C_results) <- c("Model", "Coverage", "Cov_sd", "Variance", "Var_sd", "Bias", "Bias_sd", "MAB", "MAB_sd")


Sim_C_results$Cov_full <- paste0(Sim_C_results$Coverage, bracket_mat$left_brack, Sim_C_results$Cov_sd, bracket_mat$right_brack)
Sim_C_results$Var_full <- paste0(Sim_C_results$Variance, bracket_mat$left_brack, Sim_C_results$Var_sd, bracket_mat$right_brack)
Sim_C_results$Bias_full <- paste0(Sim_C_results$Bias, bracket_mat$left_brack, Sim_C_results$Bias_sd, bracket_mat$right_brack)
Sim_C_results$MAB_full <- paste0(Sim_C_results$MAB, bracket_mat$left_brack, Sim_C_results$MAB_sd, bracket_mat$right_brack)


Sim_C_results_full <- Sim_C_results %>% select(Model, Cov_full, Var_full, Bias_full, MAB_full)

names(Sim_C_results_full) <- c("Model", "Coverage", "Variance", "Bias", "MAB")


# Sim_C_results_xtable <- xtable(t(Sim_C_results_full))
# 
# print(Sim_C_results_xtable, include.rownames = T, include.colnames = F)

Sim_C_results_xtable <- xtable((Sim_C_results_full))

print(Sim_C_results_xtable, include.rownames = F, include.colnames = T)








################################################################################
################################################################################

# SIMULATION D - Multiple Configurations

################################################################################
################################################################################


load("ITT_vals_model_new_Sim_D_C1_estimates.Rdata")
load("ITT_vals_model_new_Sim_D_C1_upper.Rdata")
load("ITT_vals_model_new_Sim_D_C1_lower.Rdata")
load("ITT_vals_model_new_Sim_D_C1_variance.Rdata")

load("ITT_vals_model_new_Sim_D_C2_estimates.Rdata")
load("ITT_vals_model_new_Sim_D_C2_upper.Rdata")
load("ITT_vals_model_new_Sim_D_C2_lower.Rdata")
load("ITT_vals_model_new_Sim_D_C2_variance.Rdata")

load("ITT_vals_model_new_Sim_D_C3_estimates.Rdata")
load("ITT_vals_model_new_Sim_D_C3_upper.Rdata")
load("ITT_vals_model_new_Sim_D_C3_lower.Rdata")
load("ITT_vals_model_new_Sim_D_C3_variance.Rdata")

load("ITT_vals_model_new_Sim_D_C4_estimates.Rdata")
load("ITT_vals_model_new_Sim_D_C4_upper.Rdata")
load("ITT_vals_model_new_Sim_D_C4_lower.Rdata")
load("ITT_vals_model_new_Sim_D_C4_variance.Rdata")

load("true_data_new_simulation_Sim_D.Rdata")

load("ITT_vals_model_new_Sim_D_X1_estimates.Rdata")
load("ITT_vals_model_new_Sim_D_X1_variance.Rdata")
load("ITT_vals_model_new_Sim_D_X2_estimates.Rdata")
load("ITT_vals_model_new_Sim_D_X2_variance.Rdata")


nsims <- ncol(true_data_new_simulation_Sim_D)
nreps <- nrow(ITT_vals_model_new_Sim_D_C1_estimates)

Model_C1_coverage <- matrix(NA, nrow = 1, ncol = nsims)
Model_C2_coverage <- matrix(NA, nrow = 1, ncol = nsims)
Model_C3_coverage <- matrix(NA, nrow = 1, ncol = nsims)
Model_C4_coverage <- matrix(NA, nrow = 1, ncol = nsims)

Model_C1_variance <- matrix(NA, nrow = 1, ncol = nsims)
Model_C2_variance <- matrix(NA, nrow = 1, ncol = nsims)
Model_C3_variance <- matrix(NA, nrow = 1, ncol = nsims)
Model_C4_variance <- matrix(NA, nrow = 1, ncol = nsims)

Model_C1_bias <- matrix(NA, nrow = 1, ncol = nsims)
Model_C2_bias <- matrix(NA, nrow = 1, ncol = nsims)
Model_C3_bias <- matrix(NA, nrow = 1, ncol = nsims)
Model_C4_bias <- matrix(NA, nrow = 1, ncol = nsims)


Model_C1_mabias <- matrix(NA, nrow = 1, ncol = nsims)
Model_C2_mabias <- matrix(NA, nrow = 1, ncol = nsims)
Model_C3_mabias <- matrix(NA, nrow = 1, ncol = nsims)
Model_C4_mabias <- matrix(NA, nrow = 1, ncol = nsims)


Model_X1_coverage <- matrix(NA, nrow = 1, ncol = nsims)
Model_X1_variance <- matrix(NA, nrow = 1, ncol = nsims)
Model_X1_bias <- matrix(NA, nrow = 1, ncol = nsims)
Model_X1_mabias <- matrix(NA, nrow = 1, ncol = nsims)

Model_X2_coverage <- matrix(NA, nrow = 1, ncol = nsims)
Model_X2_variance <- matrix(NA, nrow = 1, ncol = nsims)
Model_X2_bias <- matrix(NA, nrow = 1, ncol = nsims)
Model_X2_mabias <- matrix(NA, nrow = 1, ncol = nsims)

for (n in 1:nsims){
  
  # Coverage in each configuration
  Model_C1_coverage[1,n] <- mean(ifelse((rep(true_data_new_simulation_Sim_D[1,n],100) <= ITT_vals_model_new_Sim_D_C1_upper[,n] & 
                                           rep(true_data_new_simulation_Sim_D[1,n],100) >= ITT_vals_model_new_Sim_D_C1_lower[,n]) , 1, 0))
  
  Model_C2_coverage[1,n] <- mean(ifelse((rep(true_data_new_simulation_Sim_D[1,n],100) <= ITT_vals_model_new_Sim_D_C2_upper[,n] & 
                                           rep(true_data_new_simulation_Sim_D[1,n],100) >= ITT_vals_model_new_Sim_D_C2_lower[,n]) , 1, 0))
  
  Model_C3_coverage[1,n] <- mean(ifelse((rep(true_data_new_simulation_Sim_D[1,n],100) <= ITT_vals_model_new_Sim_D_C3_upper[,n] & 
                                           rep(true_data_new_simulation_Sim_D[1,n],100) >= ITT_vals_model_new_Sim_D_C3_lower[,n]) , 1, 0))
  
  Model_C4_coverage[1,n] <- mean(ifelse((rep(true_data_new_simulation_Sim_D[1,n],100) <= ITT_vals_model_new_Sim_D_C4_upper[,n] & 
                                           rep(true_data_new_simulation_Sim_D[1,n],100) >= ITT_vals_model_new_Sim_D_C4_lower[,n]) , 1, 0))
  
  # Mean absolute bias in each configuration
  Model_C1_mabias[1,n] <- mean(abs(rep(true_data_new_simulation_Sim_D[1,n],100) - ITT_vals_model_new_Sim_D_C1_estimates[,n]))
  
  Model_C2_mabias[1,n] <- mean(abs(rep(true_data_new_simulation_Sim_D[1,n],100) - ITT_vals_model_new_Sim_D_C2_estimates[,n]))
  
  Model_C3_mabias[1,n] <- mean(abs(rep(true_data_new_simulation_Sim_D[1,n],100) - ITT_vals_model_new_Sim_D_C3_estimates[,n]))
  
  Model_C4_mabias[1,n] <- mean(abs(rep(true_data_new_simulation_Sim_D[1,n],100) - ITT_vals_model_new_Sim_D_C4_estimates[,n]))
  
  # Mean bias
  Model_C1_bias[1,n] <- mean((rep(true_data_new_simulation_Sim_D[1,n],100) - ITT_vals_model_new_Sim_D_C1_estimates[,n]))
  
  Model_C2_bias[1,n] <- mean((rep(true_data_new_simulation_Sim_D[1,n],100) - ITT_vals_model_new_Sim_D_C2_estimates[,n]))
  
  Model_C3_bias[1,n] <- mean((rep(true_data_new_simulation_Sim_D[1,n],100) - ITT_vals_model_new_Sim_D_C3_estimates[,n]))
  
  Model_C4_bias[1,n] <- mean((rep(true_data_new_simulation_Sim_D[1,n],100) - ITT_vals_model_new_Sim_D_C4_estimates[,n]))
  
  
  # Mean Variance in each configuration
  Model_C1_variance[1,n] <- mean(ITT_vals_model_new_Sim_D_C1_variance[,n])
  
  Model_C2_variance[1,n] <- mean(ITT_vals_model_new_Sim_D_C2_variance[,n])
  
  Model_C3_variance[1,n] <- mean(ITT_vals_model_new_Sim_D_C3_variance[,n])
  
  Model_C4_variance[1,n] <- mean(ITT_vals_model_new_Sim_D_C4_variance[,n])
  
  
  
  # Model X1 Results
  ITT_vals_model_new_Sim_D_X1_upper <- ITT_vals_model_new_Sim_D_X1_estimates[,n] + 1.96*sqrt(ITT_vals_model_new_Sim_D_X1_variance[,n])
  
  ITT_vals_model_new_Sim_D_X1_lower <- ITT_vals_model_new_Sim_D_X1_estimates[,n] - 1.96*sqrt(ITT_vals_model_new_Sim_D_X1_variance[,n])
  
  Model_X1_coverage[1,n] <- mean(ifelse((rep(true_data_new_simulation_Sim_D[1,n],100) <= ITT_vals_model_new_Sim_D_X1_upper & 
                                           rep(true_data_new_simulation_Sim_D[1,n],100) >= ITT_vals_model_new_Sim_D_X1_lower) , 1, 0))
  
  Model_X1_mabias[1,n] <- mean(abs(rep(true_data_new_simulation_Sim_D[1,n],100) - ITT_vals_model_new_Sim_D_X1_estimates[,n]))
  
  Model_X1_bias[1,n] <- mean((rep(true_data_new_simulation_Sim_D[1,n],100) - ITT_vals_model_new_Sim_D_X1_estimates[,n]))
  
  Model_X1_variance[1,n] <- mean(ITT_vals_model_new_Sim_D_X1_variance[,n])
  
  # Model X2 Results
  ITT_vals_model_new_Sim_D_X2_upper <- ITT_vals_model_new_Sim_D_X2_estimates[,n] + 1.96*sqrt(ITT_vals_model_new_Sim_D_X2_variance[,n])
  
  ITT_vals_model_new_Sim_D_X2_lower <- ITT_vals_model_new_Sim_D_X2_estimates[,n] - 1.96*sqrt(ITT_vals_model_new_Sim_D_X2_variance[,n])
  
  Model_X2_coverage[1,n] <- mean(ifelse((rep(true_data_new_simulation_Sim_D[1,n],100) <= ITT_vals_model_new_Sim_D_X2_upper & 
                                           rep(true_data_new_simulation_Sim_D[1,n],100) >= ITT_vals_model_new_Sim_D_X2_lower) , 1, 0))
  
  Model_X2_mabias[1,n] <- mean(abs(rep(true_data_new_simulation_Sim_D[1,n],100) - ITT_vals_model_new_Sim_D_X2_estimates[,n]))
  
  Model_X2_bias[1,n] <- mean((rep(true_data_new_simulation_Sim_D[1,n],100) - ITT_vals_model_new_Sim_D_X2_estimates[,n]))
  
  Model_X2_variance[1,n] <- mean(ITT_vals_model_new_Sim_D_X2_variance[,n])
  
}





coverage_Sim_D <- c(mean(Model_C1_coverage),
                    mean(Model_C2_coverage),
                    mean(Model_C3_coverage),
                    mean(Model_C4_coverage),
                    mean(Model_X1_coverage),
                    mean(Model_X2_coverage))

coverage_Sim_D_sd <- c(sd(Model_C1_coverage),
                       sd(Model_C2_coverage),
                       sd(Model_C3_coverage),
                       sd(Model_C4_coverage),
                       sd(Model_X1_coverage),
                       sd(Model_X2_coverage))


bias_Sim_D <- c(mean(Model_C1_bias),
                mean(Model_C2_bias),
                mean(Model_C3_bias),
                mean(Model_C4_bias),
                mean(Model_X1_bias),
                mean(Model_X2_bias))

bias_Sim_D_sd <- c(sd(Model_C1_bias),
                   sd(Model_C2_bias),
                   sd(Model_C3_bias),
                   sd(Model_C4_bias),
                   sd(Model_X1_bias),
                   sd(Model_X2_bias))


mab_bias_Sim_D <- c(mean(Model_C1_mabias),
                    mean(Model_C2_mabias),
                    mean(Model_C3_mabias),
                    mean(Model_C4_mabias),
                    mean(Model_X1_mabias),
                    mean(Model_X2_mabias))

mab_bias_Sim_D_sd <- c(sd(Model_C1_mabias),
                       sd(Model_C2_mabias),
                       sd(Model_C3_mabias),
                       sd(Model_C4_mabias),
                       sd(Model_X1_mabias),
                       sd(Model_X2_mabias))


variance_Sim_D <- c(mean(Model_C1_variance),
                    mean(Model_C2_variance),
                    mean(Model_C3_variance),
                    mean(Model_C4_variance),
                    mean(Model_X1_variance),
                    mean(Model_X2_variance))

variance_Sim_D_sd <- c(sd(Model_C1_variance),
                       sd(Model_C2_variance),
                       sd(Model_C3_variance),
                       sd(Model_C4_variance),
                       sd(Model_X1_variance),
                       sd(Model_X2_variance))

models_names <- c("Indep", "EqVar", "DifVar", "DifCov", "BothTRT", "AugTRT")


Sim_D_results <- data.frame(models_names, 
                            round(coverage_Sim_D,3), round(coverage_Sim_D_sd,2), 
                            round(variance_Sim_D, 4), round(variance_Sim_D_sd, 4), 
                            round(bias_Sim_D, 3), round(bias_Sim_D_sd, 3), 
                            round(mab_bias_Sim_D, 3), round(mab_bias_Sim_D_sd, 3))

names(Sim_D_results) <- c("Model", "Coverage", "Cov_sd", "Variance", "Var_sd", "Bias", "Bias_sd", "MAB", "MAB_sd")


Sim_D_results$Cov_full <- paste0(Sim_D_results$Coverage, bracket_mat$left_brack, Sim_D_results$Cov_sd, bracket_mat$right_brack)
Sim_D_results$Var_full <- paste0(Sim_D_results$Variance, bracket_mat$left_brack, Sim_D_results$Var_sd, bracket_mat$right_brack)
Sim_D_results$Bias_full <- paste0(Sim_D_results$Bias, bracket_mat$left_brack, Sim_D_results$Bias_sd, bracket_mat$right_brack)
Sim_D_results$MAB_full <- paste0(Sim_D_results$MAB, bracket_mat$left_brack, Sim_D_results$MAB_sd, bracket_mat$right_brack)


Sim_D_results_full <- Sim_D_results %>% select(Model, Cov_full, Var_full, Bias_full, MAB_full)

names(Sim_D_results_full) <- c("Model", "Coverage", "Variance", "Bias", "MAB")


# Sim_D_results_xtable <- xtable(t(Sim_D_results_full))
# 
# print(Sim_D_results_xtable, include.rownames = T, include.colnames = F)

Sim_D_results_xtable <- xtable((Sim_D_results_full))

print(Sim_D_results_xtable, include.rownames = F, include.colnames = T)




