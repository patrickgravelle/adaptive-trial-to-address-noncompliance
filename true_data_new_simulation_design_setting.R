library(dplyr)
library(varhandle) # need this for to.dummy()
library(mvtnorm) # need this for rmvnorm()
library(truncnorm) # need this for rtruncnorm() and dtruncnorm()



# the variables set to vary
p_assumed <- c(0.35, 0.5, 0.65)

delta_mean <- c(0.5, 1)

ndat_true <- c(1000, 2000)

prop_true <- c(0.1, 0.2)

lambda_pe <- c(0.5, 0.75)

lambda_qu <- c(0.1, 0.2)

# # the variables set to vary
# p_assumed <- c(0.5)
# 
# delta_mean <- c(1)
# 
# ndat_true <- c(1000)
# 
# prop_true <- c(0.2)


# determine all combinations
varcombos_matrix <- expand.grid(p_assumed, delta_mean, ndat_true, prop_true, lambda_pe, lambda_qu)

# set variables to vector from matrix above
p_assumed <- varcombos_matrix$Var1

delta_mean <- varcombos_matrix$Var2

ndat_true <- varcombos_matrix$Var3

prop_true <- varcombos_matrix$Var4

lambda_pe <- varcombos_matrix$Var5

lambda_qu <- varcombos_matrix$Var6


# number of configurations
n_sims <- nrow(varcombos_matrix)

# ndat_true <- rep(2*c(923, 455, 270, 178, 227, 112, 66, 44, 102, 50, 30, 20,
#                      1073, 526, 311, 206, 269, 132, 78, 52, 120, 59, 35, 23),
#                  (n_sims / 24) )

cov_mat <- rep(2, n_sims)


# determine the true compliance p to be some percent of the assumed
p_true <- p_assumed*lambda_pe
# determine the true compliance q to be some percentage increased of the assumed but less than or equal to 1
q_true <- cbind((1 + lambda_qu)*p_assumed, rep(1, n_sims))
q_true <- apply(q_true, 1, FUN = min)

# compliance parameters
beta0_p_true <- qnorm(p_true)
beta1_p_true <- rep(0.25, n_sims)
beta2_p_true <- rep(0.25, n_sims)

beta0_q_true <- qnorm(q_true)
beta1_q_true <- rep(0.25, n_sims)
beta2_q_true <- rep(0.25, n_sims)


# SD of control outcomes
sigma0_true <- rep(1, n_sims)

# correlation coefficient
rho_true <- rep(0.5, n_sims)

# SD terms for the outcomes using Covariance Matrix 4
tau1_true <- rep(1, n_sims)
taudel_true <- rep(1, n_sims)
taueps_true <- rep(1, n_sims)

# Mean parameter Values
nu0_true <- rep(0, n_sims)
nu1_true <- rep(0, n_sims)
nu2_true <- rep(0, n_sims)
nu3_true <- rep(1, n_sims)
nu4_true <- rep(0.5, n_sims)

nu5_true <- rep(0, n_sims)
nu6_true <- rep(0, n_sims)
nu7_true <- rep(0, n_sims)
nu8_true <- rep(0, n_sims)

Nc_true <- rep(0.5, n_sims)


# slant parameter 
slant_true <- rep(5, n_sims)

# degrees of freedom
dof_true <- rep(5, n_sims)



# ensuring each model fits parameters to the same dataset as one another for each configuration
set.seed(617)
seeds <- round(runif(n_sims, 0, 100000))




true_data_new_simulation_Sim_A <- matrix(NA, nrow = 1, ncol = n_sims)



for (m in 1:n_sims){
  
  
  set.seed(seeds[m])
  # set.seed(Sys.time())
  
  ndat <- 1000000
  
  treatment <- sample(c(rep(0, (Nc_true[m]*ndat)), rep(1, ((1-Nc_true[m])*ndat))), replace = FALSE)
  
  design_mat <- data.frame(treatment)
  
  x0 <- rep(1, ndat)
  
  x1 <- rnorm(ndat, 0, 1)
  
  x2 <- rnorm(ndat, 0, 1)
  
  x3 <- rnorm(ndat, 1, 1)
  
  xmat_U1 <- cbind(x0, x1, x2)
  
  xmat_U2 <- cbind(x0, x1, x2)
  
  
  u1_true <- rnorm(ndat, (xmat_U1 %*% c(beta0_p_true[m], beta1_p_true[m], beta2_p_true[m])), 1)
  
  u2_true <- rnorm(ndat, (xmat_U2 %*% c(beta0_q_true[m], beta1_q_true[m], beta2_q_true[m])), 1)
  
  w10 <- ifelse(u1_true > 0, 1, 0)
  
  w11 <- (ifelse(u2_true > 0, 1, 0))
  
  mu0_true <- nu0_true[m] + w10*nu1_true[m] + w11*nu2_true[m] + x1*nu3_true[m] + x2*nu4_true[m] +
    exp(x1)*nu5_true[m] + x3*nu6_true[m]
  
  
  if (cov_mat[m] == 2){
    
    mu_y00_true <- mu0_true
    mu_y10_true <- mu0_true + w10*delta_mean[m] + x2*w10*nu7_true[m]
    mu_y11_true <- mu0_true + w11*delta_mean[m] + x2*w11*nu8_true[m]
    
    mu_vec_true <- cbind(mu_y00_true, mu_y10_true, mu_y11_true)
    
    
    sigma_cov_true <- array(NA, dim = c(3,3,ndat))
    
    sigma_cov_true[1,1, ] <- sigma0_true[m]^2
    sigma_cov_true[2,2, ] <- sigma0_true[m]^2
    sigma_cov_true[3,3, ] <- sigma0_true[m]^2
    
    sigma_cov_true[2,1, ] <- rho_true[m]*sigma0_true[m]^2
    sigma_cov_true[1,2, ] <- rho_true[m]*sigma0_true[m]^2
    
    sigma_cov_true[1,3, ] <- rho_true[m]*sigma0_true[m]^2
    sigma_cov_true[3,1, ] <- rho_true[m]*sigma0_true[m]^2
    
    
    sigma_cov_true[2,3, ] <- rho_true[m]*sigma0_true[m]^2
    sigma_cov_true[3,2, ] <- rho_true[m]*sigma0_true[m]^2
    
    
    y_true <- t(sapply(1:ndat, function(x) rmvnorm(1, mu_vec_true[x,], sigma_cov_true[ , , x]) ))
    
    potential_outcomes <- data.frame(y_true, w10, w11)
    
    names(potential_outcomes) <- c("y00", "y10", "y11", "w10", "w11")
    
  } else if (cov_mat[m] == 4){
    
    delta_true <- rnorm(ndat, mean = delta_mean[m], sd = tau1_true[m])
    
    y00 <- rnorm(ndat, mean = mu0_true, sd = sigma0_true[m])
    
    y10 <- y00 + w10*delta_true + x2*w10*nu7_true[m] + rnorm(ndat, mean = 0, sd = taudel_true[m])
    
    y11 <- y00 + w11*delta_true + x2*w11*nu8_true[m] + rnorm(ndat, mean = 0, sd = taueps_true[m])
    
    potential_outcomes <- data.frame(y00, y10, y11, w10, w11)
    
  }
  
  
  true_data_new_simulation_Sim_A[1,m] <- mean(potential_outcomes$y11 - potential_outcomes$y00)
  
  
  save(true_data_new_simulation_Sim_A, file = "true_data_new_simulation_Sim_A.Rdata")
  
  
  
  
}



