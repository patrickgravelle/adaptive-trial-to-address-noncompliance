library(dplyr)
library(varhandle) # need this for to.dummy()
library(mvtnorm) # need this for rmvnorm()
library(truncnorm) # need this for rtruncnorm() and dtruncnorm()



# p and q prior values
alpha_p <- 1
beta_p <- 1
alpha_q <- 1
beta_q <- 1

# mu0 prior terms

k0_mu0 <- 100
k1_mu0 <- 0


log_lik_fnct_full <- function(theta, x1, x2, y00, y10, y11, w10, w11){
  
  # write out these definitions
  theta_1 <- (y00 - theta[3] - w10*theta[4] - w11*theta[5] - x1*theta[6] - x2*theta[7] )^2
  
  theta_2 <- (y10 - theta[3] - w10*theta[4] - w11*theta[5] - x1*theta[6] - x2*theta[7] - w10*theta[8])^2
  
  theta_3 <- (y11 - theta[3] - w10*theta[4] - w11*theta[5] - x1*theta[6] - x2*theta[7] - w11*theta[8])^2
  
  theta_4 <- 2*(y00 - theta[3] - w10*theta[4] - w11*theta[5] - x1*theta[6] - x2*theta[7])*(y10 - theta[3] - w10*theta[4] - w11*theta[5] - x1*theta[6] - x2*theta[7] - w10*theta[8])
  
  theta_5 <- 2*(y00 - theta[3] - w10*theta[4] - w11*theta[5] - x1*theta[6] - x2*theta[7])*(y11 - theta[3] - w10*theta[4] - w11*theta[5] - x1*theta[6] - x2*theta[7] - w11*theta[8])
  
  theta_6 <- 2*(y10 - theta[3] - w10*theta[4] - w11*theta[5] - x1*theta[6] - x2*theta[7] - w10*theta[8])*(y11 - theta[3] - w10*theta[4] - w11*theta[5] - x1*theta[6] - x2*theta[7] - w11*theta[8])
  
  
  sigma11 <- theta[1]^2
  
  sigma12 <- theta[2]*theta[1]^2
  
  sigma13 <- theta[2]*theta[1]^2
  
  sigma23 <- theta[2]*theta[1]^2
  
  sigma22 <- theta[1]^2
  
  sigma33 <- theta[1]^2
  
  theta_big <- theta_1 * (sigma22*sigma33 - sigma23^2) + theta_2 * (sigma11*sigma33 - sigma13^2) + theta_3 * (sigma11*sigma22 - sigma12^2) +
    theta_4 * (sigma23*sigma13 - sigma12*sigma33) + theta_5 * (sigma12*sigma23 - sigma13*sigma22) + theta_6 * (sigma12*sigma13 - sigma11*sigma23)
  
  det_Z <- sigma11*sigma22*sigma33 + 2*sigma12*sigma13*sigma23 - sigma12^2 * sigma33 - sigma13^2 * sigma22 - sigma23^2 * sigma11
  
  logpost <- sum((-0.5)*log(det_Z) - 0.5*(theta_big / det_Z))
  
  return(logpost)
  
}



log_lik_fnct_WS <- function(theta, x1, x2, w10, w11){
  
  logpost <- sum( log( (pnorm((theta[4] + x1*theta[5] + x2*theta[6])))^w11 * (1-(pnorm((theta[4] + x1*theta[5] + x2*theta[6]))))^(1-w11) ) + log( (pnorm((theta[1] + x1*theta[2] + x2*theta[3])))^w10 * (1-(pnorm((theta[1] + x1*theta[2] + x2*theta[3]))))^(1-w10) ) )
  
  return(logpost)
  
}


# Prior of f(sigma0) ~ half-Cauchy

# Sigma0 conditional posterior

sigma0_cond_posterior <- function(sigma0, rho, mu0, delta, y00, y10, y11, w10, w11){
  
  # write out these definitions
  theta_1 <- (y00 - mu0)^2
  
  theta_2 <- (y10 - mu0 - w10*delta)^2
  
  theta_3 <- (y11 - mu0 - w11*delta)^2
  
  theta_4 <- 2*(y00 - mu0)*(y10 - mu0 - w10*delta)
  
  theta_5 <- 2*(y00 - mu0)*(y11 - mu0 - w11*delta)
  
  theta_6 <- 2*(y10 - mu0 - w10*delta)*(y11 - mu0 - w11*delta)
  
  
  sigma11 <- sigma0
  
  sigma12 <- rho*sigma0
  
  sigma13 <- rho*sigma0
  
  sigma23 <- rho*sigma0
  
  sigma22 <- sigma0
  
  sigma33 <- sigma0
  
  theta_big <- theta_1 * (sigma22*sigma33 - sigma23^2) + theta_2 * (sigma11*sigma33 - sigma13^2) + theta_3 * (sigma11*sigma22 - sigma12^2) +
    theta_4 * (sigma23*sigma13 - sigma12*sigma33) + theta_5 * (sigma12*sigma23 - sigma13*sigma22) + theta_6 * (sigma12*sigma13 - sigma11*sigma23)
  
  det_Z <- sigma11*sigma22*sigma33 + 2*sigma12*sigma13*sigma23 - sigma12^2 * sigma33 - sigma13^2 * sigma22 - sigma23^2 * sigma11
  
  logpost <- (-1) * log( 1 + ((sigma0)^2) / (25^2) ) + sum((-0.5)*log(det_Z) - 0.5*(theta_big / det_Z))
  
  return(logpost)
  
}



rho_cond_posterior <- function(rho, sigma0, mu0, delta, y00, y10, y11, w10, w11){
  
  # write out these definitions
  theta_1 <- (y00 - mu0)^2
  
  theta_2 <- (y10 - mu0 - w10*delta)^2
  
  theta_3 <- (y11 - mu0 - w11*delta)^2
  
  theta_4 <- 2*(y00 - mu0)*(y10 - mu0 - w10*delta)
  
  theta_5 <- 2*(y00 - mu0)*(y11 - mu0 - w11*delta)
  
  theta_6 <- 2*(y10 - mu0 - w10*delta)*(y11 - mu0 - w11*delta)
  
  
  sigma11 <- sigma0^2
  
  sigma12 <- rho*sigma0^2
  
  sigma13 <- rho*sigma0^2
  
  sigma23 <- rho*sigma0^2
  
  sigma22 <- sigma0^2
  
  sigma33 <- sigma0^2
  
  
  theta_big <- theta_1 * (sigma22*sigma33 - sigma23^2) + theta_2 * (sigma11*sigma33 - sigma13^2) + theta_3 * (sigma11*sigma22 - sigma12^2) +
    theta_4 * (sigma23*sigma13 - sigma12*sigma33) + theta_5 * (sigma12*sigma23 - sigma13*sigma22) + theta_6 * (sigma12*sigma13 - sigma11*sigma23)
  
  det_Z <- sigma11*sigma22*sigma33 + 2*sigma12*sigma13*sigma23 - sigma12^2 * sigma33 - sigma13^2 * sigma22 - sigma23^2 * sigma11
  
  logpost <- sum((-0.5)*log(det_Z) - 0.5*(theta_big / det_Z))
  
  # logpost <- (alpha - 1) * log(rho) + (beta - 1) * log(1 - rho)  + sum((-0.5)*log(det_Z) - 0.5*(theta_big / det_Z))
  
  return(logpost)
  
}



# nu0 conditional posterior

nu0_cond_posterior <- function(sigma0, rho, nu1, nu2, nu3, nu4, delta, x1, x2, y00, y10, y11, w10, w11, k0, k1){
  
  sigma11 <- sigma0^2
  
  sigma12 <- rho*sigma0^2
  
  sigma13 <- rho*sigma0^2
  
  sigma23 <- rho*sigma0^2
  
  sigma22 <- sigma0^2
  
  sigma33 <- sigma0^2
  
  
  gamma_0 <- sigma22*sigma33 - sigma23^2
  gamma_1 <- sigma11*sigma33 - sigma13^2
  gamma_2 <- sigma11*sigma22 - sigma12^2
  gamma_3 <- sigma23*sigma13 - sigma12*sigma33
  gamma_4 <- sigma12*sigma23 - sigma13*sigma22
  gamma_5 <- sigma12*sigma13 - sigma11*sigma23
  
  
  theta_denom <- sigma11*sigma22*sigma33 + 2*sigma12*sigma13*sigma23 - sigma12^2 * sigma33 - sigma13^2 * sigma22 - sigma23^2 * sigma11
  
  
  theta_C <- gamma_0 + gamma_1 + gamma_2 + 2*gamma_3 + 2*gamma_4 + 2*gamma_5
  
  theta_D <- 2 * ( -gamma_0*y00 - gamma_1*y10 - gamma_2*y11 - gamma_3*y00 - gamma_3*y10 - gamma_4*y00 -gamma_4*y11 - gamma_5*y10 - gamma_5*y11 +
                     gamma_1*w10*delta + gamma_2*w11*delta + gamma_3*w10*delta + gamma_4*w11*delta + gamma_5*w10*delta + gamma_5*w11*delta )
  
  theta_A <- theta_C
  
  theta_B <- 2*theta_C*(w10*nu1 + w11*nu2 + x1*nu3 + x2*nu4) + theta_D
  
  mu0_sd <- sqrt( 1 / (ndat*(theta_A / theta_denom) + 1/k0 ) )
  
  mu0_mean <- (-1) * ( sum(theta_B / theta_denom) - 2*k1/k0 ) / ( 2 * ( ndat*(theta_A / theta_denom) + 1/k0 ) )
  
  
  mu0_t_draw <- rnorm(1, mean = mu0_mean, sd = mu0_sd)
  
  return(mu0_t_draw)
  
}




# nu1 conditional posterior

nu1_cond_posterior <- function(sigma0, rho, nu0, nu2, nu3, nu4, delta, x1, x2, y00, y10, y11, w10, w11, k0, k1){
  
  sigma11 <- sigma0^2
  
  sigma12 <- rho*sigma0^2
  
  sigma13 <- rho*sigma0^2
  
  sigma23 <- rho*sigma0^2
  
  sigma22 <- sigma0^2
  
  sigma33 <- sigma0^2
  
  
  gamma_0 <- sigma22*sigma33 - sigma23^2
  gamma_1 <- sigma11*sigma33 - sigma13^2
  gamma_2 <- sigma11*sigma22 - sigma12^2
  gamma_3 <- sigma23*sigma13 - sigma12*sigma33
  gamma_4 <- sigma12*sigma23 - sigma13*sigma22
  gamma_5 <- sigma12*sigma13 - sigma11*sigma23
  
  
  theta_denom <- sigma11*sigma22*sigma33 + 2*sigma12*sigma13*sigma23 - sigma12^2 * sigma33 - sigma13^2 * sigma22 - sigma23^2 * sigma11
  
  
  theta_C <- gamma_0 + gamma_1 + gamma_2 + 2*gamma_3 + 2*gamma_4 + 2*gamma_5
  
  theta_D <- 2 * ( -gamma_0*y00 - gamma_1*y10 - gamma_2*y11 - gamma_3*y00 - gamma_3*y10 - gamma_4*y00 -gamma_4*y11 - gamma_5*y10 - gamma_5*y11 +
                     gamma_1*w10*delta + gamma_2*w11*delta + gamma_3*w10*delta + gamma_4*w11*delta + gamma_5*w10*delta + gamma_5*w11*delta )
  
  theta_A <- w10*theta_C
  
  theta_B <- 2*theta_C*w10*(nu0 + w11*nu2 + x1*nu3 + x2*nu4) + w10*theta_D
  
  mu0_sd <- sqrt( 1 / (sum(theta_A / theta_denom) + 1/k0 ) )
  
  mu0_mean <- (-1) * ( sum(theta_B / theta_denom) - 2*k1/k0 ) / ( 2 * ( sum(theta_A / theta_denom) + 1/k0 ) )
  
  
  mu0_t_draw <- rnorm(1, mean = mu0_mean, sd = mu0_sd)
  
  return(mu0_t_draw)
  
}




# nu2 conditional posterior

nu2_cond_posterior <- function(sigma0, rho, nu0, nu1, nu3, nu4, delta, x1, x2, y00, y10, y11, w10, w11, k0, k1){
  
  sigma11 <- sigma0^2
  
  sigma12 <- rho*sigma0^2
  
  sigma13 <- rho*sigma0^2
  
  sigma23 <- rho*sigma0^2
  
  sigma22 <- sigma0^2
  
  sigma33 <- sigma0^2
  
  
  gamma_0 <- sigma22*sigma33 - sigma23^2
  gamma_1 <- sigma11*sigma33 - sigma13^2
  gamma_2 <- sigma11*sigma22 - sigma12^2
  gamma_3 <- sigma23*sigma13 - sigma12*sigma33
  gamma_4 <- sigma12*sigma23 - sigma13*sigma22
  gamma_5 <- sigma12*sigma13 - sigma11*sigma23
  
  
  theta_denom <- sigma11*sigma22*sigma33 + 2*sigma12*sigma13*sigma23 - sigma12^2 * sigma33 - sigma13^2 * sigma22 - sigma23^2 * sigma11
  
  
  theta_C <- gamma_0 + gamma_1 + gamma_2 + 2*gamma_3 + 2*gamma_4 + 2*gamma_5
  
  theta_D <- 2 * ( -gamma_0*y00 - gamma_1*y10 - gamma_2*y11 - gamma_3*y00 - gamma_3*y10 - gamma_4*y00 -gamma_4*y11 - gamma_5*y10 - gamma_5*y11 +
                     gamma_1*w10*delta + gamma_2*w11*delta + gamma_3*w10*delta + gamma_4*w11*delta + gamma_5*w10*delta + gamma_5*w11*delta )
  
  theta_A <- w11*theta_C
  
  theta_B <- 2*theta_C*w11*(nu0 + w10*nu1 + x1*nu3 + x2*nu4) + w11*theta_D
  
  mu0_sd <- sqrt( 1 / (sum(theta_A / theta_denom) + 1/k0 ) )
  
  mu0_mean <- (-1) * ( sum(theta_B / theta_denom) - 2*k1/k0 ) / ( 2 * ( sum(theta_A / theta_denom) + 1/k0 ) )
  
  
  mu0_t_draw <- rnorm(1, mean = mu0_mean, sd = mu0_sd)
  
  return(mu0_t_draw)
  
}



# nu3 conditional posterior

nu3_cond_posterior <- function(sigma0, rho, nu0, nu1, nu2, nu4, delta, x1, x2, y00, y10, y11, w10, w11, k0, k1){
  
  sigma11 <- sigma0^2
  
  sigma12 <- rho*sigma0^2
  
  sigma13 <- rho*sigma0^2
  
  sigma23 <- rho*sigma0^2
  
  sigma22 <- sigma0^2
  
  sigma33 <- sigma0^2
  
  
  gamma_0 <- sigma22*sigma33 - sigma23^2
  gamma_1 <- sigma11*sigma33 - sigma13^2
  gamma_2 <- sigma11*sigma22 - sigma12^2
  gamma_3 <- sigma23*sigma13 - sigma12*sigma33
  gamma_4 <- sigma12*sigma23 - sigma13*sigma22
  gamma_5 <- sigma12*sigma13 - sigma11*sigma23
  
  
  theta_denom <- sigma11*sigma22*sigma33 + 2*sigma12*sigma13*sigma23 - sigma12^2 * sigma33 - sigma13^2 * sigma22 - sigma23^2 * sigma11
  
  
  theta_C <- gamma_0 + gamma_1 + gamma_2 + 2*gamma_3 + 2*gamma_4 + 2*gamma_5
  
  theta_D <- 2 * ( -gamma_0*y00 - gamma_1*y10 - gamma_2*y11 - gamma_3*y00 - gamma_3*y10 - gamma_4*y00 -gamma_4*y11 - gamma_5*y10 - gamma_5*y11 +
                     gamma_1*w10*delta + gamma_2*w11*delta + gamma_3*w10*delta + gamma_4*w11*delta + gamma_5*w10*delta + gamma_5*w11*delta )
  
  theta_A <- (x1^2)*theta_C
  
  theta_B <- 2*theta_C*x1*(nu0 + w10*nu1 + w11*nu2 + x2*nu4) + x1*theta_D
  
  mu0_sd <- sqrt( 1 / (sum(theta_A / theta_denom) + 1/k0 ) )
  
  mu0_mean <- (-1) * ( sum(theta_B / theta_denom) - 2*k1/k0 ) / ( 2 * ( sum(theta_A / theta_denom) + 1/k0 ) )
  
  
  mu0_t_draw <- rnorm(1, mean = mu0_mean, sd = mu0_sd)
  
  return(mu0_t_draw)
  
}



# nu4 conditional posterior

nu4_cond_posterior <- function(sigma0, rho, nu0, nu1, nu2, nu3, delta, x1, x2, y00, y10, y11, w10, w11, k0, k1){
  
  sigma11 <- sigma0^2
  
  sigma12 <- rho*sigma0^2
  
  sigma13 <- rho*sigma0^2
  
  sigma23 <- rho*sigma0^2
  
  sigma22 <- sigma0^2
  
  sigma33 <- sigma0^2
  
  
  gamma_0 <- sigma22*sigma33 - sigma23^2
  gamma_1 <- sigma11*sigma33 - sigma13^2
  gamma_2 <- sigma11*sigma22 - sigma12^2
  gamma_3 <- sigma23*sigma13 - sigma12*sigma33
  gamma_4 <- sigma12*sigma23 - sigma13*sigma22
  gamma_5 <- sigma12*sigma13 - sigma11*sigma23
  
  
  theta_denom <- sigma11*sigma22*sigma33 + 2*sigma12*sigma13*sigma23 - sigma12^2 * sigma33 - sigma13^2 * sigma22 - sigma23^2 * sigma11
  
  
  theta_C <- gamma_0 + gamma_1 + gamma_2 + 2*gamma_3 + 2*gamma_4 + 2*gamma_5
  
  theta_D <- 2 * ( -gamma_0*y00 - gamma_1*y10 - gamma_2*y11 - gamma_3*y00 - gamma_3*y10 - gamma_4*y00 -gamma_4*y11 - gamma_5*y10 - gamma_5*y11 +
                     gamma_1*w10*delta + gamma_2*w11*delta + gamma_3*w10*delta + gamma_4*w11*delta + gamma_5*w10*delta + gamma_5*w11*delta )
  
  theta_A <- (x2^2)*theta_C
  
  theta_B <- 2*theta_C*x2*(nu0 + w10*nu1 + w11*nu2 + x1*nu3) + x2*theta_D
  
  mu0_sd <- sqrt( 1 / (sum(theta_A / theta_denom) + 1/k0 ) )
  
  mu0_mean <- (-1) * ( sum(theta_B / theta_denom) - 2*k1/k0 ) / ( 2 * ( sum(theta_A / theta_denom) + 1/k0 ) )
  
  
  mu0_t_draw <- rnorm(1, mean = mu0_mean, sd = mu0_sd)
  
  return(mu0_t_draw)
  
}




# delta conditional posterior

delta_cond_posterior <- function(sigma0, rho, mu0, y00, y10, y11, w10, w11, v0, v1){
  
  sigma11 <- sigma0^2
  
  sigma12 <- rho*sigma0^2
  
  sigma13 <- rho*sigma0^2
  
  sigma23 <- rho*sigma0^2
  
  sigma22 <- sigma0^2
  
  sigma33 <- sigma0^2
  
  
  gamma_0 <- sigma22*sigma33 - sigma23^2
  gamma_1 <- sigma11*sigma33 - sigma13^2
  gamma_2 <- sigma11*sigma22 - sigma12^2
  gamma_3 <- sigma23*sigma13 - sigma12*sigma33
  gamma_4 <- sigma12*sigma23 - sigma13*sigma22
  gamma_5 <- sigma12*sigma13 - sigma11*sigma23
  
  
  theta_denom <- sigma11*sigma22*sigma33 + 2*sigma12*sigma13*sigma23 - sigma12^2 * sigma33 - sigma13^2 * sigma22 - sigma23^2 * sigma11
  
  
  theta_A <- gamma_1*w10 + gamma_2*w11 + 2*gamma_5*w10*w11
  
  theta_B <- 2 * ( gamma_1*mu0*w10 - gamma_1*y10*w10 - gamma_2*y11*w11 + gamma_2*mu0*w11 - gamma_3*y00*w10 + gamma_3*mu0*w10 -
                     gamma_4*y00*w11 + gamma_4*mu0*w11 - gamma_5*y10*w11 - gamma_5*y11*w10 + gamma_5*mu0*w10 + gamma_5*mu0*w11)
  
  
  delta_sd <- sqrt( 1 / (sum(theta_A / theta_denom) + 1/v0 ) )
  
  delta_mean <- (-1) * ( sum(theta_B / theta_denom) - 2*v1/v0 ) / ( 2 * (sum(theta_A / theta_denom) + 1/v0 ) )
  
  
  delta_t_draw <- rnorm(1, mean = delta_mean, sd = delta_sd)
  
  return(delta_t_draw)
  
}





ws_cond_posterior <- function(w10, w11, sigma0, rho, nu0, nu1, nu2, nu3, nu4, delta, x1, x2, y00, y10, y11, gam_p0, gam_p1, gam_p2, gam_q0, gam_q1, gam_q2){
  
  # write out these definitions
  theta_1 <- (y00 - nu0 - w10*nu1 - w11*nu2 - x1*nu3 - x2*nu4)^2
  
  theta_2 <- (y10 - nu0 - w10*nu1 - w11*nu2 - x1*nu3 - x2*nu4 - w10*delta)^2
  
  theta_3 <- (y11 - nu0 - w10*nu1 - w11*nu2 - x1*nu3 - x2*nu4 - w11*delta)^2
  
  theta_4 <- 2*(y00 - nu0 - w10*nu1 - w11*nu2 - x1*nu3 - x2*nu4)*(y10 - nu0 - w10*nu1 - w11*nu2 - x1*nu3 - x2*nu4 - w10*delta)
  
  theta_5 <- 2*(y00 - nu0 - w10*nu1 - w11*nu2 - x1*nu3 - x2*nu4)*(y11 - nu0 - w10*nu1 - w11*nu2 - x1*nu3 - x2*nu4 - w11*delta)
  
  theta_6 <- 2*(y10 - nu0 - w10*nu1 - w11*nu2 - x1*nu3 - x2*nu4 - w10*delta)*(y11 - nu0 - w10*nu1 - w11*nu2 - x1*nu3 - x2*nu4 - w11*delta)
  
  
  sigma11 <- sigma0^2
  
  sigma12 <- rho*sigma0^2
  
  sigma13 <- rho*sigma0^2
  
  sigma23 <- rho*sigma0^2
  
  sigma22 <- sigma0^2
  
  sigma33 <- sigma0^2
  
  p <- pnorm((gam_p0 + x1*gam_p1 + x2*gam_p2))
  
  q <- pnorm((gam_q0 + x1*gam_q1 + x2*gam_q2))
  
  
  theta_big <- theta_1 * (sigma22*sigma33 - sigma23^2) + theta_2 * (sigma11*sigma33 - sigma13^2) + theta_3 * (sigma11*sigma22 - sigma12^2) +
    theta_4 * (sigma23*sigma13 - sigma12*sigma33) + theta_5 * (sigma12*sigma23 - sigma13*sigma22) + theta_6 * (sigma12*sigma13 - sigma11*sigma23)
  
  det_Z <- sigma11*sigma22*sigma33 + 2*sigma12*sigma13*sigma23 - sigma12^2 * sigma33 - sigma13^2 * sigma22 - sigma23^2 * sigma11
  
  post <- (1 / sqrt(det_Z) ) * exp(-0.5*(theta_big / det_Z) ) * ( q^w11 * (1-q)^(1-w11) ) * ( p^w10 * (1-p)^(1-w10) )
  
  return(post)
  
}



###################################################
###################################################

# Start the Gibbs sampler

###################################################
###################################################




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


draws <- 200

n_iterations <- 3000

n_replications <- 100

# ensuring each model fits parameters to the same dataset as one another for each configuration
set.seed(617)
seeds <- round(runif(n_sims*n_replications, 0, 100000))






# Put together the Gibbs sampler and Metropolis Hastings Algorithms

gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2 <- matrix(NA, nrow = n_iterations, ncol = 14)
colnames(gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2) <- c("gam_p0", "gam_p1", "gam_p2", "gam_q0", "gam_q1", "gam_q2", "sigma0", "rho", "nu0", "nu1", "nu2", "nu3", "nu4", "delta")

ITT_vals_model_new_Sim_A_C2 <- matrix(NA, nrow = n_iterations, ncol = n_replications)

ITT_vals_model_new_Sim_A_C2_estimates <- matrix(NA, nrow = n_replications, ncol = n_sims)
ITT_vals_model_new_Sim_A_C2_variance <- matrix(NA, nrow = n_replications, ncol = n_sims)
ITT_vals_model_new_Sim_A_C2_lower <- matrix(NA, nrow = n_replications, ncol = n_sims)
ITT_vals_model_new_Sim_A_C2_upper <- matrix(NA, nrow = n_replications, ncol = n_sims)



for (m in 1:n_sims){
  
  for (h in 1:n_replications){
    
    set.seed(seeds[(n_replications*(m-1) + h)])
    # set.seed(Sys.time())
    
    ndat <- ndat_true[m]
    
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
    
    
    # Unset the seed so that we randomly draw again
    set.seed(Sys.time())
    
    
    
    
    # Generate the random treatment assignment
    
    # need the proportion of treatment individuals before we flip to y11
    prop_trt <- prop_true[m]
    
    # counter for the number of treatment individuals
    count <- 0
    
    for (i in 1:ndat){
      if (design_mat$treatment[i] == 0){
        design_mat$Zvec[i] <- 0
      } else if (count <= (prop_trt * ndat/2)){
        count <- count + 1
        design_mat$Zvec[i] <- 1
      } else if (count > (prop_trt * ndat/2)){
        design_mat$Zvec[i] <- 2
      }
    }
    
    
    # put together the full true data matrix
    
    true_data <- data.frame(potential_outcomes, design_mat)
    
    
    
    
    # now create an observed data matrix with NA for unobserved
    
    observed_outcome_data <- matrix(NA, nrow = ndat, ncol = 6)
    
    for (i in 1:ndat){
      if (true_data$Zvec[i] == 0){
        observed_outcome_data[i,1] <- true_data$y00[i]
      } else if (true_data$Zvec[i] == 1){
        observed_outcome_data[i,2] <- true_data$y10[i]
        observed_outcome_data[i,4] <- true_data$w10[i]
      } else if (true_data$Zvec[i] == 2){
        observed_outcome_data[i,3] <- true_data$y11[i]
        observed_outcome_data[i,5] <- true_data$w11[i]
      }
      observed_outcome_data[i,6] <- true_data$Zvec[i]
    }
    
    colnames(observed_outcome_data) <- c("y00", "y10", "y11", "w10", "w11", "Zvec")
    
    # use these indicators for deriving the mean and covariance of the bivariate normal
    z_assignments <- to.dummy(observed_outcome_data[,6], "Z")
    
    
    # Making sure we get valid probabilities and initial parameter values to start
    w_draw0 <- NULL
    while(is.null(w_draw0)){
      # take the observed outcome matrix to impute missing values
      # fill in the y's with their mean response so that we can impute the w's
      imputed_matrix <- observed_outcome_data
      
      for (j in 1:3){
        test <- observed_outcome_data[,j]
        
        imputed_matrix[,j] <- ifelse(is.na(test), mean(test[which(!(is.na(test)))]), test)
        
      }
      
      
      # need these vectors for placing correct compliance imputations into matrix
      obs_compW10 <- ifelse(is.na(imputed_matrix[,4]), 99, imputed_matrix[,4])
      obs_compW11 <- ifelse(is.na(imputed_matrix[,5]), 99, imputed_matrix[,5])
      
      
      # draw which compliance group each individual falls under
      w_draw0 <- t(rmultinom(ndat, 1, cbind(1/4, 1/4, 1/4, 1/4)))
      
      # Wi(1,0) does not comply => 0,0 or 0,1
      # Wi(1,0) does comply => 1,1 or 1,0 which is wdraw0[,3] and [,4]
      # Wi(1,1) does comply => 0,1 or 1,1 which is wdraw0[,2] and [,3]
      
      # determine bernoulli probabilities - used when zi=(1,0)
      # draw if wi(1,1) complies or not
      w_draw1 <- rbinom(ndat, 1, 1/3)
      
      # determine bernoulli probabilities - used when zi=(1,1)
      # draw if wi(1,0) complies or not
      w_draw2 <- rbinom(ndat, 1, 1/3)
      
      # a function to impute all of the drawn compliances while keeping the observed ones
      imputed_matrix[,4] <- z_assignments[,1] * ( w_draw0[,3] + w_draw0[,4] ) + 
        z_assignments[,2] * ( obs_compW10 ) + 
        z_assignments[,3] * ( w_draw2 )
      
      imputed_matrix[,5] <- z_assignments[,1] * ( w_draw0[,2] + w_draw0[,3] ) + 
        z_assignments[,2] * ( w_draw1 ) + 
        z_assignments[,3] * ( obs_compW11 )
      
      
      
      # need this to indicate which column the observed value is in
      obs_respY <- sapply(1:ndat, function(x) which(!(is.na(observed_outcome_data[x,1:3]))) )
      names(obs_respY) <- NULL
      
      
      # need these two arrays for the multivariate normal covariance for each person
      lrg_sigma11 <- array(NA, dim = c(2,2,ndat))
      
      lrg_sigma_mult <- array(NA, dim = c(2,2,ndat))
      
      
      optim_output <- optim(par = c(1, 0.5, 0, 0, 0, 0, 0, 1),    
                            fn = log_lik_fnct_full,
                            x1 = x1,
                            x2 = x2,
                            y00 = imputed_matrix[,1],
                            y10 = imputed_matrix[,2],
                            y11 = imputed_matrix[,3],
                            w10 = imputed_matrix[,4],
                            w11 = imputed_matrix[,5],
                            control=list(fnscale=-1),
                            hessian=FALSE)
      
      mu_norm <- optim_output$par
      
      
      sigma0_t <- mu_norm[1]
      rho_t <- mu_norm[2]
      nu0_t <- mu_norm[3]
      nu1_t <- mu_norm[4]
      nu2_t <- mu_norm[5]
      nu3_t <- mu_norm[6]
      nu4_t <- mu_norm[7]
      delta_t <- mu_norm[8]
      
      
      optim_output <- optim(par = c(0.5, 0, 0, 0.5, 0, 0),    
                            fn = log_lik_fnct_WS,
                            x1 = x1,
                            x2 = x2,
                            w10 = imputed_matrix[,4],
                            w11 = imputed_matrix[,5],
                            control=list(fnscale=-1),
                            hessian=FALSE)
      
      mu_norm <- optim_output$par
      
      
      gam_p0_t <- mu_norm[1]
      gam_p1_t <- mu_norm[2]
      gam_p2_t <- mu_norm[3]
      gam_q0_t <- mu_norm[4]
      gam_q1_t <- mu_norm[5]
      gam_q2_t <- mu_norm[6]
      
      
      # Make sure we have valid initial values for the parameters
      if (is.na(sigma0_t) | sigma0_t <= 0){
        sigma0_t <- runif(1, 0.5, 5)
      } 
      if (is.na(rho_t) | rho_t <= 0 | rho_t >= 1){
        rho_t <- runif(1, 0.1, 0.9)
      }
      
      
      nu0_t <- ifelse(is.na(nu0_t), runif(1, -5, 5), nu0_t)
      nu1_t <- ifelse(is.na(nu1_t), runif(1, -5, 5), nu1_t)
      nu2_t <- ifelse(is.na(nu2_t), runif(1, -5, 5), nu2_t)
      nu3_t <- ifelse(is.na(nu2_t), runif(1, -5, 5), nu3_t)
      nu4_t <- ifelse(is.na(nu2_t), runif(1, -5, 5), nu4_t)
      delta_t <- ifelse(is.na(delta_t), runif(1, -5, 5), delta_t)
      
      
      mu0_t <- nu0_t + imputed_matrix[,4]*nu1_t + imputed_matrix[,5]*nu2_t + x1*nu3_t + x2*nu4_t
      
      
      gam_p0_t <- ifelse(is.na(gam_p0_t), runif(1, -1, 1), gam_p0_t)
      gam_p1_t <- ifelse(is.na(gam_p1_t), runif(1, -1, 1), gam_p1_t)
      gam_p2_t <- ifelse(is.na(gam_p2_t), runif(1, -1, 1), gam_p2_t)
      gam_q0_t <- ifelse(is.na(gam_q0_t), runif(1, -1, 1), gam_q0_t)
      gam_q1_t <- ifelse(is.na(gam_q1_t), runif(1, -1, 1), gam_q1_t)
      gam_q2_t <- ifelse(is.na(gam_q2_t), runif(1, -1, 1), gam_q2_t)
      
      # determine the likelihood values with each combination of compliance values
      p00 <- ws_cond_posterior(0, 0, sigma0_t, rho_t, nu0_t, nu1_t, nu2_t, nu3_t, nu4_t, delta_t, x1, x2,
                               imputed_matrix[,1], imputed_matrix[,2], imputed_matrix[,3], gam_p0_t, gam_p1_t, gam_p2_t, gam_q0_t, gam_q1_t, gam_q2_t)
      
      p01 <- ws_cond_posterior(0, 1, sigma0_t, rho_t, nu0_t, nu1_t, nu2_t, nu3_t, nu4_t, delta_t, x1, x2,
                               imputed_matrix[,1], imputed_matrix[,2], imputed_matrix[,3], gam_p0_t, gam_p1_t, gam_p2_t, gam_q0_t, gam_q1_t, gam_q2_t)
      
      p11 <- ws_cond_posterior(1, 1, sigma0_t, rho_t, nu0_t, nu1_t, nu2_t, nu3_t, nu4_t, delta_t, x1, x2,
                               imputed_matrix[,1], imputed_matrix[,2], imputed_matrix[,3], gam_p0_t, gam_p1_t, gam_p2_t, gam_q0_t, gam_q1_t, gam_q2_t)
      
      p10 <- ws_cond_posterior(1, 0, sigma0_t, rho_t, nu0_t, nu1_t, nu2_t, nu3_t, nu4_t, delta_t, x1, x2,
                               imputed_matrix[,1], imputed_matrix[,2], imputed_matrix[,3], gam_p0_t, gam_p1_t, gam_p2_t, gam_q0_t, gam_q1_t, gam_q2_t)
      
      # determine the multinomial probabilities - used when zi = (0,0)
      p1 <- p00 / (p00 + p01 + p11 + p10)
      
      p2 <- p01 / (p00 + p01 + p11 + p10)
      
      p3 <- p11 / (p00 + p01 + p11 + p10)
      
      p4 <- p10 / (p00 + p01 + p11 + p10)
      # draw which compliance group each individual falls under
      w_draw0 <- tryCatch(t(sapply(1:ndat, function(x) rmultinom(1, 1, cbind(p1, p2, p3, p4)[x,]) )), error = function(e) NULL)
      
      
      
      
    }
    
    
    # Initial values for u1_t
    u1_t <- rep(0, ndat)
    
    # define mean vector for latent variable mean
    mu_U1 <- xmat_U1 %*% c(gam_p0_t, gam_p1_t, gam_p2_t)
    
    # draw new latent variable values conditional on the compliance
    
    if ((ndat-sum(imputed_matrix[,4])) > 0){
      u1_t[imputed_matrix[,4] == 0] <- rtruncnorm( (ndat-sum(imputed_matrix[,4])) , mean = mu_U1[imputed_matrix[,4] == 0], sd = 1, a = -Inf, b = 0)
    }
    
    if ((sum(imputed_matrix[,4])) > 0){
      u1_t[imputed_matrix[,4] == 1] <- rtruncnorm( (sum(imputed_matrix[,4])) , mean = mu_U1[imputed_matrix[,4] == 1], sd = 1, a = 0, b = Inf)
    }
    
    
    
    # Initial values for u2_t
    u2_t <- rep(0, ndat)
    
    # define mean vector for latent variable mean
    mu_U2 <- xmat_U2 %*% c(gam_q0_t, gam_q1_t, gam_q2_t)
    
    
    if ((ndat-sum(imputed_matrix[,5])) > 0){
      u2_t[imputed_matrix[,5] == 0] <- rtruncnorm( (ndat-sum(imputed_matrix[,5])) , mean = mu_U2[imputed_matrix[,5] == 0], sd = 1, a = -Inf, b = 0)
    }
    
    if ((sum(imputed_matrix[,5])) > 0){
      u2_t[imputed_matrix[,5] == 1] <- rtruncnorm( (sum(imputed_matrix[,5])) , mean = mu_U2[imputed_matrix[,5] == 1], sd = 1, a = 0, b = Inf)
    }
    
    
    
    
    
    seeds_loop <- round(runif(n_iterations, 0, 1000000))
    
    for (n in 1:n_iterations){
      
      set.seed(seeds_loop[n])
      
      
      # determine the likelihood values with each combination of compliance values
      p00 <- ws_cond_posterior(0, 0, sigma0_t, rho_t, nu0_t, nu1_t, nu2_t, nu3_t, nu4_t, delta_t, x1, x2,
                               imputed_matrix[,1], imputed_matrix[,2], imputed_matrix[,3], gam_p0_t, gam_p1_t, gam_p2_t, gam_q0_t, gam_q1_t, gam_q2_t)
      
      p01 <- ws_cond_posterior(0, 1, sigma0_t, rho_t, nu0_t, nu1_t, nu2_t, nu3_t, nu4_t, delta_t, x1, x2,
                               imputed_matrix[,1], imputed_matrix[,2], imputed_matrix[,3], gam_p0_t, gam_p1_t, gam_p2_t, gam_q0_t, gam_q1_t, gam_q2_t)
      
      p11 <- ws_cond_posterior(1, 1, sigma0_t, rho_t, nu0_t, nu1_t, nu2_t, nu3_t, nu4_t, delta_t, x1, x2,
                               imputed_matrix[,1], imputed_matrix[,2], imputed_matrix[,3], gam_p0_t, gam_p1_t, gam_p2_t, gam_q0_t, gam_q1_t, gam_q2_t)
      
      p10 <- ws_cond_posterior(1, 0, sigma0_t, rho_t, nu0_t, nu1_t, nu2_t, nu3_t, nu4_t, delta_t, x1, x2,
                               imputed_matrix[,1], imputed_matrix[,2], imputed_matrix[,3], gam_p0_t, gam_p1_t, gam_p2_t, gam_q0_t, gam_q1_t, gam_q2_t)
      
      # determine the multinomial probabilities - used when zi = (0,0)
      p1 <- p00 / (p00 + p01 + p11 + p10)
      
      p2 <- p01 / (p00 + p01 + p11 + p10)
      
      p3 <- p11 / (p00 + p01 + p11 + p10)
      
      p4 <- p10 / (p00 + p01 + p11 + p10)
      # draw which compliance group each individual falls under
      w_draw0 <- t(sapply(1:ndat, function(x) rmultinom(1, 1, cbind(p1, p2, p3, p4)[x,]) ))
      
      # the ifelse serves in place of more complex if else statement required
      # na's only appear sometimes in settings where those compliance statuses are impossible for individuals whose observed status makes it impossible
      # Otherwise na's do not appear
      
      # When Zi=(1,0) and complier -> Prob Wi(1,1) = 1 is:
      p1 <- p11 / (p10 + p11)
      p1 <- ifelse(is.na(p1), 0.5, p1)
      # draw if wi(1,0) complies
      w_drawZ11_comp <- sapply(1:ndat, function(x) rbinom(1, 1, p1[x]))
      
      # When Zi=(1,0) and noncomplier -> Prob Wi(1,1) = 1 is:
      p1 <- p01 / (p00 + p01)
      p1 <- ifelse(is.na(p1), 0.5, p1)
      # draw if wi(1,0) complies
      w_drawZ11_noncomp <- sapply(1:ndat, function(x) rbinom(1, 1, p1[x]))
      
      
      # When Zi=(1,1) and complier -> Prob Wi(1,0) = 1 is:
      p1 <- p11 / (p01 + p11)
      p1 <- ifelse(is.na(p1), 0.5, p1)
      # draw if wi(1,0) complies
      w_drawZ10_comp <- sapply(1:ndat, function(x) rbinom(1, 1, p1[x]))
      
      # When Zi=(1,1) and noncomplier -> Prob Wi(1,0) = 1 is:
      p1 <- p10 / (p00 + p10)
      p1 <- ifelse(is.na(p1), 0.5, p1)
      w_drawZ10_noncomp <- sapply(1:ndat, function(x) rbinom(1, 1, p1[x]))
      
      
      # a function to impute all of the drawn compliances while keeping the observed ones
      imputed_matrix[,4] <- z_assignments[,1] * ( w_draw0[,3] + w_draw0[,4] ) + 
        z_assignments[,2] * ( obs_compW10 ) + 
        z_assignments[,3] * ( w_drawZ10_comp*obs_compW11 + w_drawZ10_noncomp*(1-obs_compW11) )
      
      imputed_matrix[,5] <- z_assignments[,1] * ( w_draw0[,2] + w_draw0[,3] ) + 
        z_assignments[,2] * ( w_drawZ11_comp*obs_compW10 + w_drawZ11_noncomp*(1-obs_compW10) ) + 
        z_assignments[,3] * ( obs_compW11 )
      
      
      # Need to update mu0 here
      
      mu0_t <- nu0_t + imputed_matrix[,4]*nu1_t + imputed_matrix[,5]*nu2_t + x1*nu3_t + x2*nu4_t
      
      # Now we need to impute the missing potential response outcomes
      # first we need to define the Ymis conditional bivariate distribution
      # then when we draw missing values we'll need to place them correctly in the imputed matrix
      # here we define the conditional bivariate mean as derived in my notes for the i^th individual
      
      
      
      bvn_mu <- cbind((mu0_t + delta_t*z_assignments[,1]*imputed_matrix[,4]),
                      (mu0_t + delta_t*(z_assignments[,1]*imputed_matrix[,5] + z_assignments[,2]*imputed_matrix[,5] + z_assignments[,3]*imputed_matrix[,4]))) +
        ( cbind(rep((rho_t*(sigma0_t^2)), ndat), rep((rho_t*(sigma0_t^2)), ndat)) /
            (sigma0_t^2) * 
            (imputed_matrix[cbind(1:ndat,obs_respY)] - mu0_t - delta_t*(z_assignments[,2]*imputed_matrix[,4] + z_assignments[,3]*imputed_matrix[,5]))   )
      
      
      
      
      # define the 2x2 array of matrices to be used to compute the covariance matrix of our conditional bivariate normal
      
      
      lrg_sigma11[1,1,] <- rep((sigma0_t^2), ndat)
      lrg_sigma11[2,1,] <- rep((rho_t*(sigma0_t^2)), ndat)
      lrg_sigma11[1,2,] <- rep((rho_t*(sigma0_t^2)), ndat)
      lrg_sigma11[2,2,] <- rep((sigma0_t^2), ndat)
      
      
      
      lrg_sigma_mult[1,1,] <- rep(((rho_t*sigma0_t)^2), ndat)
      lrg_sigma_mult[2,1,] <- rep(((rho_t*sigma0_t)^2), ndat)
      lrg_sigma_mult[1,2,] <- rep(((rho_t*sigma0_t)^2), ndat)
      lrg_sigma_mult[2,2,] <- rep(((rho_t*sigma0_t)^2), ndat)
      
      
      # subtract the two matrices above to get our covariance
      bvn_sigma <- lrg_sigma11 - lrg_sigma_mult
      
      # draw the Ymis vector for the i^th individual
      ymis_i <- t(sapply(1:ndat, function(x) rmvnorm(1, bvn_mu[x,], bvn_sigma[ , , x]) ))
      
      
      # place the missing values into the correct cells
      imputed_matrix[,1] <- z_assignments[,1] * imputed_matrix[,1] + (z_assignments[,2] + z_assignments[,3]) * ymis_i[,1]
      
      imputed_matrix[,2] <- z_assignments[,1] * ymis_i[,1] + z_assignments[,2] * imputed_matrix[,2] + z_assignments[,3] * ymis_i[,2]
      
      imputed_matrix[,3] <- (z_assignments[,1] + z_assignments[,2]) * ymis_i[,2] + z_assignments[,3] * imputed_matrix[,3]
      
      
      # now save this dataframe for later use
      # gibbs_imputed_data <- rbind(gibbs_imputed_data, imputed_matrix)
      
      
      # next steps require drawing from full conditional posterior distributions
      
      # define mean vector for latent variable mean
      mu_U1 <- xmat_U1 %*% c(gam_p0_t, gam_p1_t, gam_p2_t)
      
      # draw new latent variable values conditional on the compliance
      if ((ndat-sum(imputed_matrix[,4])) > 0){
        u1_t[imputed_matrix[,4] == 0] <- rtruncnorm( (ndat-sum(imputed_matrix[,4])) , mean = mu_U1[imputed_matrix[,4] == 0], sd = 1, a = -Inf, b = 0)
      }
      
      if ((sum(imputed_matrix[,4])) > 0){
        u1_t[imputed_matrix[,4] == 1] <- rtruncnorm( (sum(imputed_matrix[,4])) , mean = mu_U1[imputed_matrix[,4] == 1], sd = 1, a = 0, b = Inf)
      }
      
      # define the new posterior mean and covariance for the covariates
      gam_p_mean <- solve(t(xmat_U1) %*% xmat_U1) %*% t(xmat_U1) %*% u1_t
      
      gam_p_cov <- solve(t(xmat_U1) %*% xmat_U1) 
      
      gam_p_t <- c(rmvnorm(1, gam_p_mean, gam_p_cov))
      
      
      # save the new parameters
      gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,1] <- gam_p_t[1]
      # Then assign it to the p_t variable use throughout
      gam_p0_t <- gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,1]
      
      
      # save the new parameters
      gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,2] <- gam_p_t[2]
      # Then assign it to the p_t variable use throughout
      gam_p1_t <- gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,2]
      
      
      # save the new parameters
      gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,3] <- gam_p_t[3]
      # Then assign it to the p_t variable use throughout
      gam_p2_t <- gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,3]
      
      
      # Need to define the covariates at the beginning and give them a matrix form for this
      # Then need to define the u1_t vector where we'll be imputing values throughout
      
      # define mean vector for latent variable mean
      mu_U2 <- xmat_U2 %*% c(gam_q0_t, gam_q1_t, gam_q2_t)
      
      
      if ((ndat-sum(imputed_matrix[,5])) > 0){
        u2_t[imputed_matrix[,5] == 0] <- rtruncnorm( (ndat-sum(imputed_matrix[,5])) , mean = mu_U2[imputed_matrix[,5] == 0], sd = 1, a = -Inf, b = 0)
      }
      
      if ((sum(imputed_matrix[,5])) > 0){
        u2_t[imputed_matrix[,5] == 1] <- rtruncnorm( (sum(imputed_matrix[,5])) , mean = mu_U2[imputed_matrix[,5] == 1], sd = 1, a = 0, b = Inf)
      }
      
      
      
      # define the new posterior mean and covariance for the covariates
      gam_q_mean <- solve(t(xmat_U2) %*% xmat_U2) %*% t(xmat_U2) %*% u2_t
      
      gam_q_cov <- solve(t(xmat_U2) %*% xmat_U2) 
      
      gam_q_t <- c(rmvnorm(1, gam_q_mean, gam_q_cov))
      
      
      # save the new parameters
      gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,4] <- gam_q_t[1]
      # Then assign it to the p_t variable use throughout
      gam_q0_t <- gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,4]
      
      
      # save the new parameters
      gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,5] <- gam_q_t[2]
      # Then assign it to the p_t variable use throughout
      gam_q1_t <- gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,5]
      
      
      # save the new parameters
      gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,6] <- gam_q_t[3]
      # Then assign it to the p_t variable use throughout
      gam_q2_t <- gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,6]
      
      
      # # Draw the next p probability and save it
      # gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,1] <- 0.7
      # # Then assign it to the p_t variable use throughout
      # p_t <- gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,1]
      # 
      # # Draw the next q probability and save it
      # gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,2] <- 0.6
      # # Then assign it to the q_t variable use throughout
      # q_t <- gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,2]
      
      
      
      
      
      #######################
      # MH For SIGMA0
      ######################
      
      draw_d <- sigma0_t^2
      sig_norm <- 0.1
      
      # here we run an MH for (sigma0)^2 and then transform back to sigma0 at the end to keep the final draw
      for (d in 1:draws){
        
        proposal <- rtruncnorm(1, a = 0.0001, mean = draw_d, sd = 2.4*sig_norm)
        
        r_top <- sigma0_cond_posterior(proposal, rho_t, mu0_t, delta_t, imputed_matrix[,1], imputed_matrix[,2], imputed_matrix[,3], imputed_matrix[,4], imputed_matrix[,5]) - log(dtruncnorm(proposal, a = 0.0001, mean = draw_d, sd = 2.4*sig_norm))
        
        r_bottom <- sigma0_cond_posterior(draw_d, rho_t, mu0_t, delta_t, imputed_matrix[,1], imputed_matrix[,2], imputed_matrix[,3], imputed_matrix[,4], imputed_matrix[,5]) - log(dtruncnorm(draw_d, a = 0.0001, mean = proposal, sd = 2.4*sig_norm))
        
        
        r <- r_top - r_bottom
        
        # if true accept the new proposal
        if (r > log(runif(1))){
          draw_d <- proposal
        } 
      }
      
      # keep the final draw from the MH for sigma0
      # must transform back from log(sigma0) to sigma0
      gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,7] <- sqrt(draw_d)
      
      # update the sigma0 parameter to reflect this
      sigma0_t <- gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,7]
      
      
      
      
      # sigma0_mat <- data.frame(1:draws, sigma0_chain)
      # names(sigma0_mat) <- c("draws", "sigma0_vals")
      # ggplot(sigma0_mat, aes(draws, sigma0_vals)) + geom_line()
      
      
      
      
      
      #######################
      # MH For RHO
      ######################
      
      
      
      draw_d <- rho_t
      
      sig_norm <- 0.025
      
      # here we run an MH for (rho)^2 and then transform back to rho at the end to keep the final draw
      for (d in 1:draws){
        
        proposal <- rtruncnorm(1, a = 0.0001, mean = draw_d, sd = 2.4*sig_norm)
        
        # we can't run any iterations when proposal = 1 because this gives x/0 from the determinant
        if (proposal < 0.9999){
          r_top <- rho_cond_posterior(proposal, sigma0_t, mu0_t, delta_t, imputed_matrix[,1], imputed_matrix[,2], imputed_matrix[,3], imputed_matrix[,4], imputed_matrix[,5]) - log(dtruncnorm(proposal, a = 0.0001, mean = draw_d, sd = 2.4*sig_norm))
          
          r_bottom <- rho_cond_posterior(draw_d, sigma0_t, mu0_t, delta_t, imputed_matrix[,1], imputed_matrix[,2], imputed_matrix[,3], imputed_matrix[,4], imputed_matrix[,5]) - log(dtruncnorm(draw_d, a = 0.0001, mean = proposal, sd = 2.4*sig_norm))
          
          
          r <- r_top - r_bottom
          
          # if true accept the new proposal
          if (r > log(runif(1))){
            draw_d <- proposal
          } 
        }
        
      }
      
      
      # if (draw_d >= 1){
      #   draw_d <- 0.95
      # }
      
      # keep the final draw from the MH for rho
      # must transform back from log(rho) to rho
      gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,8] <- draw_d
      
      # update the rho parameter to reflect this
      rho_t <- gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,8]
      
      
      
      
      # rho_mat <- data.frame(1:draws, rho_chain)
      # names(rho_mat) <- c("draws", "rho_vals")
      # ggplot(rho_mat, aes(draws, rho_vals)) + geom_line()
      
      
      # nu 0 draw
      nu0_t_draw <- nu0_cond_posterior(sigma0_t, rho_t, nu1_t, nu2_t, nu3_t, nu4_t, delta_t, x1, x2, imputed_matrix[,1], imputed_matrix[,2], imputed_matrix[,3], imputed_matrix[,4], imputed_matrix[,5], k0_mu0, k1_mu0)
      
      # Save the drawn mu0 
      gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,9] <- nu0_t_draw
      # Then update this value
      nu0_t <- gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,9]
      
      
      # nu1 draw
      nu1_t_draw <- nu1_cond_posterior(sigma0_t, rho_t, nu0_t, nu2_t, nu3_t, nu4_t, delta_t, x1, x2, imputed_matrix[,1], imputed_matrix[,2], imputed_matrix[,3], imputed_matrix[,4], imputed_matrix[,5], k0_mu0, k1_mu0)
      
      # Save the drawn mu0 
      gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,10] <- nu1_t_draw
      # Then update this value
      nu1_t <- gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,10]
      
      
      # nu 0 draw
      nu2_t_draw <- nu2_cond_posterior(sigma0_t, rho_t, nu0_t, nu1_t, nu3_t, nu4_t, delta_t, x1, x2, imputed_matrix[,1], imputed_matrix[,2], imputed_matrix[,3], imputed_matrix[,4], imputed_matrix[,5], k0_mu0, k1_mu0)
      
      # Save the drawn mu0 
      gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,11] <- nu2_t_draw
      # Then update this value
      nu2_t <- gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,11]
      
      
      
      # nu 0 draw
      nu3_t_draw <- nu3_cond_posterior(sigma0_t, rho_t, nu0_t, nu1_t, nu2_t, nu4_t, delta_t, x1, x2, imputed_matrix[,1], imputed_matrix[,2], imputed_matrix[,3], imputed_matrix[,4], imputed_matrix[,5], k0_mu0, k1_mu0)
      
      # Save the drawn mu0 
      gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,12] <- nu3_t_draw
      # Then update this value
      nu3_t <- gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,12]
      
      
      
      # nu 0 draw
      nu4_t_draw <- nu4_cond_posterior(sigma0_t, rho_t, nu0_t, nu1_t, nu2_t, nu3_t, delta_t, x1, x2, imputed_matrix[,1], imputed_matrix[,2], imputed_matrix[,3], imputed_matrix[,4], imputed_matrix[,5], k0_mu0, k1_mu0)
      
      # Save the drawn mu0 
      gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,13] <- nu4_t_draw
      # Then update this value
      nu4_t <- gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,13]
      
      
      # update the mu0_t draw
      
      mu0_t <- nu0_t + imputed_matrix[,4]*nu1_t + imputed_matrix[,5]*nu2_t + x1*nu3_t + x2*nu4_t
      
      
      # delta draw
      delta_t_draw <- delta_cond_posterior(sigma0_t, rho_t, mu0_t, imputed_matrix[,1], imputed_matrix[,2], imputed_matrix[,3], imputed_matrix[,4], imputed_matrix[,5], k0_mu0, k1_mu0)
      
      # Save the drawn delta 
      gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,14] <- delta_t_draw
      # Then update this value
      delta_t <- gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2[n,14]
      
      
      
      
      
      
      
      # w10_t <- rbinom(ndat, 1, p_t)
      # 
      # w11_t <- (rbinom(ndat, 1, q_t))^(1-w10_t)
      # 
      # 
      # mu_y00_t <- mu0_t
      # mu_y10_t <- mu0_t + w10_t*delta_t
      # mu_y11_t <- mu0_t + w11_t*delta_t
      # 
      # mu_vec_t <- cbind(mu_y00_t, mu_y10_t, mu_y11_t)
      # 
      # 
      # 
      # 
      # sigma10_t <- w10_t*sigma1_t^2 + (1-w10_t)*sigma0_t^2
      # sigma11_t <- w11_t*sigma1_t^2 + (1-w11_t)*sigma0_t^2
      # 
      # sigma_cov_t <- array(NA, dim = c(3,3,ndat))
      # 
      # sigma_cov_t[1,1, ] <- sigma0_t^2
      # sigma_cov_t[2,2, ] <- sigma10_t
      # sigma_cov_t[3,3, ] <- sigma11_t
      # 
      # sigma_cov_t[2,1, ] <- rho_t*sigma0_t*sqrt(sigma10_t)
      # sigma_cov_t[1,2, ] <- rho_t*sigma0_t*sqrt(sigma10_t)
      # 
      # sigma_cov_t[1,3, ] <- rho_t*sigma0_t*sqrt(sigma11_t)
      # sigma_cov_t[3,1, ] <- rho_t*sigma0_t*sqrt(sigma11_t)
      # 
      # 
      # sigma_cov_t[2,3, ] <- rho_t*sqrt(sigma10_t)*sqrt(sigma11_t)
      # sigma_cov_t[3,2, ] <- rho_t*sqrt(sigma10_t)*sqrt(sigma11_t)
      # 
      # 
      # y_t <- t(sapply(1:ndat, function(x) rmvnorm(1, mu_vec_t[x,], sigma_cov_t[ , , x]) ))
      # 
      
      
      
      q_t <- mean(pnorm((gam_q0_t + x1*gam_q1_t + x2*gam_q2_t)))
      
      ITT_vals_model_new_Sim_A_C2[n,h] <- q_t * delta_t
      
      
    }
    
    
  }
  
  
  ITT_vals_model_new_Sim_A_C2_estimates[,m] <- apply(ITT_vals_model_new_Sim_A_C2, 2, mean)
  ITT_vals_model_new_Sim_A_C2_variance[,m] <- apply(ITT_vals_model_new_Sim_A_C2, 2, var)
  ITT_vals_model_new_Sim_A_C2_lower[,m] <- apply(ITT_vals_model_new_Sim_A_C2, 2,quantile,c(0.025))
  ITT_vals_model_new_Sim_A_C2_upper[,m] <- apply(ITT_vals_model_new_Sim_A_C2, 2,quantile,c(0.975))
  
  
  save(ITT_vals_model_new_Sim_A_C2_estimates, file = "ITT_vals_model_new_Sim_A_C2_estimates.Rdata")
  
  save(ITT_vals_model_new_Sim_A_C2_variance, file = "ITT_vals_model_new_Sim_A_C2_variance.Rdata")
  
  save(ITT_vals_model_new_Sim_A_C2_lower, file = "ITT_vals_model_new_Sim_A_C2_lower.Rdata")
  
  save(ITT_vals_model_new_Sim_A_C2_upper, file = "ITT_vals_model_new_Sim_A_C2_upper.Rdata")
  
  
  
}


# save(gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2, file = "gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2.Rdata")
# 
# 
# save(ITT_true_vals_model_new_Sim_A_C2, file = "ITT_true_vals_model_new_Sim_A_C2.Rdata")






# gibbs_parameter_posteriors <- gibbs_parameter_posteriors_vectorized_model_new_Sim_A_C2
# 
# itt_vals_model <- ITT_vals_model_new_Sim_A_C2
# 
# # n_iterations <- 10000
# gibbs_post_mat <- data.frame(1:n_iterations, gibbs_parameter_posteriors, itt_vals_model)
# 
# names(gibbs_post_mat) <- c("draws", "p0", "p1", "p2", "q0", "q1", "q2", "sigma0", "rho", "nu0", "nu1", "nu2", "nu3", "nu4", "delta", "SP_ITT")
# 
# ggplot(gibbs_post_mat, aes(draws, p0)) + geom_line()
# 
# ggplot(gibbs_post_mat, aes(draws, p1)) + geom_line()
# 
# ggplot(gibbs_post_mat, aes(draws, p2)) + geom_line()
# 
# ggplot(gibbs_post_mat, aes(draws, q0)) + geom_line()
# 
# ggplot(gibbs_post_mat, aes(draws, q1)) + geom_line()
# 
# ggplot(gibbs_post_mat, aes(draws, q2)) + geom_line()
# 
# ggplot(gibbs_post_mat, aes(draws, sigma0)) + geom_line()
# 
# ggplot(gibbs_post_mat, aes(draws, rho)) + geom_line()
# 
# ggplot(gibbs_post_mat, aes(draws, nu0)) + geom_line()
# 
# ggplot(gibbs_post_mat, aes(draws, nu1)) + geom_line()
# 
# ggplot(gibbs_post_mat, aes(draws, nu2)) + geom_line()
# 
# ggplot(gibbs_post_mat, aes(draws, nu3)) + geom_line()
# 
# ggplot(gibbs_post_mat, aes(draws, nu4)) + geom_line()
# 
# ggplot(gibbs_post_mat, aes(draws, delta)) + geom_line()
# 
# ggplot(gibbs_post_mat, aes(draws, SP_ITT)) + geom_line()


# # var(ITT_vals_model_new_Sim_A_C2)
# 
# 
# 
# 
# draws <- 10000
# ptest <- rbeta(draws, (sum(potential_outcomes$w10) + alpha_p), (beta_p + ndat - sum(potential_outcomes$w10)))
# 
# ptest_dat <- data.frame(1:draws, ptest)
# ggplot(ptest_dat, aes(X1.draws, ptest)) + geom_line()
# 
# 
# qtest <- rbeta(draws, (alpha_q + sum(potential_outcomes$w11) - sum(potential_outcomes$w10*potential_outcomes$w11)),
#                (beta_q + ndat - sum(potential_outcomes$w11) - sum(potential_outcomes$w10) + sum(potential_outcomes$w10*potential_outcomes$w11)))
# 
# qtest_dat <- data.frame(1:draws, qtest)
# ggplot(ptest_dat, aes(X1.draws, qtest)) + geom_line()
# 
# 
# 
# 
# ############################################################################################################################
# ############################################################################################################################
# ############################################################################################################################
# ############################################################################################################################
# 
# 
# 
# 
# 
# # nu0 conditional posterior
# 
# nu0_cond_posterior <- function(sigma0, rho, nu1, nu2, delta, x1, x2, y00, y10, y11, w10, w11, k0, k1){
#   
#   sigma11 <- sigma0^2
#   
#   sigma12 <- rho*sigma0^2
#   
#   sigma13 <- rho*sigma0^2
#   
#   sigma23 <- rho*sigma0^2
#   
#   sigma22 <- sigma0^2
#   
#   sigma33 <- sigma0^2
#   
#   
#   gamma_0 <- sigma22*sigma33 - sigma23^2
#   gamma_1 <- sigma11*sigma33 - sigma13^2
#   gamma_2 <- sigma11*sigma22 - sigma12^2
#   gamma_3 <- sigma23*sigma13 - sigma12*sigma33
#   gamma_4 <- sigma12*sigma23 - sigma13*sigma22
#   gamma_5 <- sigma12*sigma13 - sigma11*sigma23
#   
#   
#   theta_denom <- sigma11*sigma22*sigma33 + 2*sigma12*sigma13*sigma23 - sigma12^2 * sigma33 - sigma13^2 * sigma22 - sigma23^2 * sigma11
#   
#   
#   theta_C <- gamma_0 + gamma_1 + gamma_2 + 2*gamma_3 + 2*gamma_4 + 2*gamma_5
#   
#   theta_D <- 2 * ( -gamma_0*y00 - gamma_1*y10 - gamma_2*y11 - gamma_3*y00 - gamma_3*y10 - gamma_4*y00 -gamma_4*y11 - gamma_5*y10 - gamma_5*y11 +
#                      gamma_1*w10*delta + gamma_2*w11*delta + gamma_3*w10*delta + gamma_4*w11*delta + gamma_5*w10*delta + gamma_5*w11*delta )
#   
#   theta_A <- theta_C
#   
#   theta_B <- 2*theta_C*(w10*nu1 + x2*w11*nu2) + theta_D
#   
#   mu0_sd <- sqrt( 1 / (ndat*(theta_A / theta_denom) + 1/k0 ) )
#   
#   mu0_mean <- (-1) * ( sum(theta_B / theta_denom) - 2*k1/k0 ) / ( 2 * ( ndat*(theta_A / theta_denom) + 1/k0 ) )
#   
#   
#   mu0_t_draw <- rnorm(10000, mean = mu0_mean, sd = mu0_sd)
#   
#   return(mu0_t_draw)
#   
# }
# 
# 
# 
# 
# # nu1 conditional posterior
# 
# nu1_cond_posterior <- function(sigma0, rho, nu0, nu2, delta, x1, x2, y00, y10, y11, w10, w11, k0, k1){
#   
#   sigma11 <- sigma0^2
#   
#   sigma12 <- rho*sigma0^2
#   
#   sigma13 <- rho*sigma0^2
#   
#   sigma23 <- rho*sigma0^2
#   
#   sigma22 <- sigma0^2
#   
#   sigma33 <- sigma0^2
#   
#   
#   gamma_0 <- sigma22*sigma33 - sigma23^2
#   gamma_1 <- sigma11*sigma33 - sigma13^2
#   gamma_2 <- sigma11*sigma22 - sigma12^2
#   gamma_3 <- sigma23*sigma13 - sigma12*sigma33
#   gamma_4 <- sigma12*sigma23 - sigma13*sigma22
#   gamma_5 <- sigma12*sigma13 - sigma11*sigma23
#   
#   
#   theta_denom <- sigma11*sigma22*sigma33 + 2*sigma12*sigma13*sigma23 - sigma12^2 * sigma33 - sigma13^2 * sigma22 - sigma23^2 * sigma11
#   
#   
#   theta_C <- gamma_0 + gamma_1 + gamma_2 + 2*gamma_3 + 2*gamma_4 + 2*gamma_5
#   
#   theta_D <- 2 * ( -gamma_0*y00 - gamma_1*y10 - gamma_2*y11 - gamma_3*y00 - gamma_3*y10 - gamma_4*y00 -gamma_4*y11 - gamma_5*y10 - gamma_5*y11 +
#                      gamma_1*w10*delta + gamma_2*w11*delta + gamma_3*w10*delta + gamma_4*w11*delta + gamma_5*w10*delta + gamma_5*w11*delta )
#   
#   theta_A <- (x1^2)*w10*theta_C
#   
#   theta_B <- 2*theta_C*x1*w10*(nu0 + x2*w11*nu2) + x1*w10*theta_D
#   
#   mu0_sd <- sqrt( 1 / (sum(theta_A / theta_denom) + 1/k0 ) )
#   
#   mu0_mean <- (-1) * ( sum(theta_B / theta_denom) - 2*k1/k0 ) / ( 2 * ( sum(theta_A / theta_denom) + 1/k0 ) )
#   
#   
#   mu0_t_draw <- rnorm(10000, mean = mu0_mean, sd = mu0_sd)
#   
#   return(mu0_t_draw)
#   
# }
# 
# 
# 
# 
# # nu2 conditional posterior
# 
# nu2_cond_posterior <- function(sigma0, rho, nu0, nu1, delta, x1, x2, y00, y10, y11, w10, w11, k0, k1){
#   
#   sigma11 <- sigma0^2
#   
#   sigma12 <- rho*sigma0^2
#   
#   sigma13 <- rho*sigma0^2
#   
#   sigma23 <- rho*sigma0^2
#   
#   sigma22 <- sigma0^2
#   
#   sigma33 <- sigma0^2
#   
#   
#   gamma_0 <- sigma22*sigma33 - sigma23^2
#   gamma_1 <- sigma11*sigma33 - sigma13^2
#   gamma_2 <- sigma11*sigma22 - sigma12^2
#   gamma_3 <- sigma23*sigma13 - sigma12*sigma33
#   gamma_4 <- sigma12*sigma23 - sigma13*sigma22
#   gamma_5 <- sigma12*sigma13 - sigma11*sigma23
#   
#   
#   theta_denom <- sigma11*sigma22*sigma33 + 2*sigma12*sigma13*sigma23 - sigma12^2 * sigma33 - sigma13^2 * sigma22 - sigma23^2 * sigma11
#   
#   
#   theta_C <- gamma_0 + gamma_1 + gamma_2 + 2*gamma_3 + 2*gamma_4 + 2*gamma_5
#   
#   theta_D <- 2 * ( -gamma_0*y00 - gamma_1*y10 - gamma_2*y11 - gamma_3*y00 - gamma_3*y10 - gamma_4*y00 -gamma_4*y11 - gamma_5*y10 - gamma_5*y11 +
#                      gamma_1*w10*delta + gamma_2*w11*delta + gamma_3*w10*delta + gamma_4*w11*delta + gamma_5*w10*delta + gamma_5*w11*delta )
#   
#   theta_A <- (x2^2)*w11*theta_C
#   
#   theta_B <- 2*theta_C*x2*w11*(nu0 + x1*w10*nu1) + x2*w11*theta_D
#   
#   mu0_sd <- sqrt( 1 / (sum(theta_A / theta_denom) + 1/k0 ) )
#   
#   mu0_mean <- (-1) * ( sum(theta_B / theta_denom) - 2*k1/k0 ) / ( 2 * ( sum(theta_A / theta_denom) + 1/k0 ) )
#   
#   
#   mu0_t_draw <- rnorm(10000, mean = mu0_mean, sd = mu0_sd)
#   
#   return(mu0_t_draw)
#   
# }
# 
# 
# 
# # delta conditional posterior
# 
# delta_cond_posterior <- function(sigma0, rho, mu0, y00, y10, y11, w10, w11, v0, v1){
#   
#   sigma11 <- sigma0^2
#   
#   sigma12 <- rho*sigma0^2
#   
#   sigma13 <- rho*sigma0^2
#   
#   sigma23 <- rho*sigma0^2
#   
#   sigma22 <- sigma0^2
#   
#   sigma33 <- sigma0^2
#   
#   
#   gamma_0 <- sigma22*sigma33 - sigma23^2
#   gamma_1 <- sigma11*sigma33 - sigma13^2
#   gamma_2 <- sigma11*sigma22 - sigma12^2
#   gamma_3 <- sigma23*sigma13 - sigma12*sigma33
#   gamma_4 <- sigma12*sigma23 - sigma13*sigma22
#   gamma_5 <- sigma12*sigma13 - sigma11*sigma23
#   
#   
#   theta_denom <- sigma11*sigma22*sigma33 + 2*sigma12*sigma13*sigma23 - sigma12^2 * sigma33 - sigma13^2 * sigma22 - sigma23^2 * sigma11
#   
#   
#   theta_A <- gamma_1*w10 + gamma_2*w11 + 2*gamma_5*w10*w11
#   
#   theta_B <- 2 * ( gamma_1*mu0*w10 - gamma_1*y10*w10 - gamma_2*y11*w11 + gamma_2*mu0*w11 - gamma_3*y00*w10 + gamma_3*mu0*w10 -
#                      gamma_4*y00*w11 + gamma_4*mu0*w11 - gamma_5*y10*w11 - gamma_5*y11*w10 + gamma_5*mu0*w10 + gamma_5*mu0*w11)
#   
#   
#   delta_sd <- sqrt( 1 / (sum(theta_A / theta_denom) + 1/v0 ) )
#   
#   delta_mean <- (-1) * ( sum(theta_B / theta_denom) - 2*v1/v0 ) / ( 2 * (sum(theta_A / theta_denom) + 1/v0 ) )
#   
#   
#   delta_t_draw <- rnorm(10000, mean = delta_mean, sd = delta_sd)
#   
#   return(delta_t_draw)
#   
# }
# 
# 
# 
# 
# draws <- 10000
# 
# nu0_t_draw <- nu0_cond_posterior(sigma0_true, rho_true, nu1_true, nu2_true, delta_mean, x1, x2, potential_outcomes$y00, potential_outcomes$y10, potential_outcomes$y11, potential_outcomes$w10, potential_outcomes$w11, k0_mu0, k1_mu0)
# 
# nu0_mat <- data.frame(1:draws, nu0_t_draw)
# names(nu0_mat) <- c("draws", "nu0_vals")
# ggplot(nu0_mat, aes(draws, nu0_vals)) + geom_line()
# 
# 
# # nu1
# 
# nu1_t_draw <- nu1_cond_posterior(sigma0_true, rho_true, nu0_true, nu2_true, delta_mean, x1, x2, potential_outcomes$y00, potential_outcomes$y10, potential_outcomes$y11, potential_outcomes$w10, potential_outcomes$w11, k0_mu0, k1_mu0)
# 
# nu1_mat <- data.frame(1:draws, nu1_t_draw)
# names(nu1_mat) <- c("draws", "nu1_vals")
# ggplot(nu1_mat, aes(draws, nu1_vals)) + geom_line()
# 
# 
# # nu1
# 
# nu2_t_draw <- nu2_cond_posterior(sigma0_true, rho_true, nu0_true, nu1_true, delta_mean, x1, x2, potential_outcomes$y00, potential_outcomes$y10, potential_outcomes$y11, potential_outcomes$w10, potential_outcomes$w11, k0_mu0, k1_mu0)
# 
# nu2_mat <- data.frame(1:draws, nu2_t_draw)
# names(nu2_mat) <- c("draws", "nu2_vals")
# ggplot(nu2_mat, aes(draws, nu2_vals)) + geom_line()
# 
# 
# # update the mu0_t draw
# 
# mu0_true <- nu0_true + x1*potential_outcomes$w10*nu1_true + x2*potential_outcomes$w11*nu2_true
# # delta
# delta_t_draw <- delta_cond_posterior(sigma0_true, rho_true, mu0_true, potential_outcomes$y00, potential_outcomes$y10, potential_outcomes$y11, potential_outcomes$w10, potential_outcomes$w11, v0_delta, v1_delta)
# 
# 
# delta_mat <- data.frame(1:draws, delta_t_draw)
# names(delta_mat) <- c("draws", "delta_vals")
# ggplot(delta_mat, aes(draws, delta_vals)) + geom_line()
# 
# 
# 
# 
# 
# 
