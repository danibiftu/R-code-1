###################### R code for simulation of latent trait model with Bayesian rank likelihood method########################
# Load necessary libraries
library(LaplacesDemon)
library(geoR)
library(mvtnorm)
library(MCMCpack)
library(msm)
# Set seed for reproducibility
set.seed(12345)
# Define constants
n <- 300  # sample size
groups <- 11
beta <- matrix(c(-0.9, 0.4, 0.6, 0.2, 0.1), ncol = 1)
# Initialize variables
epsilon <- matrix(NA, n, 1)
x_m <- list()
# Create 11 design matrices with the covariates
for (i in 1:groups) {
  x_m[[i]] <- matrix(NA, nrow = n, ncol = length(beta))
  for (j in 1:length(beta)) {
    x_m[[i]][, 1] <- rnorm(n, 15, 0.1)
    x_m[[i]][, 2] <- rpois(n, 2)
    x_m[[i]][, 3] <- rnorm(n, 23, 0.2)
    x_m[[i]][, 4] <- rpois(n, 2)
    x_m[[i]][, 5] <- rnorm(n, 16, 0.2)
  }
}
# Generate epsilon
epsilon <- rnorm(n, 0, 1)
# Initialize random effects and latent variables
sigma2_random_effects <- 0.83
lamda2 <- 0.5
omega <- 5
random_effects <- numeric(groups)
latents <- latents_cat <- matrix(NA, n, groups)
for (j in 1:groups) {
  random_effects[j] <- rnorm(1, 0, sigma2_random_effects)
  latents[, j] <- x_m[[j]] %*% beta + random_effects[j] + epsilon
}
# Convert latents to categorical data
latents_cat <- ceiling(latents)
for (i in 1:nrow(latents_cat)) {
  for (j in 1:ncol(latents_cat)) {
    if (latents_cat[i, j] >= -1 & latents_cat[i, j] < 2) {
      latents_cat[i, j] <- (-1 + 2) / 2
    } else if (latents_cat[i, j] >= 2 & latents_cat[i, j] < 5) {
      latents_cat[i, j] <- (0 + 3) / 2
    } else if (latents_cat[i, j] >= 3 & latents_cat[i, j] < 6) {
      latents_cat[i, j] <- (3 + 6) / 2
    } else if (latents_cat[i, j] >= 6 & latents_cat[i, j] < 9) {
      latents_cat[i, j] <- (6 + 9) / 2
    }
  }
}
# Plot the first group of latents_cat
plot(latents_cat[, 1], type = "l")
# Prepare data for MCMC
m <- 11
n_obs <- n - 1
Y <- latents_cat[1:n_obs, ]
X <- array(NA, c(n_obs, length(beta), m))
X_init <- array(NA, c(n, length(beta), m))
for (i in 1:groups) {
  x_m[[i]] <- scale(x_m[[i]], center = TRUE, scale = TRUE)
  X[,, i] <- x_m[[i]][1:n_obs, ]
  X_init[,, i] <- x_m[[i]]
}
# Define prior hyperparameters
lamda2_prior <- 2
omega_upsilon_prior <- 15
BETA <- matrix(0, ncol(X[,,1]), 1)
z <- matrix(0.5, nrow = n_obs, ncol = m)
a <- b <- matrix(NA, nrow = n_obs, ncol = m)
# Initialize MCMC variables
n_iter <- 50000
BETA_POST_4 <- NULL
upsilon_POST_4 <- NULL
sigmasq_upsilon_POST_4 <- NULL
z_all_4 <- NULL
# Run Gibbs Sampler
for (i in 1:n_iter) {
  omega_upsilon_posterior <- rep(omega_upsilon_prior + 1, m)
  lamda2_posterior <- lamda2_prior + (upsilon_POST_5^2) / omega_upsilon_posterior
  sigmasq_upsilon_post <- rinvchisq(m, omega_upsilon_posterior, lamda2_posterior)
  mean_upsilon_posterior <- array(NA, c(n_obs, 1, m))
  for (t in 1:n_obs) {
    for (j in 1:m) {
      mean_upsilon_posterior[t, , j] <- z[t, j] - X[t, , j] %*% BETA
    }
  }
  var_upsilon_posterior <- 1 / (n_obs + (1 / sigmasq_upsilon_post))
  mean_upsilon_posterior_sum <- apply(mean_upsilon_posterior, 3, sum) / (n_obs + 1 / sigmasq_upsilon_post)
  upsilon_posterior <- rnorm(m, mean_upsilon_posterior_sum, sqrt(var_upsilon_posterior))
  upsilon_rep <- array(rep(upsilon_posterior, each = n_obs), c(n_obs, 1, m))
  for (t in 1:n_obs) {
    for (j in 1:m) {
      zt <- z[t,]
      a[t, j] <- max(zt[Y[t, ] < Y[t, j]])
      b[t, j] <- min(zt[Y[t, j] < Y[t, ]])
      z[t, j] <- rtnorm(1, X[t, , j] %*% BETA + upsilon_rep[t, 1, j], 1, a[t, j], b[t, j])
    }
  }
  muj1 <- array(NA, c(length(beta), length(beta), m))
  muj1_sum <- matrix(0, length(beta), length(beta))
  muj2 <- array(NA, c(length(beta), 1, m))
  muj2_sum <- matrix(0, length(beta), 1) 
  for (j in 1:m) {
    muj1[,, j] <- t(X[,, j]) %*% X[,, j]
    muj1_sum <- muj1_sum + muj1[,, j]
    muj2[,, j] <- t(X[,, j]) %*% (z[, j] - upsilon_rep[, 1, j])
    muj2_sum <- muj2_sum + muj2[,, j]
  }
  BETA <- mvrnorm(1, solve(muj1_sum) %*% muj2_sum, solve(muj1_sum))
  BETA_POST_4 <- rbind(BETA_POST_4, BETA)
  upsilon_POST_4 <- rbind(upsilon_POST_4, upsilon_posterior)
  sigmasq_upsilon_POST_4 <- rbind(sigmasq_upsilon_POST_4, sigmasq_upsilon_post)
  print(paste("Iteration", i))
  z_all_4 <- rbind(z_all_4, z)
}
# Sensitivity Analysis for BETA
# Load required packages
install.packages("viridis")
library(bayesplot)
library(viridis)
# Prepare the data for beta chains
beta1 <- cbind(BETA_POST_1[1:50000, 1], BETA_POST_2[1:50000, 1], BETA_POST_3[1:50000, 1], BETA_POST_4[1:50000, 1], BETA_POST_5[1:50000, 1])
beta2 <- cbind(BETA_POST_1[1:50000, 2], BETA_POST_2[1:50000, 2], BETA_POST_3[1:50000, 2], BETA_POST_4[1:50000, 2], BETA_POST_5[1:50000, 2])
beta3 <- cbind(BETA_POST_1[1:50000, 3], BETA_POST_2[1:50000, 3], BETA_POST_3[1:50000, 3], BETA_POST_4[1:50000, 3], BETA_POST_5[1:50000, 3])
beta4 <- cbind(BETA_POST_1[1:50000, 4], BETA_POST_2[1:50000, 4], BETA_POST_3[1:50000, 4], BETA_POST_4[1:50000, 4], BETA_POST_5[1:50000, 4])
beta5 <- cbind(BETA_POST_1[1:50000, 5], BETA_POST_2[1:50000, 5], BETA_POST_3[1:50000, 5], BETA_POST_4[1:50000, 5], BETA_POST_5[1:50000, 5])
chains_beta <- array(cbind(beta1, beta2, beta3, beta4, beta5), c(50000, 5, 5))
# Set column and parameter names
dimnames(chains_beta)[[2]] <- c("chain1", "chain2", "chain3", "chain4", "chain5")
dimnames(chains_beta)[[3]] <- c("\u03B21", "\u03B22", "\u03B23", "\u03B24", "\u03B25")
# Plot MCMC sensitivity trace for beta
color_scheme_set("viridis")
mcmc_trace(chains_beta, pars = c("\u03B21", "\u03B22", "\u03B23", "\u03B24", "\u03B25"), 
           facet_args = list(ncol = 1, strip.position = "left"))
# Prepare the data for upsilon chains
upsilon1 <- cbind(upsilon_POST_1[1:50000, 1], upsilon_POST_2[1:50000, 1], upsilon_POST_3[1:50000, 1], upsilon_POST_4[1:50000, 1], upsilon_POST_5[1:50000, 1])
upsilon2 <- cbind(upsilon_POST_1[1:50000, 2], upsilon_POST_2[1:50000, 2], upsilon_POST_3[1:50000, 2], upsilon_POST_4[1:50000, 2], upsilon_POST_5[1:50000, 2])
upsilon3 <- cbind(upsilon_POST_1[1:50000, 3], upsilon_POST_2[1:50000, 3], upsilon_POST_3[1:50000, 3], upsilon_POST_4[1:50000, 3], upsilon_POST_5[1:50000, 3])
upsilon4 <- cbind(upsilon_POST_1[1:50000, 4], upsilon_POST_2[1:50000, 4], upsilon_POST_3[1:50000, 4], upsilon_POST_4[1:50000, 4], upsilon_POST_5[1:50000, 4])
upsilon5 <- cbind(upsilon_POST_1[1:50000, 5], upsilon_POST_2[1:50000, 5], upsilon_POST_3[1:50000, 5], upsilon_POST_4[1:50000, 5], upsilon_POST_5[1:50000, 5])
upsilon6 <- cbind(upsilon_POST_1[1:50000, 6], upsilon_POST_2[1:50000, 6], upsilon_POST_3[1:50000, 6], upsilon_POST_4[1:50000, 6], upsilon_POST_5[1:50000, 6])
upsilon7 <- cbind(upsilon_POST_1[1:50000, 7], upsilon_POST_2[1:50000, 7], upsilon_POST_3[1:50000, 7], upsilon_POST_4[1:50000, 7], upsilon_POST_5[1:50000, 7])
upsilon8 <- cbind(upsilon_POST_1[1:50000, 8], upsilon_POST_2[1:50000, 8], upsilon_POST_3[1:50000, 8], upsilon_POST_4[1:50000, 8], upsilon_POST_5[1:50000, 8])
upsilon9 <- cbind(upsilon_POST_1[1:50000, 9], upsilon_POST_2[1:50000, 9], upsilon_POST_3[1:50000, 9], upsilon_POST_4[1:50000, 9], upsilon_POST_5[1:50000, 9])
upsilon10 <- cbind(upsilon_POST_1[1:50000, 10], upsilon_POST_2[1:50000, 10], upsilon_POST_3[1:50000, 10], upsilon_POST_4[1:50000, 10], upsilon_POST_5[1:50000, 10])
upsilon11 <- cbind(upsilon_POST_1[1:50000, 11], upsilon_POST_2[1:50000, 11], upsilon_POST_3[1:50000, 11], upsilon_POST_4[1:50000, 11], upsilon_POST_5[1:50000, 11])
chains_upsilons <- array(cbind(upsilon1, upsilon2, upsilon3, upsilon4, upsilon5, 
                               upsilon6, upsilon7, upsilon8, upsilon9, upsilon10, upsilon11), c(50000, 5, 11))
# Set column and parameter names for upsilons
dimnames(chains_upsilons)[[2]] <- c("chain1", "chain2", "chain3", "chain4", "chain5")
dimnames(chains_upsilons)[[3]] <- c("\u03C51", "\u03C52", "\u03C53", "\u03C54", "\u03C55", 
                                    "\u03C56", "\u03C57", "\u03C58", "\u03C59", "\u03C510", "\u03C511")
# Plot MCMC sensitivity trace for upsilons
mcmc_trace(chains_upsilons, pars = c("\u03C51", "\u03C52", "\u03C53", "\u03C54", "\u03C55", 
                                     "\u03C56", "\u03C57", "\u03C58", "\u03C59", "\u03C510", "\u03C511"), 
           facet_args = list(ncol = 2, strip.position = "left"))
# Compute RMSE & MAE for beta  - CHAIN 1
posterior_means_beta_1 <- matrix(colMeans(BETA_POST_1), 5, 1)
mse_beta_1 <- matrix(NA, 5, 1)
for (j in 1:5) {
  mse_beta_1[j] <- (posterior_means_beta_1[j, 1] - beta[j, 1])^2
}
MSE_beta_1 <- mean(mse_beta_1)
RMSE_beta_1 <- sqrt(MSE_beta_1)
mae_beta_1 <- matrix(NA, 5, 1)
for (j in 1:5) {
  mae_beta_1[j] <- abs(posterior_means_beta_1[j, 1] - beta[j, ])
}
MAE_beta_1 <- mean(mae_beta_1)
# Compute RMSE & MAE for upsilon - CHAIN 1
posterior_means_upsilon_1 <- matrix(colMeans(upsilon_POST_1), m, 1)
mse_upsilon_1 <- matrix(NA, m, 1)
for (j in 1:m) {
  mse_upsilon_1[j] <- (posterior_means_upsilon_1[j, 1] - upsilon[j, 1])^2
}
MSE_upsilon_1 <- mean(mse_upsilon_1)
RMSE_upsilon_1 <- sqrt(MSE_upsilon_1)
mae_upsilon_1 <- matrix(NA, m, 1)
for (j in 1:m) {
  mae_upsilon_1[j] <- abs(posterior_means_upsilon_1[j, 1] - upsilon[j, 1])
}
MAE_upsilon_1 <- mean(mae_upsilon_1)
# Compute RMSE & MAE for beta - CHAIN 2
posterior_means_beta_2 <- matrix(colMeans(BETA_POST_2), 5, 1)
mse_beta_2 <- matrix(NA, 5, 1)
for (j in 1:5) {
  mse_beta_2[j] <- (posterior_means_beta_2[j, 1] - beta[j, 1])^2
}
MSE_beta_2 <- mean(mse_beta_2)
RMSE_beta_2 <- sqrt(MSE_beta_2)
mae_beta_2 <- matrix(NA, 5, 1)
for (j in 1:5) {
  mae_beta_2[j] <- abs(posterior_means_beta_2[j, 1] - beta[j, ])
}
MAE_beta_2 <- mean(mae_beta_2)
# Compute RMSE & MAE for upsilon - CHAIN 2
posterior_means_upsilon_2 <- matrix(colMeans(upsilon_POST_2), m, 1)
mse_upsilon_2 <- matrix(NA, m, 1)
for (j in 1:m) {
  mse_upsilon_2[j] <- (posterior_means_upsilon_2[j, 1] - upsilon[j, 1])^2
}
MSE_upsilon_2 <- mean(mse_upsilon_2)
RMSE_upsilon_2 <- sqrt(MSE_upsilon_2)

mae_upsilon_2 <- matrix(NA, m, 1)
for (j in 1:m) {
  mae_upsilon_2[j] <- abs(posterior_means_upsilon_2[j, 1] - upsilon[j, 1])
}
MAE_upsilon_2 <- mean(mae_upsilon_2)
# Compute RMSE & MAE for beta - CHAIN 3
posterior_means_beta_3 <- matrix(colMeans(BETA_POST_3), 5, 1)
mse_beta_3 <- matrix(NA, 5, 1)
for (j in 1:5) {
  mse_beta_3[j] <- (posterior_means_beta_3[j, 1] - beta[j, 1])^2
}
MSE_beta_3 <- mean(mse_beta_3)
RMSE_beta_3 <- sqrt(MSE_beta_3)
mae_beta_3 <- matrix(NA, 5, 1)
for (j in 1:5) {
  mae_beta_3[j] <- abs(posterior_means_beta_3[j, 1] - beta[j, ])
}
MAE_beta_3 <- mean(mae_beta_3)
# Compute RMSE & MAE for upsilon - CHAIN 3
posterior_means_upsilon_3 <- matrix(colMeans(upsilon_POST_3), m, 1)
mse_upsilon_3 <- matrix(NA, m, 1)
for (j in 1:m) {
  mse_upsilon_3[j] <- (posterior_means_upsilon_3[j, 1] - upsilon[j, 1])^2
}
MSE_upsilon_3 <- mean(mse_upsilon_3)
RMSE_upsilon_3 <- sqrt(MSE_upsilon_3)
mae_upsilon_3 <- matrix(NA, m, 1)
for (j in 1:m) {
  mae_upsilon_3[j] <- abs(posterior_means_upsilon_3[j, 1] - upsilon[j, 1])
}
MAE_upsilon_3 <- mean(mae_upsilon_3)
# Compute RMSE & MAE for beta - CAINH 4
posterior_means_beta_4 <- matrix(colMeans(BETA_POST_4), 5, 1)
mse_beta_4 <- matrix(NA, 5, 1)
for (j in 1:5) {
  mse_beta_4[j] <- (posterior_means_beta_4[j, 1] - beta[j, 1])^2
}
MSE_beta_4 <- mean(mse_beta_4)
RMSE_beta_4 <- sqrt(MSE_beta_4)
mae_beta_4 <- matrix(NA, 5, 1)
for (j in 1:5) {
  mae_beta_4[j] <- abs(posterior_means_beta_4[j, 1] - beta[j, ])
}
MAE_beta_4 <- mean(mae_beta_4)
# Compute RMSE & MAE for upsilon - Chain 4
posterior_means_upsilon_4 <- matrix(colMeans(upsilon_POST_4), m, 1)
mse_upsilon_4 <- matrix(NA, m, 1)
for (j in 1:m) {
  mse_upsilon_4[j] <- (posterior_means_upsilon_4[j, 1] - upsilon[j, 1])^2
}
MSE_upsilon_4 <- mean(mse_upsilon_4)
RMSE_upsilon_4 <- sqrt(MSE_upsilon_4)
mae_upsilon_4 <- matrix(NA, m, 1)
for (j in 1:m) {
  mae_upsilon_4[j] <- abs(posterior_means_upsilon_4[j, 1] - upsilon[j, 1])
}
MAE_upsilon_4 <- mean(mae_upsilon_4)
# Compute RMSE & MAE for beta - CHAIN 5
posterior_means_beta_5 <- matrix(colMeans(BETA_POST_5), 5, 1)
mse_beta_5 <- matrix(NA, 5, 1)
for (j in 1:5) {
  mse_beta_5[j] <- (posterior_means_beta_5[j, 1] - beta[j, 1])^2
}
MSE_beta_5 <- mean(mse_beta_5)
RMSE_beta_5 <- sqrt(MSE_beta_5)
mae_beta_5 <- matrix(NA, 5, 1)
for (j in 1:5) {
  mae_beta_5[j] <- abs(posterior_means_beta_5[j, 1] - beta[j, ])
}
MAE_beta_5 <- mean(mae_beta_5)
# Compute RMSE & MAE for upsilon - CHAIN 5
posterior_means_upsilon_5 <- matrix(colMeans(upsilon_POST_5), m, 1)
mse_upsilon_5 <- matrix(NA, m, 1)
for (j in 1:m) {
  mse_upsilon_5[j] <- (posterior_means_upsilon_5[j, 1] - upsilon[j, 1])^2
}
MSE_upsilon_5 <- mean(mse_upsilon_5)
RMSE_upsilon_5 <- sqrt(MSE_upsilon_5)
mae_upsilon_5 <- matrix(NA, m, 1)
for (j in 1:m) {
  mae_upsilon_5[j] <- abs(posterior_means_upsilon_5[j, 1] - upsilon[j, 1])
}
MAE_upsilon_5 <- mean(mae_upsilon_5)
########## Gelman-Rubin Convergence Test ##########
# Load required libraries
library(stableGR)
library(coda)
library(bayesplot)
# Calculate PSRF for BETA parameters
PSRF_1_beta_test <- stable.GR(BETA_POST_1)
PSRF_2_beta_test <- stable.GR(BETA_POST_2)
PSRF_3_beta_test <- stable.GR(BETA_POST_3)
PSRF_4_beta_test <- stable.GR(BETA_POST_4)
PSRF_5_beta_test <- stable.GR(BETA_POST_5)
# Calculate PSRF for upsilon parameters
PSRF_1_upsilon_test <- stable.GR(upsilon_POST_1)
PSRF_2_upsilon_test <- stable.GR(upsilon_POST_2)
PSRF_3_upsilon_test <- stable.GR(upsilon_POST_3)
PSRF_4_upsilon_test <- stable.GR(upsilon_POST_4)
PSRF_5_upsilon_test <- stable.GR(upsilon_POST_5)
# 95% Credible interval & coverage probability of beta chain 1
Beta_mcmc_object_1 <- as.mcmc(BETA_POST_1)
Beta_credible_interval_1 <- HPDinterval(Beta_mcmc_object_1, prob = 0.95)
coverage_prob_beta_1 <- matrix(NA, 5, 1)
for (j in 1:5) {
  if (beta[j, 1] >= Beta_credible_interval_1[j, 1] && beta[j, 1] <= Beta_credible_interval_1[j, 2]) {
    coverage_prob_beta_1[j, 1] <- 1
  } else {
    coverage_prob_beta_1[j, 1] <- 0
  }
}
mean_coverage_prob_beta_1 <- mean(coverage_prob_beta_1)
# 95% Credible interval & coverage probability of upsilon chain 1
upsilon_mcmc_object_1 <- as.mcmc(upsilon_POST_1)
upsilon_credible_interval_1 <- HPDinterval(upsilon_mcmc_object_1, prob = 0.95)
coverage_prob_upsilon_1 <- matrix(NA, m, 1)
for (j in 1:m) {
  if (upsilon[j, 1] >= upsilon_credible_interval_1[j, 1] && upsilon[j, 1] <= upsilon_credible_interval_1[j, 2]) {
    coverage_prob_upsilon_1[j, 1] <- 1
  } else {
    coverage_prob_upsilon_1[j, 1] <- 0
  }
}
mean_coverage_prob_upsilon_1 <- mean(coverage_prob_upsilon_1)
# 95% Credible interval & coverage probability of beta chain 2
Beta_mcmc_object_2 <- as.mcmc(BETA_POST_2)
Beta_credible_interval_2 <- HPDinterval(Beta_mcmc_object_2, prob = 0.95)
coverage_prob_beta_2 <- matrix(NA, 5, 1)
for (j in 1:5) {
  if (beta[j, 1] >= Beta_credible_interval_2[j, 1] && beta[j, 1] <= Beta_credible_interval_2[j, 2]) {
    coverage_prob_beta_2[j, 1] <- 1
  } else {
    coverage_prob_beta_2[j, 1] <- 0
  }
}
mean_coverage_prob_beta_2 <- mean(coverage_prob_beta_2)
# 95% Credible interval & coverage probability of upsilon chain 2
upsilon_mcmc_object_2 <- as.mcmc(upsilon_POST_2)
upsilon_credible_interval_2 <- HPDinterval(upsilon_mcmc_object_2, prob = 0.95)
coverage_prob_upsilon_2 <- matrix(NA, m, 1)
for (j in 1:m) {
  if (upsilon[j, 1] >= upsilon_credible_interval_2[j, 1] && upsilon[j, 1] <= upsilon_credible_interval_2[j, 2]) {
    coverage_prob_upsilon_2[j, 1] <- 1
  } else {
    coverage_prob_upsilon_2[j, 1] <- 0
  }
}
mean_coverage_prob_upsilon_2 <- mean(coverage_prob_upsilon_2)
# 95% Credible interval & coverage probability of beta chain 3
Beta_mcmc_object_3 <- as.mcmc(BETA_POST_3)
Beta_credible_interval_3 <- HPDinterval(Beta_mcmc_object_3, prob = 0.95)
coverage_prob_beta_3 <- matrix(NA, 5, 1)
for (j in 1:5) {
  if (beta[j, 1] >= Beta_credible_interval_3[j, 1] && beta[j, 1] <= Beta_credible_interval_3[j, 2]) {
    coverage_prob_beta_3[j, 1] <- 1
  } else {
    coverage_prob_beta_3[j, 1] <- 0
  }
}
mean_coverage_prob_beta_3 <- mean(coverage_prob_beta_3)
# 95% Credible interval & coverage probability of upsilon chain 3
upsilon_mcmc_object_3 <- as.mcmc(upsilon_POST_3)
upsilon_credible_interval_3 <- HPDinterval(upsilon_mcmc_object_3, prob = 0.95)
coverage_prob_upsilon_3 <- matrix(NA, m, 1)
for (j in 1:m) {
  if (upsilon[j, 1] >= upsilon_credible_interval_3[j, 1] && upsilon[j, 1] <= upsilon_credible_interval_3[j, 2]) {
    coverage_prob_upsilon_3[j, 1] <- 1
  } else {
    coverage_prob_upsilon_3[j, 1] <- 0
  }
}
mean_coverage_prob_upsilon_3 <- mean(coverage_prob_upsilon_3)
# 95% Credible interval & coverage probability of beta chain 4
Beta_mcmc_object_4 <- as.mcmc(BETA_POST_4)
Beta_credible_interval_4 <- HPDinterval(Beta_mcmc_object_4, prob = 0.95)
coverage_prob_beta_4 <- matrix(NA, 5, 1)
for (j in 1:5) {
  if (beta[j, 1] >= Beta_credible_interval_4[j, 1] && beta[j, 1] <= Beta_credible_interval_4[j, 2]) {
    coverage_prob_beta_4[j, 1] <- 1
  } else {
    coverage_prob_beta_4[j, 1] <- 0
  }
}
mean_coverage_prob_beta_4 <- mean(coverage_prob_beta_4)
# 95% Credible interval & coverage probability of upsilon chain 4
upsilon_mcmc_object_4 <- as.mcmc(upsilon_POST_4)
upsilon_credible_interval_4 <- HPDinterval(upsilon_mcmc_object_4, prob = 0.95)
coverage_prob_upsilon_4 <- matrix(NA, m, 1)
for (j in 1:m) {
  if (upsilon[j, 1] >= upsilon_credible_interval_4[j, 1] && upsilon[j, 1] <= upsilon_credible_interval_4[j, 2]) {
    coverage_prob_upsilon_4[j, 1] <- 1
  } else {
    coverage_prob_upsilon_4[j, 1] <- 0
  }
}
mean_coverage_prob_upsilon_4 <- mean(coverage_prob_upsilon_4)
# 95% Credible interval & coverage probability of beta chain 5
Beta_mcmc_object_5 <- as.mcmc(BETA_POST_5)
Beta_credible_interval_5 <- HPDinterval(Beta_mcmc_object_5, prob = 0.95)
coverage_prob_beta_5 <- matrix(NA, 5, 1)
for (j in 1:5) {
  if (beta[j, 1] >= Beta_credible_interval_5[j, 1] && beta[j, 1] <= Beta_credible_interval_5[j, 2]) {
    coverage_prob_beta_5[j, 1] <- 1
  } else {
    coverage_prob_beta_5[j, 1] <- 0
  }
}
mean_coverage_prob_beta_5 <- mean(coverage_prob_beta_5)
# 95% Credible interval & coverage probability of upsilon chain 5
upsilon_mcmc_object_5 <- as.mcmc(upsilon_POST_5)
upsilon_credible_interval_5 <- HPDinterval(upsilon_mcmc_object_5, prob = 0.95)
coverage_prob_upsilon_5 <- matrix(NA, m, 1)
for (j in 1:m) {
  if (upsilon[j, 1] >= upsilon_credible_interval_5[j, 1] && upsilon[j, 1] <= upsilon_credible_interval_5[j, 2]) {
    coverage_prob_upsilon_5[j, 1] <- 1
  } else {
    coverage_prob_upsilon_5[j, 1] <- 0
  }
}
mean_coverage_prob_upsilon_5 <- mean(coverage_prob_upsilon_5)
# Checking the convergence after burning period
z_gibbis_burn <- z_all_4[seq(1, nrow(z_all_4[10001:4950000,]), 100),]
Beta_gibbis_burn <- BETA_POST_4[seq(1, nrow(BETA_POST_4[10001:50000,]), 10),]
upsilon_gibbis_burn <- upsilon_POST_4[seq(1, nrow(upsilon_POST_4[10001:50000,]), 10),]
sigma_gibbis_burn <- sigmasq_upsilon_POST_4[seq(1, nrow(sigmasq_upsilon_POST_4[10001:50000,]), 10),]
# Trace plots for Beta
par(mfrow=c(5,1), mar=c(3,3,2,1))
traceplot(as.mcmc(Beta_gibbis_burn[,1]), main=expression(paste(beta[1])))
traceplot(as.mcmc(Beta_gibbis_burn[,2]), main=expression(paste(beta[2])))
traceplot(as.mcmc(Beta_gibbis_burn[,3]), main=expression(paste(beta[3])))
traceplot(as.mcmc(Beta_gibbis_burn[,4]), main=expression(paste(beta[4])))
traceplot(as.mcmc(Beta_gibbis_burn[,5]), main=expression(paste(beta[5])))
# Density plots for Beta
par(mfrow=c(3,2), mar=c(7,3,3,1))
densplot(as.mcmc(Beta_gibbis_burn[,1]), main=expression(paste(beta[1])))
densplot(as.mcmc(Beta_gibbis_burn[,2]), main=expression(paste(beta[2])))
densplot(as.mcmc(Beta_gibbis_burn[,3]), main=expression(paste(beta[3])))
densplot(as.mcmc(Beta_gibbis_burn[,4]), main=expression(paste(beta[4])))
densplot(as.mcmc(Beta_gibbis_burn[,5]), main=expression(paste(beta[5])))
# Autocorrelation plots for Beta
par(mfrow=c(3,2), mar=c(7,4,3,1))
acf(as.mcmc(Beta_gibbis_burn[,1]), ylab="Autocorrelation", ci=F, lag.max=4000, main=expression(paste(beta[1])))
acf(as.mcmc(Beta_gibbis_burn[,2]), ylab="Autocorrelation", ci=F, lag.max=4000, main=expression(paste(beta[2])))
acf(as.mcmc(Beta_gibbis_burn[,3]), ylab="Autocorrelation", ci=F, lag.max=4000, main=expression(paste(beta[3])))
acf(as.mcmc(Beta_gibbis_burn[,4]), ylab="Autocorrelation", ci=F, lag.max=4000, main=expression(paste(beta[4])))
acf(as.mcmc(Beta_gibbis_burn[,5]), ylab="Autocorrelation", ci=F, lag.max=4000, main=expression(paste(beta[5])))
# Posterior means of Beta
mean(Beta_gibbis_burn[,1])
mean(Beta_gibbis_burn[,2])
mean(Beta_gibbis_burn[,3]) 
mean(Beta_gibbis_burn[,4])
mean(Beta_gibbis_burn[,5])
# Posterior standard deviations of Beta
sd(Beta_gibbis_burn[,1])
sd(Beta_gibbis_burn[,2])
sd(Beta_gibbis_burn[,3]) 
sd(Beta_gibbis_burn[,4])
sd(Beta_gibbis_burn[,5])
# Trace plots for upsilon
par(mfrow=c(6,1),mar=c(3,3,2,2))
traceplot(as.mcmc(upsilon_gibbis_burn[,1]), main=expression(paste(upsilon[1])))
traceplot(as.mcmc(upsilon_gibbis_burn[,2]), main=expression(paste(upsilon[2])))
traceplot(as.mcmc(upsilon_gibbis_burn[,3]), main=expression(paste(upsilon[3])))
traceplot(as.mcmc(upsilon_gibbis_burn[,4]), main=expression(paste(upsilon[4])))
traceplot(as.mcmc(upsilon_gibbis_burn[,5]), main=expression(paste(upsilon[5])))
traceplot(as.mcmc(upsilon_gibbis_burn[,6]), main=expression(paste(upsilon[6])))
par(mfrow=c(5,1),mar=c(3,3,2,2))
traceplot(as.mcmc(upsilon_gibbis_burn[,7]), main=expression(paste(upsilon[7])))
traceplot(as.mcmc(upsilon_gibbis_burn[,8]), main=expression(paste(upsilon[8])))
traceplot(as.mcmc(upsilon_gibbis_burn[,9]), main=expression(paste(upsilon[9])))
traceplot(as.mcmc(upsilon_gibbis_burn[,10]), main=expression(paste(upsilon[10])))
traceplot(as.mcmc(upsilon_gibbis_burn[,11]), main=expression(paste(upsilon[11])))
# Density plots for upsilon
par(mfrow=c(3,2),mar=c(7,3,3,1))
densplot(as.mcmc(upsilon_gibbis_burn[,1]), main=expression(paste(upsilon[1])))
densplot(as.mcmc(upsilon_gibbis_burn[,2]), main=expression(paste(upsilon[2])))
densplot(as.mcmc(upsilon_gibbis_burn[,3]), main=expression(paste(upsilon[3])))
densplot(as.mcmc(upsilon_gibbis_burn[,4]), main=expression(paste(upsilon[4])))
densplot(as.mcmc(upsilon_gibbis_burn[,5]), main=expression(paste(upsilon[5])))
densplot(as.mcmc(upsilon_gibbis_burn[,6]), main=expression(paste(upsilon[6])))
par(mfrow=c(3,2),mar=c(7,3,3,1))
densplot(as.mcmc(upsilon_gibbis_burn[,7]), main=expression(paste(upsilon[7])))
densplot(as.mcmc(upsilon_gibbis_burn[,8]), main=expression(paste(upsilon[8])))
densplot(as.mcmc(upsilon_gibbis_burn[,9]), main=expression(paste(upsilon[9])))
densplot(as.mcmc(upsilon_gibbis_burn[,10]), main=expression(paste(upsilon[10])))
densplot(as.mcmc(upsilon_gibbis_burn[,11]), main=expression(paste(upsilon[11])))
# Autocorrelation plots for upsilon
par(mfrow=c(3,2),mar=c(7,4,3,1))
acf(as.mcmc(upsilon_gibbis_burn[,1]), ylab="Autocorrelation", ci=F, lag.max=4000, main=expression(paste(upsilon[1])))
acf(as.mcmc(upsilon_gibbis_burn[,2]), ylab="Autocorrelation", ci=F, lag.max=4000, main=expression(paste(upsilon[2])))
acf(as.mcmc(upsilon_gibbis_burn[,3]), ylab="Autocorrelation", ci=F, lag.max=4000, main=expression(paste(upsilon[3])))
acf(as.mcmc(upsilon_gibbis_burn[,4]), ylab="Autocorrelation", ci=F, lag.max=4000, main=expression(paste(upsilon[4])))
acf(as.mcmc(upsilon_gibbis_burn[,5]), ylab="Autocorrelation", ci=F, lag.max=4000, main=expression(paste(upsilon[5])))
acf(as.mcmc(upsilon_gibbis_burn[,6]), ylab="Autocorrelation", ci=F, lag.max=4000, main=expression(paste(upsilon[6])))
par(mfrow=c(3,2),mar=c(7,4,3,1))
acf(as.mcmc(upsilon_gibbis_burn[,7]), ylab="Autocorrelation", ci=F, lag.max=4000, main=expression(paste(upsilon[7])))
acf(as.mcmc(upsilon_gibbis_burn[,8]), ylab="Autocorrelation", ci=F, lag.max=4000, main=expression(paste(upsilon[8])))
acf(as.mcmc(upsilon_gibbis_burn[,9]), ylab="Autocorrelation", ci=F, lag.max=4000, main=expression(paste(upsilon[9])))
acf(as.mcmc(upsilon_gibbis_burn[,10]), ylab="Autocorrelation", ci=F, lag.max=4000, main=expression(paste(upsilon[10])))
acf(as.mcmc(upsilon_gibbis_burn[,11]), ylab="Autocorrelation", ci=F, lag.max=4000, main=expression(paste(upsilon[11])))
# Posterior means of upsilon
mean(upsilon_gibbis_burn[,1])
mean(upsilon_gibbis_burn[,2])
mean(upsilon_gibbis_burn[,3])
mean(upsilon_gibbis_burn[,4])
mean(upsilon_gibbis_burn[,5])
mean(upsilon_gibbis_burn[,6])
mean(upsilon_gibbis_burn[,7])
mean(upsilon_gibbis_burn[,8])
mean(upsilon_gibbis_burn[,9])
mean(upsilon_gibbis_burn[,10])
mean(upsilon_gibbis_burn[,11])
#Posterior standard deviations of upsilon
sd(upsilon_gibbis_burn[,1])
sd(upsilon_gibbis_burn[,2])
sd(upsilon_gibbis_burn[,3])
sd(upsilon_gibbis_burn[,4])
sd(upsilon_gibbis_burn[,5])
sd(upsilon_gibbis_burn[,6])
sd(upsilon_gibbis_burn[,7])
sd(upsilon_gibbis_burn[,8])
sd(upsilon_gibbis_burn[,9])
sd(upsilon_gibbis_burn[,10])
sd(upsilon_gibbis_burn[,11])
# Load required libraries
library(reshape2)
library(viridis)
library(ggplot2)
library(bayesplot)
# Renaming column names
new_colnames <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")
colnames(g) <- new_colnames
g <- data.frame(upsilon_POST_4)
# Calculate and plot the variation between groups
data <- var(g[sapply(g, is.numeric)])
data1 <- melt(data)
colnames(data1) <- c("Var1", "Var2", "Variations")
x_labels <- y_labels <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")
ggplot(data1, aes(x = Var1, y = Var2, fill = Variations)) +
  geom_tile() + scale_fill_viridis(discrete = FALSE) + geom_tile() +
  labs(title = "", x = "Regions", y = "Regions") +
  scale_x_discrete(labels = x_labels) +
  scale_y_discrete(labels = y_labels)
# Posterior predictive check
last_obs_1 <- array(NA, c(ncol(X), m, 1))
for (j in 1:m) {
  last_obs_1[, j, 1] <- X_init[(n_obs + 1), , j]
}
z_pred_outsample_1 <- array(NA, c(n_iter, 1, m))
for (i in 1:1) {
  for (j in 1:m) {
    z_pred_outsample_1[, i, j] <- rnorm(n_iter, BETA_POST_5 %*% last_obs_1[, j, i] + upsilon_POST_1[, j], 1)
  }
}
Y_init <- latent_cat
posterior_meansz1 <- matrix(NA, 1, m)
posterior_meansz1 <- t(colMeans(z_pred_outsample_1, dims = 1))
# Plot the real and predicted data
par(mfrow = c(1, 2))
plot(latents[100, ], type = "b", main = "Real data", ylab = "Dependent variable")
plot(posterior_meansz1[, 1], type = "b", main = "Predicted data", ylab = "Latent variable")
# Combine real and predicted data in a single plot
d1 <- data.frame(x = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11), y = latents[100, ])
d2 <- data.frame(x = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11), y = posterior_meansz1[, 1])
par(mar = c(5, 4, 4, 4) + 0.3)
plot(d1$x, d1$y, pch = 16, col = "blue", xlab = "", ylab = "", cex = 2, main = "")
par(new = TRUE)
plot(d2$x, d2$y, pch = 17, col = "green", cex = 2, axes = FALSE, xlab = "Groups", ylab = "Observed variable Y")
axis(side = 4, at = pretty(range(d2$y)))
# Estimates for the classical model
mydata <- read.csv('C:\\Users\\user\\Desktop\\data.csv')
attach(mydata)
library(brms)
library(rstan)
require(brms)
require(rstan)
library(StanHeaders)
# Fit the classical model
model <- brm(formula = Response ~ V1 + V2 + V3 + V4 + V5, data = mydata, family = cumulative)
summary(model)
# Extract parameter estimates and posterior samples
posterior_samples <- posterior_samples(model)
point_estimates <- apply(posterior_samples[, 9:13], 2, mean)
# Calculate MSE and MAE of the classical model
x <- matrix(point_estimates, 5, 1)
y <- posterior_samples[, 9:13]
mse <- matrix(NA, 4000, 5)
for (i in 1:4000) {
  for (j in 1:5) {
    mse[i, j] <- (y[i, j] - x[j, 1])^2
  }
}
mse1 <- mean(mse)
RMSE <- sqrt(mse1)
mae <- matrix(NA, 4000, 5)
for (i in 1:4000) {
  for (j in 1:5) {
    mae[i, j] <- abs(y[i, j] - x[j, 1])
  }
}
mae1 <- mean(mae)
# Calculate coverage probability of the classical model
point_estimates <- matrix(point_estimates, 5, 1)
param_intervals <- posterior_interval(model, prob = 0.95)
param_intervals_2 <- param_intervals[9:13]
coverage_prob <- matrix(NA, 5, 1)
for (j in 1:5) {
  if (point_estimates[j, 1] >= param_intervals_2[j, 1] && point_estimates[j, 1] <= param_intervals_2[j, 2]) {
    coverage_prob[j, 1] <- 1
  } else {
    coverage_prob[j, 1] <- 0
  }
}
mean(coverage_prob)