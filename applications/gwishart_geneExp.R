# ========
# Genetic 50 Miss
# ========

# Load MCMC output for the 50-gene experiment with missing data
# setwd("/Users/gianlucamastrantonio/Library/CloudStorage/GoogleDrive-mastrantonio.gluca@gmail.com/Shared drives/Size and Shape/dataset/results")
load("./applications/out/GeneExp50_missmcmc_out_lambda.Rdata")

# Load required libraries for data manipulation and plotting
library(ggplot2)
library(tidyverse)
library(magrittr)

set.seed(1)

# ========
# G-Wishart
# ========

# Set G-Wishart prior parameters
par1 <- 3
par2 <- 1

# Compile and load the C++ code for G-Wishart sampling using Rcpp
withr::with_makevars(
  new = c(PKG_LIBS = "$(LAPACK_LIBS) $(BLAS_LIBS)"),
  code = Rcpp::sourceCpp(file = "./applications/gwishart/ggm_new5.cpp")
)

# Wrapper function for updating the adjacency matrix (G) using the C++ function
# adj: adjacency matrix
# edge_prob: prior edge probability
# df_0: prior degrees of freedom
# U: data cross-product matrix
# n: sample size
update_G <- function(adj, edge_prob, df_0, U, n) {
  p <- nrow(adj)
  return(update_G_Rcpp(
    adj = adj,
    edge_prob_mat = matrix(edge_prob, nrow = p, ncol = p),
    df = df_0 + n,
    df_0 = df_0,
    rate = (1 / par2) * diag(p) + U,
    # rate = (1.0/(2.0*10.0^(-1)))*diag(p) + U,
    n_edge = p, # Number of single edge updates attempted in this MCMC update
    seed = sample.int(n = .Machine$integer.max, size = 1),
    loc_bal = FALSE
  ))
}

# Wrapper function for sampling from the G-Wishart distribution using the C++ function
# adj: adjacency matrix
# df: degrees of freedom
# rate: rate matrix
rgwish <- function(adj, df = par1, rate = NULL) {
  p <- nrow(adj)
  if (is.null(rate)) rate <- diag(p)
  return(rgwish_Rcpp(
    adj = adj,
    df = df,
    rate = rate,
    seed = sample.int(n = .Machine$integer.max, size = 1)
  ))
}

# Main function for running the MCMC model
# data: observed data matrix
# n_iter: number of MCMC iterations
# percmiss: percentage of missing data
# index_miss_mat: matrix of missing value indices
test <- function(data, n_iter = 1000, percmiss, index_miss_mat = NULL) {
  n <- nrow(data)
  p <- ncol(data)
  # Estimate means (mu) instead of assuming zero mean
  mu <- rep(0, p)
  cdata <- t(t(data) - mu)
  U <- t(cdata) %*% (cdata)

  # Define a single MCMC step for updating the graph and precision matrix
  MCMC_step <- function(adj) {
    update_G(
      adj = adj, edge_prob = 0.5, df_0 = par1, U = U, n = n # parameters can be passed as input
    )
  }

  # Initialize adjacency matrix
  adj <- matrix(0L, nrow = p, ncol = p)
  # Burn-in phase: update adjacency and precision matrices
  for (s in 1:6000) { # Burn-in iterations
    oo <- MCMC_step(adj)
    adj <- oo[[1]]
    K <- oo[[2]]
    # Update mean and covariance for the data
    Vp <- solve(n * K + diag(1 / 1000, p))
    mu <- mvtnorm::rmvnorm(1, mean = c(Vp %*% rowSums(K %*% t(data))), sigma = Vp)
    cdata <- t(t(data) - c(mu))
    U <- t(cdata) %*% (cdata)
    Sigma <- solve(K)
    # Impute missing values if present
    if (percmiss > 0) {
      for (imiss in 1:nrow(index_miss_mat)) {
        icur <- index_miss_mat[imiss, ]
        varcond <- Sigma[icur[1], icur[1]] - Sigma[icur[1], -icur[1]] %*% solve(Sigma[-icur[1], -icur[1]]) %*% Sigma[icur[1], -icur[1]]
        mucond <- mu[icur[1]] + Sigma[icur[1], -icur[1]] %*% solve(Sigma[-icur[1], -icur[1]]) %*% (t(data)[-icur[1], icur[2]] - mu[-icur[1]])
        data[icur[2], icur[1]] <- rnorm(1, mucond, sqrt(varcond))
      }
    }
  }

  # Set up progress bar for MCMC sampling
  pb <- txtProgressBar(max = n_iter)

  # Arrays to store sampled adjacency and precision matrices
  adjout <- array(0, dim = c(p, p, n_iter + 1))
  Kout <- array(0, dim = c(p, p, n_iter))
  if (percmiss > 0) {
    ysim <- array(NA, dim = c(nrow(index_miss_mat), n_iter))
  } else {
    ysim <- NULL
  }
  adjout[, , 1] <- adj

  # Main MCMC loop
  for (s in 1:n_iter) {
    oo <- MCMC_step(adjout[, , s])
    adjout[, , s + 1] <- oo[[1]]
    Kout[, , s] <- oo[[2]]
    Vp <- solve(n * Kout[, , s] + diag(1 / 1000, p))
    mu <- mvtnorm::rmvnorm(1, mean = c(Vp %*% rowSums(Kout[, , s] %*% t(data))), sigma = Vp)
    cdata <- t(t(data) - c(mu))
    U <- t(cdata) %*% (cdata)
    Sigma <- solve(Kout[, , s])
    # Impute missing values if present and store them
    if (percmiss > 0) {
      for (imiss in 1:nrow(index_miss_mat)) {
        icur <- index_miss_mat[imiss, ]
        varcond <- Sigma[icur[1], icur[1]] - Sigma[icur[1], -icur[1]] %*% solve(Sigma[-icur[1], -icur[1]]) %*% Sigma[icur[1], -icur[1]]
        mucond <- mu[icur[1]] + Sigma[icur[1], -icur[1]] %*% solve(Sigma[-icur[1], -icur[1]]) %*% (t(data)[-icur[1], icur[2]] - mu[-icur[1]])
        data[icur[2], icur[1]] <- rnorm(1, mucond, sqrt(varcond))
        ysim[imiss, s] <- data[icur[2], icur[1]]
      }
    }
    # Update progress bar
    setTxtProgressBar(pb, s)
  }

  close(pb)
  # Return sampled adjacency matrices, precision matrices, and imputed values
  return(list(zmat = adjout, lambda = Kout, ymiss = ysim))
}

# Run the MCMC model for the gene expression data
mcmc_gwish <- test(t(y), n_iter = 4000, percmiss = 10, index_miss_mat = index_miss_mat)

# Compute mean imputed values for comparison
miss_sw <- colMeans(missing_out)
miss_gw <- rowMeans(mcmc_gwish$ymiss[, -c(1:2000)])

# Plot comparison between two imputation methods
plot(miss_sw, miss_gw)
abline(a = 0, b = 1, col = 2)

# Compute RMSE for both imputation methods
sqrt(mean((miss_sw - true_y)^2))
sqrt(mean((miss_gw - true_y)^2))

# Prepare arrays for CRPS calculation
missing_gw <- t(mcmc_gwish$ymiss[, -c(1:2000)])
crps_sw <- c()
crps_gw <- c()

# Compute CRPS for each missing value
for (i in 1:ncol(missing_out))
{
  dd <- as.matrix(missing_out[, i])
  crps_sw[i] <- mean(abs(true_y[i] - missing_out[, i])) - 1 / ncol(missing_out)^2 * sum(dd)
  dd <- as.matrix(missing_gw[, i])
  crps_gw[i] <- mean(abs(true_y[i] - missing_gw[, i])) - 1 / ncol(missing_gw)^2 * sum(dd)
}
plot(crps_sw, crps_gw)
abline(a = 0, b = 1, col = 2)

# Summarize CRPS results
mean(crps_sw)
mean(crps_gw)

# Inspect first sampled adjacency and precision matrices
mcmc_gwish$zmat[, , 1]
round(mcmc_gwish$lambda[1:5, 1:5, 1], 4)
mcmc_gwish$zmat[1:5, 1:5, 1]

# ========
# R
# ========

# Calculate the number of edges in each sampled graph for both methods
nsim <- dim(zmat_out)[1]
k <- dim(zmat_out)[2]
n_edge_sw <- rep(NA, nsim)
n_edge_gw <- rep(NA, nsim)
n_el <- (k^2 - k) / 2

for (i in 1:nsim)
{
  n_edge_sw[i] <- n_el - sum(c(zmat_out[i, , ] == 0))
  n_edge_gw[i] <- n_el - sum(c(mcmc_gwish$zmat[, , i][lower.tri(mcmc_gwish$zmat[, , i])] == 0))
}

# Prepare data for plotting edge counts
data_plot <- data.frame(n_edge = c(n_edge_sw, n_edge_gw), type = rep(c("SW", "GW"), each = nsim))
ggplot(data_plot, aes(y = n_edge))
# par(mfrow = c(1, 2))
# plot(n_edge_sw, type = "l")
# plot(n_edge_gw, type = "l")
# par(mfrow = c(1, 1))

# Save workspace image for later analysis
save.image("./applications/out/GeneExp50_Confronti_GW.Rdata")

# ========
#
# ========

# Compute mean adjacency matrix (posterior edge probabilities)
zmap <- apply(zmat_out, c(2, 3), mean)
zmap[upper.tri(zmap)] <- t(zmap[lower.tri(zmap)])

library(tidyverse)
library(magrittr)

# Plot heatmap of posterior edge probabilities (commented out)
# pdf("/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/Applicazioni/Overleaf/SparsePrecisionPatterns/gene50pattern.pdf")
# zmap %>%
#  as.data.frame() %>%
#  set_colnames(value = 1:k) %>%
#  rownames_to_column(var = "V1") %>%
#  pivot_longer(-V1, names_to = "V2", values_to = "val") %>%
#  mutate(V1 = as.numeric(V1), V2 = as.numeric(V2)) %>%
#  ggplot(aes(x = V1, y = rev(V2), fill = round(val))) +
#  geom_tile() +
#  scale_fill_distiller(palette = "YlOrRd")
# dev.off()

# ========
# Repeat for 100-gene experiment
# ========

# Load MCMC output for the 100-gene experiment with missing data
load("./applications/out/GeneExp100_missmcmc_out_lambda.Rdata")

# Set G-Wishart prior parameters
par1 <- 3
par2 <- 1

# Compile and load the C++ code for G-Wishart sampling using Rcpp
withr::with_makevars(
  new = c(PKG_LIBS = "$(LAPACK_LIBS) $(BLAS_LIBS)"),
  code = Rcpp::sourceCpp(file = "./applications/gwishart/ggm_new5.cpp")
)

# Wrapper function for updating the adjacency matrix (G) using the C++ function
update_G <- function(adj, edge_prob, df_0, U, n) {
  p <- nrow(adj)
  return(update_G_Rcpp(
    adj = adj,
    edge_prob_mat = matrix(edge_prob, nrow = p, ncol = p),
    df = df_0 + n,
    df_0 = df_0,
    rate = (1 / par2) * diag(p) + U,
    # rate = (1.0/(2.0*10.0^(-1)))*diag(p) + U,
    n_edge = p, # Number of single edge updates attempted in this MCMC update
    seed = sample.int(n = .Machine$integer.max, size = 1),
    loc_bal = FALSE
  ))
}

# Wrapper function for sampling from the G-Wishart distribution using the C++ function
rgwish <- function(adj, df = par1, rate = NULL) {
  p <- nrow(adj)
  if (is.null(rate)) rate <- diag(p)
  return(rgwish_Rcpp(
    adj = adj,
    df = df,
    rate = rate,
    seed = sample.int(n = .Machine$integer.max, size = 1)
  ))
}

# Main function for running the MCMC model (same as above)
test <- function(data, n_iter = 1000, percmiss, index_miss_mat = NULL) {
  n <- nrow(data)
  p <- ncol(data)
  mu <- rep(0, p)
  cdata <- t(t(data) - mu)
  U <- t(cdata) %*% (cdata)

  MCMC_step <- function(adj) {
    update_G(
      adj = adj, edge_prob = 0.5, df_0 = par1, U = U, n = n
    )
  }

  adj <- matrix(0L, nrow = p, ncol = p)
  for (s in 1:6000) { # Burn-in iterations
    oo <- MCMC_step(adj)
    adj <- oo[[1]]
    K <- oo[[2]]
    Vp <- solve(n * K + diag(1 / 1000, p))
    mu <- mvtnorm::rmvnorm(1, mean = c(Vp %*% rowSums(K %*% t(data))), sigma = Vp)
    cdata <- t(t(data) - c(mu))
    U <- t(cdata) %*% (cdata)
    Sigma <- solve(K)
    if (percmiss > 0) {
      for (imiss in 1:nrow(index_miss_mat)) {
        icur <- index_miss_mat[imiss, ]
        varcond <- Sigma[icur[1], icur[1]] - Sigma[icur[1], -icur[1]] %*% solve(Sigma[-icur[1], -icur[1]]) %*% Sigma[icur[1], -icur[1]]
        mucond <- mu[icur[1]] + Sigma[icur[1], -icur[1]] %*% solve(Sigma[-icur[1], -icur[1]]) %*% (t(data)[-icur[1], icur[2]] - mu[-icur[1]])
        data[icur[2], icur[1]] <- rnorm(1, mucond, sqrt(varcond))
      }
    }
  }

  pb <- txtProgressBar(max = n_iter)
  adjout <- array(0, dim = c(p, p, n_iter + 1))
  Kout <- array(0, dim = c(p, p, n_iter))
  if (percmiss > 0) {
    ysim <- array(NA, dim = c(nrow(index_miss_mat), n_iter))
  } else {
    ysim <- NULL
  }
  adjout[, , 1] <- adj
  for (s in 1:n_iter) {
    oo <- MCMC_step(adjout[, , s])
    adjout[, , s + 1] <- oo[[1]]
    Kout[, , s] <- oo[[2]]
    Vp <- solve(n * Kout[, , s] + diag(1 / 1000, p))
    mu <- mvtnorm::rmvnorm(1, mean = c(Vp %*% rowSums(Kout[, , s] %*% t(data))), sigma = Vp)
    cdata <- t(t(data) - c(mu))
    U <- t(cdata) %*% (cdata)
    Sigma <- solve(Kout[, , s])
    if (percmiss > 0) {
      for (imiss in 1:nrow(index_miss_mat)) {
        icur <- index_miss_mat[imiss, ]
        varcond <- Sigma[icur[1], icur[1]] - Sigma[icur[1], -icur[1]] %*% solve(Sigma[-icur[1], -icur[1]]) %*% Sigma[icur[1], -icur[1]]
        mucond <- mu[icur[1]] + Sigma[icur[1], -icur[1]] %*% solve(Sigma[-icur[1], -icur[1]]) %*% (t(data)[-icur[1], icur[2]] - mu[-icur[1]])
        data[icur[2], icur[1]] <- rnorm(1, mucond, sqrt(varcond))
        ysim[imiss, s] <- data[icur[2], icur[1]]
      }
    }
    setTxtProgressBar(pb, s)
  }

  close(pb)
  return(list(zmat = adjout, lambda = Kout, ymiss = ysim))
}

# Run the MCMC model for the 100-gene data
mcmc_gwish <- test(t(y), n_iter = 4000, percmiss = 10, index_miss_mat = index_miss_mat)

# Compute mean imputed values for comparison
miss_sw <- colMeans(missing_out)
miss_gw <- rowMeans(mcmc_gwish$ymiss[, -c(1:2000)])

# Compute RMSE for both imputation methods
sqrt(mean((miss_sw - true_y)^2))
sqrt(mean((miss_gw - true_y)^2))

# Prepare arrays for CRPS calculation
missing_gw <- t(mcmc_gwish$ymiss[, -c(1:2000)])
crps_sw <- c()
crps_gw <- c()

# Compute CRPS for each missing value
for (i in 1:ncol(missing_out))
{
  dd <- as.matrix(missing_out[, i])
  crps_sw[i] <- mean(abs(true_y[i] - missing_out[, i])) - 1 / ncol(missing_out)^2 * sum(dd)
  dd <- as.matrix(missing_gw[, i])
  crps_gw[i] <- mean(abs(true_y[i] - missing_gw[, i])) - 1 / ncol(missing_gw)^2 * sum(dd)
}
# plot(crps_sw, crps_gw)
# abline(a = 0, b = 1, col = 2)

# Summarize CRPS results
summary(crps_sw)
summary(crps_gw)

# Inspect first sampled adjacency and precision matrices
mcmc_gwish$zmat[, , 1]
round(mcmc_gwish$lambda[1:5, 1:5, 1], 4)
mcmc_gwish$zmat[1:5, 1:5, 1]

# ========
# R
# ========

# Calculate the number of edges in each sampled graph for both methods
nsim <- dim(zmat_out)[1]
k <- dim(zmat_out)[2]
n_edge_sw <- rep(NA, nsim)
n_edge_gw <- rep(NA, nsim)
n_el <- (k^2 - k) / 2

for (i in 1:nsim)
{
  n_edge_sw[i] <- n_el - sum(c(zmat_out[i, , ] == 0))
  n_edge_gw[i] <- n_el - sum(c(mcmc_gwish$zmat[, , i][lower.tri(mcmc_gwish$zmat[, , i])] == 0))
}



save.image("./applications/out/GeneExp100_Confronti_GW.Rdata")
