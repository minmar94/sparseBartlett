library(expm) # For matrix exponential and logarithm operations
library(scoringRules) # For scoring rules such as CRPS
library(ggplot2) # For plotting
library(tidyverse) # For data manipulation and visualization
library(magrittr) # For pipe operations

# Allow script to accept command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Set working directory for simulation outputs
wd <- "./simulations/"

# Load additional required packages for C++ integration and environment management
require(tidyverse)
require(magrittr)
require(BH)
require(Rcpp)
require(withr)

# Set simulation parameters
par1 <- 3
data_max <- 20

# Define a function to compute the Continuous Ranked Probability Score (CRPS)
compute_crps <- function(obs, sim) {
  dd <- as.matrix(dist(sim))
  return(mean(abs(obs - sim)) - 1 / (2 * length(sim)^2) * sum(dd))
}

setwd(wd)

# Example argument structure: c(k, percmiss, sparsityType, perczero, decay)
args <- c(10, 0.1, "Banded1", "0.25", "1.0")

# Loop over different simulation settings for dimension (k)
for (k_sim in c(10, 25)[c(1, 2)]) {
  par2 <- 1
  args[1] <- k_sim

  # Loop over different missing data percentages
  for (percmiss_sim in c("0.1", "0.0")[c(1, 2)]) {
    args[2] <- percmiss_sim

    # Loop over different sparsity patterns
    for (pattern_sim in c("FullRandom", "Banded1", "Banded2", "Banded3")[rev(c(1, 2, 3, 4))]) {
      args[3] <- pattern_sim

      # Set missingness weights depending on pattern
      if (pattern_sim == "FullRandom") {
        w_miss <- rev(c("0.0", "0.5", "0.25", "0.75"))
      } else {
        w_miss <- c("NA")
      }

      # Loop over different proportions of zeros
      for (perczero_sim in w_miss) {
        args[4] <- perczero_sim

        # Find all result files for the current scenario
        fls_replica <- list.files(path = "out_normal", pattern = "LikeNormal", full.names = TRUE)

        # Filter files by current simulation parameters
        fls_replica_k <- fls_replica[grep(paste0("k=", args[1]), fls_replica)]
        name_fold <- paste("k=", args[1], sep = "")
        fls_replica_k <- fls_replica_k[grep(paste0("percmiss=", args[2]), fls_replica_k)]
        name_fold <- paste(name_fold, paste("percmiss=", args[2], sep = ""), sep = " ")
        fls_replica_k <- fls_replica_k[grep(args[3], fls_replica_k)]
        name_fold <- paste(name_fold, paste("pattern=", args[3], sep = ""), sep = " ")
        print(args[3])

        # Further filter for random patterns
        if (args[3] == "Random") {
          fls_replica_k <- fls_replica_k[grep(paste0("perczero=", args[4]), fls_replica_k)]
          name_fold <- paste(name_fold, paste("perczero=", args[4], sep = ""), sep = " ")
          fls_replica_k <- fls_replica_k[grep(paste0("decay=", args[5]), fls_replica_k)]
          name_fold <- paste(name_fold, paste("decay=", args[5], sep = ""), sep = " ")
        }
        if (args[3] == "FullRandom") {
          fls_replica_k <- fls_replica_k[grep(paste0("perczero=", args[4]), fls_replica_k)]
          name_fold <- paste(name_fold, paste("perczero=", args[4], sep = ""), sep = " ")
        }

        # Create output directory for current scenario
        dir.create(file.path("output", name_fold), showWarnings = FALSE)

        # Compile and load the C++ code for G-Wishart sampling
        withr::with_makevars(
          new = c(PKG_LIBS = "$(LAPACK_LIBS) $(BLAS_LIBS)"),
          code = Rcpp::sourceCpp(file = "./gwishart/ggm_new5.cpp")
        )

        # Wrapper function to update the graph structure in the G-Wishart sampler
        update_G <- function(adj, edge_prob, df_0, U, n, seed_values) {
          p <- nrow(adj)
          return(update_G_Rcpp(
            adj = adj,
            edge_prob_mat = matrix(edge_prob, nrow = p, ncol = p),
            df = df_0 + n,
            df_0 = df_0,
            rate = (1 / par2) * diag(p) + U,
            n_edge = p,
            seed = sample.int(n = .Machine$integer.max, size = 1),
            loc_bal = FALSE
          ))
        }

        # Wrapper function to sample from the G-Wishart distribution
        rgwish <- function(adj, df = par1, rate = NULL, seed_values) {
          p <- nrow(adj)
          if (is.null(rate)) rate <- diag(p)
          return(rgwish_Rcpp(
            adj = adj,
            df = df,
            rate = rate,
            seed = sample.int(n = .Machine$integer.max, size = 1)
          ))
        }

        # Main function to run the G-Wishart MCMC for a given dataset
        test <- function(data, n_iter = 1000, percmiss, index_miss_mat = NULL, seed_values = seed_values) {
          n <- nrow(data)
          p <- ncol(data)
          mu <- rep(0, p)
          cdata <- t(t(data) - mu)
          U <- t(cdata) %*% (cdata)

          # Inner function for a single MCMC step
          MCMC_step <- function(adj) {
            update_G(
              adj = adj, edge_prob = 0.5, df_0 = par1, U = U, n = n, seed_values = seed_values
            )
          }

          adj <- matrix(0L, nrow = p, ncol = p)
          for (s in 1:8000) { # Burn-in iterations
            oo <- MCMC_step(adj)
            adj <- oo[[1]]
            K <- oo[[2]]
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

          pb <- txtProgressBar(max = n_iter)
          adjout <- array(0, dim = c(p, p, n_iter + 1))
          Kout <- array(0, dim = c(p, p, n_iter))
          if (percmiss > 0) {
            ysim <- array(NA, dim = c(nrow(index_miss_mat), n_iter))
          } else {
            ysim <- NULL
          }
          adjout[, , 1] <- adj
          # Main MCMC sampling phase
          for (s in 1:n_iter) {
            oo <- MCMC_step(adjout[, , s])
            adjout[, , s + 1] <- oo[[1]]
            Kout[, , s] <- oo[[2]]
            Vp <- solve(n * Kout[, , s] + diag(1 / 1000, p))
            mu <- mvtnorm::rmvnorm(1, mean = c(Vp %*% rowSums(Kout[, , s] %*% t(data))), sigma = Vp)
            cdata <- t(t(data) - c(mu))
            U <- t(cdata) %*% (cdata)
            Sigma <- solve(Kout[, , s])
            # Impute missing values if present
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

        # Run the G-Wishart MCMC for each simulation replicate, handling errors by retrying with a new seed
        gw_replica <- lapply(fls_replica_k, \(i){
          print(i)
          load(i)
          seed_values <- as.numeric(gsub(".*seed=(\\d+).*", "\\1", i))
          seed_values <- 1
          set.seed(seed_values)
          if (as.numeric(args[2]) > 0) {
            success <- FALSE
            while (!success) {
              tryCatch(
                {
                  result <- test(t(y),
                    n_iter = 2000, percmiss = as.numeric(args[2]),
                    index_miss_mat = index_miss_mat, seed_values = seed_values
                  )
                  success <- TRUE
                },
                error = function(e) {
                  cat("Error occurred. Incrementing seed values and retrying...\n")
                  seed_values <- seed_values + 1
                }
              )
            }
            result
          } else {
            success <- FALSE
            while (!success) {
              tryCatch(
                {
                  result <- test(t(y),
                    n_iter = 2000, percmiss = as.numeric(args[2]),
                    seed_values = seed_values
                  )
                  success <- TRUE
                },
                error = function(e) {
                  cat("Error occurred. Incrementing seed values and retrying...\n")
                  seed_values <- seed_values + 1
                  set.seed(seed_values)
                }
              )
            }
            result
          }
        })

        # Functions to compute Kullback-Leibler divergences
        kl_d <- function(lambda, post_samp) {
          return(sum(diag(solve(post_samp) %*% (lambda))) - k - log(det(solve(post_samp) %*% (lambda))))
        }
        kl_inv_d <- function(lambda, post_samp) {
          return(sum(diag((post_samp) %*% solve(lambda))) - k - log(det((post_samp) %*% solve(lambda))))
        }

        # Compute posterior means and medians for Lambda and Z matrices for each replicate
        load(fls_replica_k[1])
        Lambdas_gw <- map(gw_replica, \(l) apply(l[[2]][, , ], c(1, 2), mean))
        Lambdas_gw2 <- Lambdas_gw
        for (m in 1:length(gw_replica)) {
          mat_app <- gw_replica[[m]]$lambda
          num_matrices <- dim(mat_app)[3]
          dim_mat <- dim(mat_app)[1]
          log_array <- array(0, dim = c(dim_mat, dim_mat, num_matrices)) # Initialize log array

          for (iii in 1:num_matrices) {
            log_array[, , iii] <- logm(mat_app[, , iii])
          }
          log_median <- apply(log_array, c(1, 2), mean)
          Lambdas_gw2[[m]] <- expm(log_median)
        }
        Zetas_gw <- map(gw_replica, \(z) apply(z[[1]][, , ], c(1, 2), mean))
        if (as.numeric(args[2]) > 0) Y_gw <- map(gw_replica, \(z) rowMeans(z[[3]][, ]))

        # Initialize lists and matrices for storing results and metrics
        Lambdas_siw2 <- Lambdas_siw <- Zetas_siw <- list()
        KLdiscr_siw <- KLdiscr_gw <- rmse_siw <- rmse_gw <- c()
        KLdiscr_siw_inv <- KLdiscr_gw_inv <- rmse_siw <- rmse_gw <- c()
        BWdiscr_siw <- BWdiscr_gw <- rmse_siw <- rmse_gw <- c()
        crps_sw <- crps_gw <- matrix(NA, ncol = nrow(index_miss_mat), nrow = min(length(fls_replica_k), data_max))
        spec_sw <- spec_gw <- matrix(NA, ncol = 2000, nrow = min(length(fls_replica_k), data_max))
        sens_sw <- sens_gw <- matrix(NA, ncol = 2000, nrow = min(length(fls_replica_k), data_max))
        acc_sw <- acc_gw <- matrix(NA, ncol = 2000, nrow = min(length(fls_replica_k), data_max))
        Y_siw <- matrix(NA, nrow = length(fls_replica_k), ncol = nrow(index_miss_mat))
        zerr_siw <- zerr_gw <- sparse_perc_siw <- sparse_perc_gw <- qsparse_est_siw <- qsparse_est_gw <- c()
        k <- as.numeric(args[1])

        # Loop over each simulation replicate to compute metrics and prepare plots
        for (m in 1:min(length(fls_replica_k), data_max)) {
          load(fls_replica_k[m])
          Lambdas_siw[[m]] <- apply(lambda_out, c(2, 3), mean)
          mat_app <- lambda_out
          num_matrices <- dim(mat_app)[1]
          dim_mat <- dim(mat_app)[2]
          log_array <- array(0, dim = c(dim_mat, dim_mat, num_matrices))
          for (iii in 1:num_matrices) {
            log_array[, , iii] <- logm(mat_app[iii, , ])
          }
          log_median <- apply(log_array, c(1, 2), mean)
          Lambdas_siw2[[m]] <- expm(log_median)
          Zetas_siw[[m]] <- apply(zmat_out, c(2, 3), mean)
          Zetas_gw[[m]][upper.tri(Zetas_gw[[m]])] <- 1
          diag(Zetas_gw[[m]]) <- 1
          n_el <- (k^2 - k) / 2

          # Compute RMSE and CRPS for missing data imputation if applicable
          if (as.numeric(args[2]) > 0) {
            Y_siw[m, ] <- colMeans(missing_out)
            ytrue <- sapply(1:nrow(index_miss_mat), \(i) y[index_miss_mat[i, 1], index_miss_mat[i, 2]])
            rmse_siw[m] <- sqrt(mean((ytrue - Y_siw[m, ])^2))
            rmse_gw[m] <- sqrt(mean((ytrue - Y_gw[[m]])^2))
            missing_gw <- t(gw_replica[[m]]$ymiss)
            crps_sw[m, ] <- sapply(1:nrow(gw_replica[[m]]$ymiss[, ]), \(i) compute_crps(true_y[i], missing_out[, i]))
            crps_gw[m, ] <- sapply(1:nrow(gw_replica[[m]]$ymiss[, ]), \(i) compute_crps(true_y[i], missing_gw[, i]))
          }

          # Compute KL divergence and Bures-Wasserstein distance metrics
          KLdiscr_siw[m] <- sum(diag(solve(Lambdas_siw[[m]]) %*% (lambda))) - k - log(det(solve(Lambdas_siw[[m]]) %*% (lambda)))
          KLdiscr_gw[m] <- sum(diag(solve(Lambdas_gw[[m]]) %*% (lambda))) - k - log(det(solve(Lambdas_gw[[m]]) %*% (lambda)))
          KLdiscr_siw_inv[m] <- sum(diag((Lambdas_siw[[m]]) %*% solve(lambda))) - k - log(det((Lambdas_siw[[m]]) %*% solve(lambda)))
          KLdiscr_gw_inv[m] <- sum(diag((Lambdas_gw[[m]]) %*% solve(lambda))) - k - log(det((Lambdas_gw[[m]]) %*% solve(lambda)))
          app_eig <- eigen(Lambdas_siw[[m]])
          lambda_eig <- app_eig$vectors %*% diag(app_eig$values^0.5) %*% t(app_eig$vectors)
          app_molt_eig <- eigen(lambda_eig %*% lambda %*% lambda_eig)
          molt_lambda_eig <- app_molt_eig$vectors %*% diag(app_molt_eig$values^0.5) %*% t(app_molt_eig$vectors)
          BWdiscr_siw[m] <- sum(diag(Lambdas_siw[[m]] + lambda - 2 * molt_lambda_eig))
          app_eig <- eigen(Lambdas_gw[[m]])
          lambda_eig <- app_eig$vectors %*% diag(app_eig$values^0.5) %*% t(app_eig$vectors)
          app_molt_eig <- eigen(lambda_eig %*% lambda %*% lambda_eig)
          molt_lambda_eig <- app_molt_eig$vectors %*% diag(app_molt_eig$values^0.5) %*% t(app_molt_eig$vectors)
          BWdiscr_gw[m] <- sum(diag(Lambdas_gw[[m]] + lambda - 2 * molt_lambda_eig))

          # Compute errors and sparsity metrics for Z matrices
          zerr_siw[m] <- sum(abs(z_mat[lower.tri(z_mat)] - Zetas_siw[[m]][lower.tri(Zetas_siw[[m]])]))
          zerr_gw[m] <- sum(abs(z_mat[lower.tri(z_mat)] - Zetas_gw[[m]][lower.tri(Zetas_gw[[m]])]))
          sparse_perc_siw[m] <- sum(round(Zetas_siw[[m]]) == 0) / sum(z_mat == 0)
          sparse_perc_gw[m] <- sum(round(Zetas_gw[[m]]) == 0) / sum(z_mat == 0)
          qsparse_est_siw[m] <- length(intersect(which(z_mat == 0), which(round(Zetas_siw[[m]]) == 0))) / sum(z_mat == 0)
          qsparse_est_gw[m] <- length(intersect(which(z_mat == 0), which(round(Zetas_gw[[m]]) == 0))) / sum(z_mat == 0)

          n_el <- (k^2 - k) / 2

          # Compute accuracy, sensitivity, and specificity for edge recovery
          for (i in 1:dim(zmat_out)[1]) {
            tt <- table(factor(z_mat[lower.tri(z_mat[, ])], levels = c(0, 1)), factor(zmat_out[i, , ][lower.tri(zmat_out[i, , ])], levels = c(0, 1)))
            acc_sw[m, i] <- (tt[1, 1] + tt[2, 2]) / sum(tt)
            sens_sw[m, i] <- tt[1, 1] / (tt[1, 1] + tt[1, 2])
            spec_sw[m, i] <- (tt[2, 2]) / (tt[2, 2] + tt[2, 1])


            tt <- table(z_mat[lower.tri(z_mat[, ])], gw_replica[[m]]$zmat[, , -c(1)][, , i][lower.tri(zmat_out[i, , ])])
            if (perczero_sim == "0.0") {
              tt <- rbind(c(0, 0), tt)
            }
            acc_gw[m, i] <- (tt[1, 1] + tt[2, 2]) / sum(tt)
            sens_gw[m, i] <- tt[1, 1] / (tt[1, 1] + tt[1, 2])
            spec_gw[m, i] <- (tt[2, 2]) / (tt[2, 2] + tt[2, 1])
          }

          # Prepare heatmap plots for Lambda matrices (true, SIW, GW, and their differences)
          plist <- list(
            lambda, Lambdas_siw[[m]], Lambdas_gw[[m]],
            Lambdas_siw[[m]] - lambda, Lambdas_gw[[m]] - lambda, Lambdas_siw[[m]] - Lambdas_gw[[m]]
          ) %>%
            map(\(lambda){
              lambda %>%
                as.data.frame() %>%
                rownames_to_column(var = "rn") %>%
                pivot_longer(-rn, names_to = "cn", values_to = "val") %>%
                mutate(
                  rn = factor(rn, levels = rev(1:nrow(lambda))),
                  cn = factor(cn, levels = paste0("V", 1:nrow(lambda)))
                ) %>%
                ggplot(aes(cn, rn, fill = val)) +
                geom_raster() +
                scale_fill_distiller(palette = "RdBu")
            })
        }

        # Save all results for the current scenario
        save.image(file = paste0(getwd(), "/for_figures_normal/LikeNormal_", paste0(args, collapse = "_"), ".RData"))
      }
    }
  }
}
