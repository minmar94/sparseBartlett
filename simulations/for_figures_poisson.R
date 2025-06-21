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

# Loop over different model types (here only "LikePoisson")
for (patt in c("LikePoisson")) {
  # Loop over different values of k (number of variables)
  for (k_sim in c(10, 25)[c(2, 1)]) {
    par2 <- 1
    args[1] <- k_sim
    # Loop over different missing data percentages
    for (percmiss_sim in c("0.1", "0.0")[c(2, 1)]) {
      args[2] <- percmiss_sim
      # Loop over different sparsity patterns
      for (pattern_sim in c("FullRandom", "Banded1", "Banded2", "Banded3")[c(1, 2, 3, 4)]) {
        args[3] <- pattern_sim
        # Set possible proportions of zeros depending on the pattern
        if (pattern_sim == "FullRandom") {
          w_miss <- c("0.0", "0.5", "0.25", "0.75")
        } else {
          w_miss <- c("NA")
        }
        # Loop over different proportions of zeros
        for (perczero_sim in w_miss) {
          args[4] <- perczero_sim
          # List all result files for the current scenario
          fls_replica <- list.files(path = "out_poisson", pattern = patt, full.names = T)
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
          print(fls_replica_k)
          # Load the first result file to initialize variables
          load(fls_replica_k[1])
          Lambdas_siw <- Zetas_siw <- list()
          KLdiscr_siw <- rmse_siw <- c()
          KLdiscr_siw_inv <- rmse_siw <- c()
          BWdiscr_siw <- rmse_siw <- c()
          crps_sw <- matrix(NA, ncol = nrow(index_miss_mat), nrow = min(length(fls_replica_k), data_max))
          spec_sw <- matrix(NA, ncol = dim(zmat_out)[1], nrow = min(length(fls_replica_k), data_max))
          sens_sw <- matrix(NA, ncol = dim(zmat_out)[1], nrow = min(length(fls_replica_k), data_max))
          acc_sw <- matrix(NA, ncol = dim(zmat_out)[1], nrow = min(length(fls_replica_k), data_max))
          Y_siw <- matrix(NA, nrow = length(fls_replica_k), ncol = nrow(index_miss_mat))
          zerr_siw <- sparse_perc_siw <- qsparse_est_siw <- c()
          k <- as.numeric(args[1])
          # Loop over each simulation replicate
          for (m in 1:min(length(fls_replica_k), data_max)) {
            load(fls_replica_k[m])
            # Compute posterior mean for Lambda and Z matrices
            Lambdas_siw[[m]] <- apply(lambda_out, c(2, 3), mean)
            Zetas_siw[[m]] <- apply(zmat_out, c(2, 3), mean)
            n_el <- (k^2 - k) / 2
            # Compute RMSE and CRPS for missing data imputation if applicable
            if (as.numeric(args[2]) > 0) {
              Y_siw[m, ] <- colMeans(missing_out)
              ytrue <- sapply(1:nrow(index_miss_mat), \(i) y[index_miss_mat[i, 1], index_miss_mat[i, 2]])
              rmse_siw[m] <- sqrt(mean((ytrue - Y_siw[m, ])^2))
              crps_sw[m, ] <- sapply(1:ncol((missing_out)), \(i) compute_crps(true_y[i], missing_out[, i]))
            }
            # Compute KL divergence and Bures-Wasserstein distance metrics
            KLdiscr_siw[m] <- sum(diag(solve(Lambdas_siw[[m]]) %*% (lambda))) - k - log(det(solve(Lambdas_siw[[m]]) %*% (lambda)))
            KLdiscr_siw_inv[m] <- sum(diag((Lambdas_siw[[m]]) %*% solve(lambda))) - k - log(det((Lambdas_siw[[m]]) %*% solve(lambda)))
            app_eig <- eigen(Lambdas_siw[[m]])
            lambda_eig <- app_eig$vectors %*% diag(app_eig$values^0.5) %*% t(app_eig$vectors)
            app_molt_eig <- eigen(lambda_eig %*% lambda %*% lambda_eig)
            molt_lambda_eig <- app_molt_eig$vectors %*% diag(app_molt_eig$values^0.5) %*% t(app_molt_eig$vectors)
            BWdiscr_siw[m] <- sum(diag(Lambdas_siw[[m]] + lambda - 2 * molt_lambda_eig))
            # Compare Z: error and sparsity metrics
            zerr_siw[m] <- sum(abs(z_mat[lower.tri(z_mat)] - Zetas_siw[[m]][lower.tri(Zetas_siw[[m]])]))
            sparse_perc_siw[m] <- sum(round(Zetas_siw[[m]]) == 0) / sum(z_mat == 0)
            qsparse_est_siw[m] <- length(intersect(which(z_mat == 0), which(round(Zetas_siw[[m]]) == 0))) / sum(z_mat == 0)
            n_el <- (k^2 - k) / 2
            # Compute accuracy, sensitivity, and specificity for edge recovery
            for (i in 1:dim(zmat_out)[1]) {
              tt <- table(factor(z_mat[lower.tri(z_mat[, ])], levels = c(0, 1)), factor(zmat_out[i, , ][lower.tri(zmat_out[i, , ])], levels = c(0, 1)))
              acc_sw[m, i] <- (tt[1, 1] + tt[2, 2]) / sum(tt)
              sens_sw[m, i] <- tt[1, 1] / (tt[1, 1] + tt[1, 2])
              spec_sw[m, i] <- (tt[2, 2]) / (tt[2, 2] + tt[2, 1])
            }
            # Prepare heatmap plots for Lambda matrices (true, SIW, and their difference)
            plist <- list(
              lambda, Lambdas_siw[[m]],
              Lambdas_siw[[m]] - lambda
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
          save.image(file = paste0(getwd(), "/for_figures_poisson/", patt, "_", name_fold, ".rdata"))
        }
      }
    }
  }
}
# End
