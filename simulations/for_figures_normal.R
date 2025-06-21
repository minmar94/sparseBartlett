library(expm)
library(scoringRules)
library(ggplot2)
library(tidyverse)
library(magrittr)
# !/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

wd <- "./simulations/"
# Carico i pacchetti
# Wishart(1.0*(k+3-1), (1.0/(2.0*10.0^(-5))).*Matrix(I(k))),
require(tidyverse)
require(magrittr)
require(BH)
require(Rcpp)
require(withr)

# 1.0*(k+3-1)
# (1.0/(k+3-1)).
par1 <- 3
data_max <- 20

### crps function
compute_crps <- function(obs, sim) {
  dd <- as.matrix(dist(sim))
  return(mean(abs(obs - sim)) - 1 / (2 * length(sim)^2) * sum(dd))
}
setwd(wd)


# args = c(k, percmiss, tipoSparsità, perczero, decay)# Sono stringhe, per cui ci vanno i valori esatti con cui hai simulato!
args <- c(10, 0.1, "Banded1", "0.25", "1.0")
for (k_sim in c(10, 25)[c(1, 2)]) {
  par2 <- 1

  args[1] <- k_sim
  # for(percmiss_sim in c("0.0","0.1")){
  for (percmiss_sim in c("0.1", "0.0")[c(1, 2)])
  {
    args[2] <- percmiss_sim

    for (pattern_sim in c("FullRandom", "Banded1", "Banded2", "Banded3")[rev(c(1, 2, 3, 4))])
    {
      args[3] <- pattern_sim
      if (pattern_sim == "FullRandom") {
        w_miss <- rev(c("0.0", "0.5", "0.25", "0.75"))
      } else {
        w_miss <- c("NA")
      }
      for (perczero_sim in w_miss)
      {
        args[4] <- perczero_sim
        # for (decay_sim in c("6.0", "3.0", "1.0"))
        # {
        #  args[5] <- decay_sim
        #  print(args)
        # }

        # if (pattern_sim == "Random") {
        #  for (perczero_sim in c("0.0", "0.5", "0.25", "0.75"))
        #  {
        #    args[4] <- perczero_sim
        #    for (decay_sim in c("6.0", "3.0", "1.0"))
        #    {
        #      args[5] <- decay_sim
        #      print(args)
        #    }
        #  }
        # }
        # args <- c("10",     "0.1"   ,   "Random" ,"0.0"   ,   "6.0"  )


        # PER ESEGUIRE DA TERMINALE #
        # Rscript --vanilla confronto_gwishart.R 10 0.1 Banded1 0.25 1.0


        # Carico WS #
        # prendo dalla cartella Simulations che contiene le workspace, solo le ws relative alla Normale


        fls_replica <- list.files(path = "out_normal", pattern = "LikeNormal", full.names = T)

        # Qui mi prendo solo quelle per k = 10
        fls_replica_k <- fls_replica[grep(paste0("k=", args[1]), fls_replica)]
        name_fold <- paste("k=", args[1], sep = "")
        # Qui sottoseleziono solo quelle con percmiss = 0
        fls_replica_k <- fls_replica_k[grep(paste0("percmiss=", args[2]), fls_replica_k)]
        name_fold <- paste(name_fold, paste("percmiss=", args[2], sep = ""), sep = " ")
        # Qui solo quelle con pattern di sparsità random
        fls_replica_k <- fls_replica_k[grep(args[3], fls_replica_k)]
        name_fold <- paste(name_fold, paste("pattern=", args[3], sep = ""), sep = " ")
        print(args[3])
        if (args[3] == "Random") {
          fls_replica_k <- fls_replica_k[grep(paste0("perczero=", args[4]), fls_replica_k)]
          name_fold <- paste(name_fold, paste("perczero=", args[4], sep = ""), sep = " ")
          fls_replica_k <- fls_replica_k[grep(paste0("decay=", args[5]), fls_replica_k)]
          name_fold <- paste(name_fold, paste("decay=", args[5], sep = ""), sep = " ")
        }
        if (args[3] == "FullRandom") {
          fls_replica_k <- fls_replica_k[grep(paste0("perczero=", args[4]), fls_replica_k)]
          name_fold <- paste(name_fold, paste("perczero=", args[4], sep = ""), sep = " ")
          # fls_replica_k <- fls_replica_k[grep(paste0("decay=", args[5]), fls_replica_k)]
          # name_fold <- paste(name_fold, paste("decay=", args[5], sep = ""), sep = " ")
        }
        dir.create(file.path("output", name_fold), showWarnings = FALSE)
        # Potrei andare avanti a piacere. La logica sta nel selezionare in qualche modo le classi di modelli per cui è ragionevole confrontare il nostro approccio con la G-Wishart. Una volta finiti tutti gli esperimenti, ogni nostro singolo scenario di simulazione dovrebbe contenere 50 repliche (e.g. 50 ws).
        print(fls_replica_k)

        #### Definisco cose per la GWishart rubate e modificare a Van de Boom (magari esplode visto che si chiama così) ####
        # setwd("/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/lavori/SSsim")
        withr::with_makevars(
          new = c(PKG_LIBS = "$(LAPACK_LIBS) $(BLAS_LIBS)"),
          code = Rcpp::sourceCpp(file = "./gwishart/ggm_new5.cpp")
        )
        # withr::with_makevars(
        #  new = c(PKG_LIBS = "-llapack -lblas"),
        #  code = Rcpp::sourceCpp(file = "ggm_new5.cpp")
        # )
        # setwd("/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/lavori/")


        # Update G wishart
        # Questo è un wrapper per la funzione di Rcpp
        # Da capire come tunare df_0 e edge_prob
        update_G <- function(adj, edge_prob, df_0, U, n, seed_values) {
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

        # Wrapper per la generazione della G wishart
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

        # Main function del modello
        test <- function(data, n_iter = 1000, percmiss, index_miss_mat = NULL, seed_values = seed_values) {
          n <- nrow(data)
          p <- ncol(data)
          # NOSTRA MODIFICA: gli facciamo stimare le medie, lui di base stima a media 0 #
          mu <- rep(0, p)
          cdata <- t(t(data) - mu)
          U <- t(cdata) %*% (cdata)

          # MCMC
          MCMC_step <- function(adj) {
            update_G(
              adj = adj, edge_prob = 0.5, df_0 = par1, U = U, n = n, seed_values = seed_values # parametri che si potrebbero passare in input
            )
          }

          adj <- matrix(0L, nrow = p, ncol = p)
          for (s in 1:8000) { # Burn-in iterations
            oo <- MCMC_step(adj)
            adj <- oo[[1]]
            K <- oo[[2]]
            # Aggiorniamo la matrice di varianza/covarianza per la media (controllare che la nostra prior su mu sia coerente)
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

          # Come sopra, ma MCMC a convergenza
          # Salvo le matrici in un array di dimensione pari a k x k x num iter.
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

            # adj_cum <- adj_cum + adj
            setTxtProgressBar(pb, s)
          }

          close(pb)
          return(list(zmat = adjout, lambda = Kout, ymiss = ysim))
        }



        # Confronto ---------------------------------------------------------------
        # Per ciascuna replica dello scenario considerato, stimo la G wishart
        # Finche k è piccolo e le ws non sono assai, si può girare anche in locale (siamo nell'ordine dei minuti)
        # Ogni tanto dice che rgwish non converge, ma non mi sembra comporti problemi di stima o altre criticità perché continua a funzionare tutto e pure abbastanza bene
        gw_replica <- lapply(fls_replica_k, \(i){
          print(i)
          load(i)

          # Extract seed values
          seed_values <- as.numeric(gsub(".*seed=(\\d+).*", "\\1", i))
          seed_values <- 1
          set.seed(seed_values)
          if (as.numeric(args[2]) > 0) {
            # test(t(y), n_iter = 2000, percmiss = as.numeric(args[2]), index_miss_mat = index_miss_mat, seed_values = seed_values)



            success <- FALSE # Flag to track whether the function call succeeded
            while (!success) {
              tryCatch(
                {
                  # Attempt to run the function
                  result <- test(t(y),
                    n_iter = 2000, percmiss = as.numeric(args[2]),
                    index_miss_mat = index_miss_mat, seed_values = seed_values
                  )

                  # If no error occurs, mark as success
                  success <- TRUE
                },
                error = function(e) {
                  # If an error occurs, increment the seed values and try again
                  cat("Error occurred. Incrementing seed values and retrying...\n")
                  seed_values <- seed_values + 1
                }
              )
            }
            result
          } else {
            # test(t(y), n_iter = 2000, percmiss = as.numeric(args[2]), seed_values = seed_values)

            success <- FALSE # Flag to track whether the function call succeeded
            while (!success) {
              tryCatch(
                {
                  # Attempt to run the function
                  result <- test(t(y),
                    n_iter = 2000, percmiss = as.numeric(args[2]),
                    seed_values = seed_values
                  )

                  # If no error occurs, mark as success
                  success <- TRUE
                },
                error = function(e) {
                  # If an error occurs, increment the seed values and try again
                  cat("Error occurred. Incrementing seed values and retrying...\n")
                  seed_values <- seed_values + 1
                  set.seed(seed_values)
                }
              )
            }
            result
          }
        })
        ## QUI
        kl_d <- function(lambda, post_samp) {
          return(sum(diag(solve(post_samp) %*% (lambda))) - k - log(det(solve(post_samp) %*% (lambda))))
        }
        kl_inv_d <- function(lambda, post_samp) {
          return(sum(diag((post_samp) %*% solve(lambda))) - k - log(det((post_samp) %*% solve(lambda))))
        }



        # Salvo la media a posteriori di Lambda per ciascuna replica. Quindi questi oggetti sono liste contenenti la media stimata di lambda e la media di Z relative a ciascuna replica
        load(fls_replica_k[1])

        Lambdas_gw <- map(gw_replica, \(l) apply(l[[2]][, , ], c(1, 2), mean))
        Lambdas_gw2 <- Lambdas_gw
        # KL Diverìgence
        for (m in 1:length(gw_replica))
        {
          mat_app <- gw_replica[[m]]$lambda
          num_matrices <- dim(mat_app)[3]
          dim_mat <- dim(mat_app)[1]
          log_array <- array(0, dim = c(dim_mat, dim_mat, num_matrices)) # Initialize log array

          for (iii in 1:num_matrices) {
            log_array[, , iii] <- logm(mat_app[, , iii]) # Apply matrix logarithm to each SPD matrix
          }
          # Compute the element-wise median across all matrices
          log_median <- apply(log_array, c(1, 2), mean)

          # Compute matrix exponential to map back to SPD space
          Lambdas_gw2[[m]] <- expm(log_median)
        }




        Zetas_gw <- map(gw_replica, \(z) apply(z[[1]][, , ], c(1, 2), mean))
        if (as.numeric(args[2]) > 0) Y_gw <- map(gw_replica, \(z) rowMeans(z[[3]][, ]))



        # Ricarico (poco efficientemente) tutte le ws e faccio cose, tra cui:
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
        for (m in 1:min(length(fls_replica_k), data_max)) {
          load(fls_replica_k[m])

          # Salvo le nostre stime
          Lambdas_siw[[m]] <- apply(lambda_out, c(2, 3), mean)
          mat_app <- lambda_out
          num_matrices <- dim(mat_app)[1]
          dim_mat <- dim(mat_app)[2]
          log_array <- array(0, dim = c(dim_mat, dim_mat, num_matrices)) # Initialize log array

          for (iii in 1:num_matrices) {
            log_array[, , iii] <- logm(mat_app[iii, , ]) # Apply matrix logarithm to each SPD matrix
          }
          # Compute the element-wise median across all matrices
          log_median <- apply(log_array, c(1, 2), mean)

          # Compute matrix exponential to map back to SPD space
          Lambdas_siw2[[m]] <- expm(log_median)


          Zetas_siw[[m]] <- apply(zmat_out, c(2, 3), mean)

          Zetas_gw[[m]][upper.tri(Zetas_gw[[m]])] <- 1
          diag(Zetas_gw[[m]]) <- 1

          n_el <- (k^2 - k) / 2
          # RMSE #
          if (as.numeric(args[2]) > 0) {
            Y_siw[m, ] <- colMeans(missing_out)
            ytrue <- sapply(1:nrow(index_miss_mat), \(i) y[index_miss_mat[i, 1], index_miss_mat[i, 2]])
            rmse_siw[m] <- sqrt(mean((ytrue - Y_siw[m, ])^2))
            rmse_gw[m] <- sqrt(mean((ytrue - Y_gw[[m]])^2))



            missing_gw <- t(gw_replica[[m]]$ymiss)

            crps_sw[m, ] <- sapply(1:nrow(gw_replica[[m]]$ymiss[, ]), \(i) compute_crps(true_y[i], missing_out[, i]))
            crps_gw[m, ] <- sapply(1:nrow(gw_replica[[m]]$ymiss[, ]), \(i) compute_crps(true_y[i], missing_gw[, i]))
          }

          # Grafici brutali (Gian sii fiero!)
          # Ho un file per ws con i traceplot di tutti i parametri (nero G wishart, rosso noi)
          pdf(paste0("output/", name_fold, "/Traceplots_", m, ".pdf"))
          par(mfrow = c(1, 2))
          for (s in 1:(nrow(lambda) - 1)) {
            for (h in (s + 1):nrow(lambda)) {
              plot(gw_replica[[m]][[2]][h, s, ],
                type = "l",
                ylim = c(
                  min(c(lambda_out[, h, s], gw_replica[[m]][[2]][h, s, ])),
                  max(c(lambda_out[, h, s], gw_replica[[m]][[2]][h, s, ]))
                ),
                main = paste0("L[", h, ", ", s, "], z = ", z_mat[h, s])
              )
              abline(h = lambda[h, s], lwd = 2, col = "darkgreen", lty = 2)
              plot(lambda_out[, h, s],
                type = "l",
                ylim = c(
                  min(c(lambda_out[, h, s], gw_replica[[m]][[2]][h, s, ])),
                  max(c(lambda_out[, h, s], gw_replica[[m]][[2]][h, s, ]))
                ),
                main = paste0("L[", h, ", ", s, "], z = ", z_mat[h, s]), col = 2
              )
              # lines(lambda_out[,h,s], col = 2)
              abline(h = lambda[h, s], lwd = 2, col = "darkgreen", lty = 2)
            }
          }
          dev.off()





          # KLdiscr_siw[m] <- sum(diag(solve(Lambdas_siw2[[m]]) %*% (lambda))) - k - log(det(solve(Lambdas_siw2[[m]]) %*% (lambda)))
          # KLdiscr_gw[m] <- sum(diag(solve(Lambdas_gw2[[m]]) %*% (lambda))) - k - log(det(solve(Lambdas_gw2[[m]]) %*% (lambda)))

          # KLdiscr_siw_inv[m] <- sum(diag((Lambdas_siw2[[m]]) %*% solve(lambda))) - k - log(det((Lambdas_siw2[[m]]) %*% solve(lambda)))
          # KLdiscr_gw_inv[m] <- sum(diag((Lambdas_gw2[[m]]) %*% solve(lambda))) - k - log(det((Lambdas_gw2[[m]]) %*% solve(lambda)))


          # Metrica suggerita qui https://www.tandfonline.com/doi/epdf/10.1080/01621459.2023.2267777?needAccess=true (pag. 9, colonna I)

          KLdiscr_siw[m] <- sum(diag(solve(Lambdas_siw[[m]]) %*% (lambda))) - k - log(det(solve(Lambdas_siw[[m]]) %*% (lambda)))
          KLdiscr_gw[m] <- sum(diag(solve(Lambdas_gw[[m]]) %*% (lambda))) - k - log(det(solve(Lambdas_gw[[m]]) %*% (lambda)))

          KLdiscr_siw_inv[m] <- sum(diag((Lambdas_siw[[m]]) %*% solve(lambda))) - k - log(det((Lambdas_siw[[m]]) %*% solve(lambda)))
          KLdiscr_gw_inv[m] <- sum(diag((Lambdas_gw[[m]]) %*% solve(lambda))) - k - log(det((Lambdas_gw[[m]]) %*% solve(lambda)))

          # kl_app <- rep(NA, dim(gw_replica[[m]]$lambda)[3])
          # for (isim in 1:dim(gw_replica[[m]]$lambda)[3])
          # {
          #  kl_app[isim] <- kl_d(lambda, lambda_out[isim, , ])
          # }
          # KLdiscr_siw[m] <- mean(kl_app)

          # for (isim in 1:dim(gw_replica[[m]]$lambda)[3])
          # {
          #  kl_app[isim] <- kl_d(lambda, gw_replica[[m]]$lambda[, , isim])
          # }
          # KLdiscr_gw[m] <- mean(kl_app)

          # for (isim in 1:dim(gw_replica[[m]]$lambda)[3])
          # {
          #  kl_app[isim] <- kl_inv_d(lambda, lambda_out[isim, , ])
          # }
          # KLdiscr_siw_inv[m] <- mean(kl_app)

          # for (isim in 1:dim(gw_replica[[m]]$lambda)[3])
          # {
          #  kl_app[isim] <- kl_inv_d(lambda, gw_replica[[m]]$lambda[, , isim])
          # }
          # KLdiscr_gw_inv[m] <- mean(kl_app)


          ## https://www.sciencedirect.com/science/article/pii/S0723086918300021
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


          # Compare Z
          # Errore su Z: siamo bravi a stimare la sparsità?
          zerr_siw[m] <- sum(abs(z_mat[lower.tri(z_mat)] - Zetas_siw[[m]][lower.tri(Zetas_siw[[m]])]))
          zerr_gw[m] <- sum(abs(z_mat[lower.tri(z_mat)] - Zetas_gw[[m]][lower.tri(Zetas_gw[[m]])]))

          # Dovrebbe essere giusto, ma da ricontrollare e vedere se si applica a tutti gli scenari
          # percentuale di 0 presi
          sparse_perc_siw[m] <- sum(round(Zetas_siw[[m]]) == 0) / sum(z_mat == 0)
          sparse_perc_gw[m] <- sum(round(Zetas_gw[[m]]) == 0) / sum(z_mat == 0)
          # percentuale di 0 esatti presi
          qsparse_est_siw[m] <- length(intersect(which(z_mat == 0), which(round(Zetas_siw[[m]]) == 0))) / sum(z_mat == 0)
          qsparse_est_gw[m] <- length(intersect(which(z_mat == 0), which(round(Zetas_gw[[m]]) == 0))) / sum(z_mat == 0)


          n_el <- (k^2 - k) / 2

          for (i in 1:dim(zmat_out)[1])
          {
            tt <- table(z_mat[lower.tri(z_mat[, ])], zmat_out[i, , ][lower.tri(zmat_out[i, , ])])
            if (perczero_sim == "0.0") {
              tt <- rbind(c(0, 0), tt)
            }
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



          # Secondo round di grafici brutali
          # top left: lambda vera
          # top center: lambda nostra stimata
          # top right: lambda stimata da G wishart
          # bottom left: differenza tra lambda stimata nostra e vera
          # bottom center: differenza tra lambda stimata con G wishart e vera
          # bottom right: differenza tra lambda nostra a con G wishart
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
          pdf(paste0("output/", name_fold, "/Lambdas_", m, ".pdf"), width = 14, height = 8)
          gridExtra::grid.arrange(grobs = plist, nrow = 2)
          par(mfrow = c(1, 1))
          if (args[2] != "0.0") {
            plot(rmse_siw, rmse_gw, pch = 20)
            abline(a = 0, b = 1, col = 2, lwd = 2)
          }


          dev.off()
        }

        pdf(paste0("Confronti/", name_fold, ".pdf"), width = 14, height = 8)
        par(mfrow = c(2, 2))
        plot(KLdiscr_siw, KLdiscr_gw, pch = 20, xlim = range(c(KLdiscr_siw, KLdiscr_gw)), ylim = range(c(KLdiscr_siw, KLdiscr_gw)), main = paste("SIW ", round(mean(KLdiscr_siw), 3), " GIW ", round(mean(KLdiscr_gw), 3)))
        abline(a = 0, b = 1, col = 2, lwd = 2)

        plot(KLdiscr_siw_inv, KLdiscr_gw_inv, pch = 20, xlim = range(c(KLdiscr_siw_inv, KLdiscr_gw_inv)), ylim = range(c(KLdiscr_siw_inv, KLdiscr_gw_inv)), main = paste("SIW ", round(mean(KLdiscr_siw_inv), 3), " GIW ", round(mean(KLdiscr_gw_inv), 3)))
        abline(a = 0, b = 1, col = 2, lwd = 2)

        plot(BWdiscr_siw, BWdiscr_gw, pch = 20, xlim = range(c(BWdiscr_siw, BWdiscr_gw)), ylim = range(c(BWdiscr_siw, BWdiscr_gw)), main = paste("SIW ", round(mean(BWdiscr_siw), 3), " GIW ", round(mean(BWdiscr_gw), 3)))
        abline(a = 0, b = 1, col = 2, lwd = 2)

        par(mfrow = c(1, 1))
        if (args[2] != "0.0") {
          plot(rmse_siw, rmse_gw, pch = 20, xlim = range(c(rmse_siw, rmse_gw)), ylim = range(c(rmse_siw, rmse_gw)), main = paste("SIW ", round(mean(rmse_siw), 3), " GIW ", round(mean(rmse_gw), 3)))
          abline(a = 0, b = 1, col = 2, lwd = 2)

          plot(c(crps_sw), c(crps_gw), pch = 20, xlim = range(c(rmse_siw, rmse_gw)), ylim = range(c(rmse_siw, rmse_gw)), main = paste("SIW ", round(mean(c(crps_sw)), 3), " GIW ", round(mean(c(crps_gw)), 3)))
          abline(a = 0, b = 1, col = 2, lwd = 2)
        }

        class_index_sw <- tibble(
          acc = rowMeans(acc_sw),
          sens = rowMeans(sens_sw),
          spec = rowMeans(spec_sw), model = "sw"
        )
        class_index_gw <- tibble(
          acc = rowMeans(acc_gw),
          sens = rowMeans(sens_gw),
          spec = rowMeans(spec_gw), model = "gw"
        )

        p <- bind_rows(class_index_sw, class_index_gw) %>%
          pivot_longer(-model, names_to = "metric", values_to = "value") %>%
          ggplot(aes(x = metric, y = value, fill = model)) +
          geom_boxplot(position = "dodge") +
          stat_summary(fun.y = mean, geom = "point", shape = 20, size = 14, position = position_dodge(width = 0.75))

        print(p)

        dev.off()


        save.image(file = paste0("Confronti//CompareGwishart_LikeNormal_", paste0(args, collapse = "_"), ".RData"))
      }
    }
  }
}
