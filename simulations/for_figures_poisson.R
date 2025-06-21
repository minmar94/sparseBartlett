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

for (patt in c("LikePoisson"))
{
  for (k_sim in c(10, 25)[c(2, 1)]) {
    par2 <- 1

    args[1] <- k_sim
    # for(percmiss_sim in c("0.0","0.1")){
    for (percmiss_sim in c("0.1", "0.0")[c(2, 1)])
    {
      args[2] <- percmiss_sim

      for (pattern_sim in c("FullRandom", "Banded1", "Banded2", "Banded3")[c(1, 2, 3, 4)])
      {
        args[3] <- pattern_sim
        if (pattern_sim == "FullRandom") {
          w_miss <- c("0.0", "0.5", "0.25", "0.75")
        } else {
          w_miss <- c("NA")
        }
        for (perczero_sim in w_miss)
        {
          args[4] <- perczero_sim


          fls_replica <- list.files(path = "out_poisson", pattern = patt, full.names = T)

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
          for (m in 1:min(length(fls_replica_k), data_max)) {
            load(fls_replica_k[m])

            # Salvo le nostre stime
            Lambdas_siw[[m]] <- apply(lambda_out, c(2, 3), mean)
            Zetas_siw[[m]] <- apply(zmat_out, c(2, 3), mean)

            # Zetas_gw[[m]][upper.tri(Zetas_gw[[m]])] <- 1
            # diag(Zetas_gw[[m]]) <- 1

            n_el <- (k^2 - k) / 2
            # RMSE #
            if (as.numeric(args[2]) > 0) {
              Y_siw[m, ] <- colMeans(missing_out)
              ytrue <- sapply(1:nrow(index_miss_mat), \(i) y[index_miss_mat[i, 1], index_miss_mat[i, 2]])
              rmse_siw[m] <- sqrt(mean((ytrue - Y_siw[m, ])^2))
              # rmse_gw[m] <- sqrt(mean((ytrue - Y_gw[[m]])^2))



              # missing_gw <- t(gw_replica[[m]]$ymiss)

              crps_sw[m, ] <- sapply(1:ncol((missing_out)), \(i) compute_crps(true_y[i], missing_out[, i]))
              # crps_gw[m, ] <- sapply(1:nrow(gw_replica[[m]]$ymiss[, ]), \(i) compute_crps(true_y[i], missing_gw[, i]))
            }

            # Grafici brutali (Gian sii fiero!)
            # Ho un file per ws con i traceplot di tutti i parametri (nero G wishart, rosso noi)
            # pdf(paste0("output/", name_fold, "/Traceplots_", m, ".pdf"))
            # par(mfrow = c(1, 2))
            # for (s in 1:(nrow(lambda) - 1)) {
            #  for (h in (s + 1):nrow(lambda)) {
            #    plot(gw_replica[[m]][[2]][h, s, ],
            #      type = "l",
            #      ylim = c(
            #        min(c(lambda_out[, h, s], gw_replica[[m]][[2]][h, s, ])),
            #        max(c(lambda_out[, h, s], gw_replica[[m]][[2]][h, s, ]))
            #      ),
            #      main = paste0("L[", h, ", ", s, "], z = ", z_mat[h, s])
            #    )
            #    abline(h = lambda[h, s], lwd = 2, col = "darkgreen", lty = 2)
            #    plot(lambda_out[, h, s],
            #      type = "l",
            #      ylim = c(
            #        min(c(lambda_out[, h, s], gw_replica[[m]][[2]][h, s, ])),
            #        max(c(lambda_out[, h, s], gw_replica[[m]][[2]][h, s, ]))
            #      ),
            #      main = paste0("L[", h, ", ", s, "], z = ", z_mat[h, s]), col = 2
            #    )
            #    # lines(lambda_out[,h,s], col = 2)
            #    abline(h = lambda[h, s], lwd = 2, col = "darkgreen", lty = 2)
            #  }
            # }
            # dev.off()


            # KL Diverìgence
            # Metrica suggerita qui https://www.tandfonline.com/doi/epdf/10.1080/01621459.2023.2267777?needAccess=true (pag. 9, colonna I)
            KLdiscr_siw[m] <- sum(diag(solve(Lambdas_siw[[m]]) %*% (lambda))) - k - log(det(solve(Lambdas_siw[[m]]) %*% (lambda)))
            # KLdiscr_gw[m] <- sum(diag(solve(Lambdas_gw[[m]]) %*% (lambda))) - k - log(det(solve(Lambdas_gw[[m]]) %*% (lambda)))

            KLdiscr_siw_inv[m] <- sum(diag((Lambdas_siw[[m]]) %*% solve(lambda))) - k - log(det((Lambdas_siw[[m]]) %*% solve(lambda)))
            #  KLdiscr_gw_inv[m] <- sum(diag((Lambdas_gw[[m]]) %*% solve(lambda))) - k - log(det((Lambdas_gw[[m]]) %*% solve(lambda)))


            ## https://www.sciencedirect.com/science/article/pii/S0723086918300021
            app_eig <- eigen(Lambdas_siw[[m]])
            lambda_eig <- app_eig$vectors %*% diag(app_eig$values^0.5) %*% t(app_eig$vectors)
            app_molt_eig <- eigen(lambda_eig %*% lambda %*% lambda_eig)
            molt_lambda_eig <- app_molt_eig$vectors %*% diag(app_molt_eig$values^0.5) %*% t(app_molt_eig$vectors)
            BWdiscr_siw[m] <- sum(diag(Lambdas_siw[[m]] + lambda - 2 * molt_lambda_eig))


            # app_eig <- eigen(Lambdas_gw[[m]])
            # lambda_eig <- app_eig$vectors %*% diag(app_eig$values^0.5) %*% t(app_eig$vectors)
            # app_molt_eig <- eigen(lambda_eig %*% lambda %*% lambda_eig)
            # molt_lambda_eig <- app_molt_eig$vectors %*% diag(app_molt_eig$values^0.5) %*% t(app_molt_eig$vectors)
            # BWdiscr_gw[m] <- sum(diag(Lambdas_gw[[m]] + lambda - 2 * molt_lambda_eig))


            # Compare Z
            # Errore su Z: siamo bravi a stimare la sparsità?
            zerr_siw[m] <- sum(abs(z_mat[lower.tri(z_mat)] - Zetas_siw[[m]][lower.tri(Zetas_siw[[m]])]))
            #  zerr_gw[m] <- sum(abs(z_mat[lower.tri(z_mat)] - Zetas_gw[[m]][lower.tri(Zetas_gw[[m]])]))

            # Dovrebbe essere giusto, ma da ricontrollare e vedere se si applica a tutti gli scenari
            # percentuale di 0 presi
            sparse_perc_siw[m] <- sum(round(Zetas_siw[[m]]) == 0) / sum(z_mat == 0)
            # sparse_perc_gw[m] <- sum(round(Zetas_gw[[m]]) == 0) / sum(z_mat == 0)
            # percentuale di 0 esatti presi
            qsparse_est_siw[m] <- length(intersect(which(z_mat == 0), which(round(Zetas_siw[[m]]) == 0))) / sum(z_mat == 0)
            # qsparse_est_gw[m] <- length(intersect(which(z_mat == 0), which(round(Zetas_gw[[m]]) == 0))) / sum(z_mat == 0)


            n_el <- (k^2 - k) / 2

            for (i in 1:dim(zmat_out)[1])
            {
              tt <- table(factor(z_mat[lower.tri(z_mat[, ])], levels = c(0, 1)), factor(zmat_out[i, , ][lower.tri(zmat_out[i, , ])], levels = c(0, 1)))


              # if (perczero_sim == "0.0") {
              #  tt <- rbind(c(0, 0), tt)
              # }
              # if(dim(tt)[1]==1)
              # {

              # }



              acc_sw[m, i] <- (tt[1, 1] + tt[2, 2]) / sum(tt)
              sens_sw[m, i] <- tt[1, 1] / (tt[1, 1] + tt[1, 2])
              spec_sw[m, i] <- (tt[2, 2]) / (tt[2, 2] + tt[2, 1])

              # tt <- table(z_mat[lower.tri(z_mat[, ])], gw_replica[[m]]$zmat[, , -c(1)][, , i][lower.tri(zmat_out[i, , ])])
              # if (perczero_sim == "0.0") {
              #  tt <- rbind(c(0, 0), tt)
              # }
              # acc_gw[m, i] <- (tt[1, 1] + tt[2, 2]) / sum(tt)
              # sens_gw[m, i] <- tt[1, 1] / (tt[1, 1] + tt[1, 2])
              # spec_gw[m, i] <- (tt[2, 2]) / (tt[2, 2] + tt[2, 1])
            }



            # Secondo round di grafici brutali
            # top left: lambda vera
            # top center: lambda nostra stimata
            # top right: lambda stimata da G wishart
            # bottom left: differenza tra lambda stimata nostra e vera
            # bottom center: differenza tra lambda stimata con G wishart e vera
            # bottom right: differenza tra lambda nostra a con G wishart
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
            # pdf(paste0("output/", name_fold, "/Lambdas_", m, ".pdf"), width = 14, height = 8)
            # gridExtra::grid.arrange(grobs = plist, nrow = 2)
            # par(mfrow = c(1, 1))
            # if (args[2] != "0.0") {
            #  plot(rmse_siw, rmse_gw, pch = 20)
            #  abline(a = 0, b = 1, col = 2, lwd = 2)
            # }


            # dev.off()
          }

          # pdf(paste0("glm/", patt, "_", name_fold, ".pdf"), width = 14, height = 8)


          # class_index_sw <- tibble(
          #  acc = rowMeans(acc_sw),
          #  sens = rowMeans(sens_sw),
          #  spec = rowMeans(spec_sw), model = "sw"
          # )


          # p <- class_index_sw %>%
          #  pivot_longer(-model, names_to = "metric", values_to = "value") %>%
          #  ggplot(aes(x = metric, y = value, fill = model)) +
          #  geom_boxplot(position = "dodge") +
          #  stat_summary(fun.y = mean, geom = "point", shape = 20, size = 14, position = position_dodge(width = 0.75))

          # print(p)

          # dev.off()


          save.image(file = paste0(getwd(), "/for_figures_poisson/", patt, "_", name_fold, ".Rdata"))
        }
      }
    }
  }
}












# install.packages("image.darknet", repos = "https://bnosac.github.io/drat")
