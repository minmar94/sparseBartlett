require(tidyverse)
require(magrittr)

knitr::opts_chunk$set(echo = F, message = F, fig.width = 16, warning = F)


#' ## Poisson

#+ fls_pois
fls_pois <- list.files(path = "Simulations", pattern = args[1], full.names = T)

tabout_pois <- map_dfr(1:length(fls_pois), \(l){
  load(fls_pois[l])
  
  fls_pois_split <- str_split(fls_pois[l], "_")[[1]]
  percmiss_oracle <- parse_number(fls_pois_split[grep("percmiss", fls_pois_split)])
  k <- parse_number(fls_pois_split[grep("k=", fls_pois_split)])
  perczero_oracle <- as.numeric(str_extract(fls_pois_split[grep("perczero", fls_pois_split)],
                                            "(?<=perczero=)[0-9.]+(?=decay)"))
  percdecay_oracle <- as.numeric(str_extract(fls_pois_split[grep("decay", fls_pois_split)],
                                             "(?<=decay=)[0-9.]+(?=seed)"))
  # Distanza tra matrici
  lambda_est <- apply(lambda_out, c(2, 3), mean)
  Frob_dist <- SMFilter::FDist2(lambda, lambda_est)
  KLdiscr <- sum(diag(lambda_est%*%solve(lambda))) - k - log(det(lambda_est%*%solve(lambda)))
  # Perc sparsity
  # psparse_true <- mean(!z_mat[lower.tri(z_mat)])
  # zmat_est <- matrix(data_gg$meanzeta, nrow = k, byrow = T)
  zmat_est <- apply(zmat_out, c(2,3), mean)
  zmat_est <- round(zmat_est)
  psparse_est <- mean(!zmat_est)

  # Quali zeri?
  qsparse_true = sum(z_mat == 0)
  qsparse_est <- length(intersect(which(z_mat == 0), which(zmat_est == 0)))
  
  
  # Errore stima missing
  if(percmiss_oracle>0){
    ymiss <- y %>% as.data.frame() %>%
      set_colnames(value = 1:200) %>% 
      rownames_to_column(var = "V1") %>% 
      pivot_longer(-V1, names_to = "V2", values_to = "val") %>% 
      mutate(V1 = as.numeric(V1), V2 = as.numeric(V2)) %>% 
      right_join(as.data.frame(index_miss_mat)) %>% 
      arrange(V2, V1) %>% 
      mutate(missest = colMeans(missing_out))
    
    rmse <- sqrt(mean(ymiss$val - ymiss$missest)^2)
    
  }else{
    rmse <- NA
  }
  
  tibble(
    Model = "Poisson",
    k = k,
    MissingTrue = percmiss_oracle,
    ZeroTrue = perczero_oracle,
    Decay = percdecay_oracle,
    # ZeroSim = psparse_true,
    ZeroEst = psparse_est,
    QualiZeroTrue = qsparse_true,
    QualiZero = qsparse_est,
    FrobLambda = Frob_dist,
    KlD = KLdiscr,
    RmseMiss = rmse,
    Time = elapsed
  )
})


#+ tabout_pois
knitr::kable(tabout_pois, digits = 3, align = "c") %>% 
kableExtra::kable_styling(bootstrap_options = "striped", full_width = F)

 
#' ## Normal

#+ fls_norm
fls_norm <- list.files(path = "Simulations", pattern = "LikeNorm", full.names = T)

tabout_norm <- map_dfr(1:length(fls_norm), \(l){
  load(fls_norm[l])
  
  fls_norma_split <- str_split(fls_norm[l], "_")[[1]]
  percmiss_oracle <- parse_number(fls_norma_split[grep("percmiss", fls_norma_split)])
  k <- parse_number(fls_norma_split[grep("k=", fls_norma_split)])
  perczero_oracle <- as.numeric(str_extract(fls_norma_split[grep("perczero", fls_norma_split)],
                                            "(?<=perczero=)[0-9.]+(?=decay)"))
  percdecay_oracle <- as.numeric(str_extract(fls_norma_split[grep("decay", fls_norma_split)],
                                             "(?<=decay=)[0-9.]+(?=mcmc)"))
  # Distanza tra matrici
  lambda_est <- apply(lambda_out, c(2, 3), mean)
  Frob_dist <- SMFilter::FDist2(lambda, lambda_est)
  KLdiscr <- sum(diag(lambda_est%*%solve(lambda))) - k - log(det(lambda_est%*%solve(lambda)))
  # Perc sparsity
  psparse_true <- mean(!z_mat)
  # psparse_est <- mean(!data_gg$meanzeta)
  # zmat_est <- matrix(data_gg$meanzeta, nrow = k, byrow = T)
  zmat_est <- apply(zmat_out, c(2,3), mean)
  zmat_est <- round(zmat_est)
  psparse_est <- mean(!zmat_est)
  # Quali zeri?
  qsparse_true = sum(z_mat == 0)
  qsparse_est <- length(intersect(which(z_mat == 0), which(zmat_est == 0)))
  
  
  # Errore stima missing
  if(percmiss_oracle>0){
    ymiss <- y %>% as.data.frame() %>%
      set_colnames(value = 1:200) %>% 
      rownames_to_column(var = "V1") %>% 
      pivot_longer(-V1, names_to = "V2", values_to = "val") %>% 
      mutate(V1 = as.numeric(V1), V2 = as.numeric(V2)) %>% 
      right_join(as.data.frame(index_miss_mat)) %>% 
      arrange(V2, V1) %>% 
      mutate(missest = colMeans(missing_out))
    
    rmse <- sqrt(mean(ymiss$val - ymiss$missest)^2)
    
  }else{
    rmse <- NA
  }
  
  tibble(
    Model = "Normal",
    k = k,
    MissingTrue = percmiss_oracle,
    ZeroTrue = perczero_oracle,
    Decay = percdecay_oracle,
    # ZeroSim = psparse_true,
    ZeroEst = psparse_est,
    QualiZeroTrue = qsparse_true,
    QualiZero = qsparse_est,
    FrobLambda = Frob_dist,
    KlD = KLdiscr,
    RmseMiss = rmse,
    Time = elapsed
  )
})


#+ tabout_norm
knitr::kable(tabout_norm, digits = 3, align = "c") %>% 
  kableExtra::kable_styling(bootstrap_options = "striped", full_width = F)


#' ## Bernoulli

#+ fls_ber
fls_ber <- list.files(path = "Simulations", pattern = "LikeBer", full.names = T)

tabout_ber <- map_dfr(1:length(fls_ber), \(l){
  load(fls_ber[l])
  
  fls_ber_split <- str_split(fls_ber[l], "_")[[1]]
  percmiss_oracle <- parse_number(fls_ber_split[grep("percmiss", fls_ber_split)])
  k <- parse_number(fls_ber_split[grep("k=", fls_ber_split)])
  perczero_oracle <- as.numeric(str_extract(fls_ber_split[grep("perczero", fls_ber_split)],
                                            "(?<=perczero=)[0-9.]+(?=decay)"))
  percdecay_oracle <- as.numeric(str_extract(fls_ber_split[grep("decay", fls_ber_split)],
                                             "(?<=decay=)[0-9.]+(?=mcmc)"))
  # Distanza tra matrici
  lambda_est <- apply(lambda_out, c(2, 3), mean)
  Frob_dist <- SMFilter::FDist2(lambda, lambda_est)
  KLdiscr <- sum(diag(lambda_est%*%solve(lambda))) - k - log(det(lambda_est%*%solve(lambda)))
  # Perc sparsity
  psparse_true <- mean(!z_mat)
  # psparse_est <- mean(!data_gg$meanzeta)
  # zmat_est <- matrix(data_gg$meanzeta, nrow = k, byrow = T)
  zmat_est <- apply(zmat_out, c(2,3), mean)
  zmat_est <- round(zmat_est)
  psparse_est <- mean(!zmat_est)
  # Quali zeri?
  qsparse_true = sum(z_mat == 0)
  qsparse_est <- length(intersect(which(z_mat == 0), which(zmat_est == 0)))
  
  
  # Errore stima missing
  if(percmiss_oracle>0){
    ymiss <- y %>% as.data.frame() %>%
      set_colnames(value = 1:200) %>% 
      rownames_to_column(var = "V1") %>% 
      pivot_longer(-V1, names_to = "V2", values_to = "val") %>% 
      mutate(V1 = as.numeric(V1), V2 = as.numeric(V2)) %>% 
      right_join(as.data.frame(index_miss_mat)) %>% 
      arrange(V2, V1) %>% 
      mutate(missest = colMeans(missing_out))
    
    rmse <- sqrt(mean(ymiss$val - ymiss$missest)^2)
    
  }else{
    rmse <- NA
  }
  
  tibble(
    Model = "Bernoulli",
    k = k,
    MissingTrue = percmiss_oracle,
    ZeroTrue = perczero_oracle,
    Decay = percdecay_oracle,
    # ZeroSim = psparse_true,
    ZeroEst = psparse_est,
    QualiZeroTrue = qsparse_true,
    QualiZero = qsparse_est,
    FrobLambda = Frob_dist,
    KlD = KLdiscr,
    RmseMiss = rmse,
    Time = elapsed
  )
})


#+ tabout_ber
knitr::kable(tabout_ber, digits = 3, align = "c") %>% 
  kableExtra::kable_styling(bootstrap_options = "striped", full_width = F)
  
# data_gg # qui ci sono i dati di sintesi
# index_miss_mat # matrice con indici dove ci sono i missing
# lambda # lambda vera
# y # dati usati per la stima
# z_mat # matrice con indici di sparsità vera
# 
# lambda_out # catene di lambda
# zmat_out # catene di z
# missing_out # catene dei missing
# 
# elapsed # tempo per la stima




mu = c(1,2,3)

matrix(mu, ncol=3, nrow=100, byrow=T)