#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Carico i pacchetti
require(tidyverse)
require(magrittr)

# Draw lambda
draw_lambda <- function(mat, is.corr = F){
  
  k <- nrow(mat)
  if(is.corr) mat <- cov2cor(solve(mat))
  
  
  mat <- mat %>% as.data.frame() %>% 
    rownames_to_column(var = "Rows") %>% 
    pivot_longer(-Rows, names_to = "Cols", values_to = "val") %>% 
    mutate(Rows = as.numeric(Rows),  Cols = as.numeric(parse_number(Cols)))
  
  if(!is.corr){
    mat <- mat %>% filter(Rows!=Cols)
    limits <- c(min(mat$val), max(mat$val))
  }else{
    limits <- c(-1, 1)
  }
  
  
  p <- mat %>%
    ggplot() +
    geom_raster(aes(Rows, -Cols, fill = val)) +
    scale_fill_distiller(palette = "YlOrRd", direction = 1, limits = limits) +
    scale_x_continuous(expand = c(0,0), breaks = 1:k, labels = 1:k) + 
    scale_y_continuous(expand = c(0,0), breaks = -(1:k), labels = 1:k) +
    theme_bw() +
    theme(text = element_text(size = 12), legend.position = "top", panel.grid = element_blank())
  
  return(p)
  
}

# Draw zeta
draw_zeta <- function(mat){
  
  k <- nrow(mat)
  
  matout <- mat %>% as.data.frame() %>% 
    rownames_to_column(var = "Rows") %>% 
    pivot_longer(-Rows, names_to = "Cols", values_to = "val") %>% 
    mutate(Zetas = c(round(mat)), Rows = as.numeric(Rows),  Cols = as.numeric(parse_number(Cols)))
  
  p <- matout %>%
    ggplot() +
    geom_raster(aes(Rows, -Cols, fill = factor(Zetas))) +
    scale_fill_manual(values = rev(c("black", "grey80"))) +
    scale_x_continuous(expand = c(0,0), breaks = 1:k, labels = 1:k) + 
    scale_y_continuous(expand = c(0,0), breaks = -(1:k), labels = 1:k) +
    labs(fill= expression(paste(hat(Z))), x = "", y = "") +
    theme_bw() +
    theme(text = element_text(size = 12), legend.position = "top", panel.grid = element_blank())
  
  return(p)
  
}

modtype <- args[1]

if(modtype == "Gene"){
  
  # Normal: gene expression -------------------------------------------------
  k <- args[2]
  load(paste0("out/GeneExp", k, "_Confronti_GW.Rdata")) # cambiare l'RData
  
  # meanzeta for gw and sb
  # sb
  meanzeta <- apply(zmat_out, c(2, 3), mean)
  meanzeta[upper.tri(meanzeta)] <- t(meanzeta)[upper.tri(meanzeta)]
  
  #gw
  meanzeta_gw <- apply(mcmc_gwish$zmat[,,2001:4000], c(1, 2), mean)
  diag(meanzeta_gw) <- 1
  
  # sparsity
  1-sum(round(meanzeta)[lower.tri(meanzeta)]==1)/(k*(k-1)/2) # sb
  1-sum(round(meanzeta_gw)[lower.tri(meanzeta_gw)]==1)/(k*(k-1)/2) # gw
  
  # meanlambda for gw and sb
  # sb
  meanlambda <- apply(lambda_out, c(2, 3), mean)
  
  #gw
  meanlambda_gw <- apply(mcmc_gwish$lambda[,,2001:4000], c(1, 2), mean)
  
  # pdf(pate0("images/EstimatedLambda_Gene",k,"SB.pdf"), width = 8, height = 6)
  # draw_lambda(meanlambda, is.corr = F) + labs(fill= expression(paste(hat(Lambda))), x = "", y = "")
  # dev.off()
  
  # pdf(paste0("images/EstimatedLambda_Gene",k,"GW.pdf"), width = 8, height = 6)
  # draw_lambda(meanlambda_gw, is.corr = F) + labs(fill= expression(paste(hat(Lambda))), x = "", y = "")
  # dev.off()
  
  # pdf(paste0("images/EstimatedZ_Gene", k, "SB.pdf"), width = 8, height = 6)
  # draw_zeta(meanzeta)
  # dev.off()
  
  # pdf(paste0("images/EstimatedZ_Gene", k, "GW.pdf"), width = 8, height = 6)
  # draw_zeta(round(meanzeta_gw))
  # dev.off()
  
  
  ## Rolling ARI and average
  zchain_l <- lapply(1:2000, \(i) c(zmat_out[i,,][lower.tri(meanzeta)]))
  zchain_lgw <- lapply(2001:4000, \(i) c(mcmc_gwish$zmat[,,i][lower.tri(meanzeta)]))
  # aris <- map2_dbl(zchain_l, zchain_lgw, mclust::adjustedRandIndex)
  # mclust::adjustedRandIndex(c(round(meanzeta[lower.tri(meanzeta)])), c(round(meanzeta_gw[lower.tri(meanzeta)])))
  
  
  pdf(paste0("images/CompareZ_Gene", k, "GW.pdf"), width = 8, height = 6)
  round(meanzeta_gw) %>%
    as.data.frame() %>%
    rownames_to_column(var = "Rows") %>%
    pivot_longer(-Rows, names_to = "Cols", values_to = "val") %>%
    mutate(Zetas = c(round(meanzeta)), Rows = as.numeric(Rows),  Cols = as.numeric(parse_number(Cols)),
           FacGW = ifelse(val == 0 & Zetas == 1, "SB = 1 V GW = 0", ifelse(val == 1  & Zetas == 0, "SB = 0 V GW = 1",
                                                                           ifelse(val == 1 & Zetas == 1, "SB = 1 V GW = 1", "SB = 0 V GW = 0")))) %>%
    ggplot() +
    geom_raster(aes(Rows, -Cols, fill = factor(FacGW))) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = c( "grey90", "red3", "green3",  "black")) +
    labs(fill= "", x = "", y = "") +
    theme_bw() +
    theme(text = element_text(size = 18), legend.position = "top", panel.grid = element_blank(), axis.ticks = element_blank(),
          axis.text = element_blank(), legend.key.width = unit(.5,"cm"))
  dev.off()
  
  round(meanzeta_gw) %>%
    as.data.frame() %>%
    rownames_to_column(var = "Rows") %>%
    pivot_longer(-Rows, names_to = "Cols", values_to = "val") %>%
    mutate(Zetas = c(round(meanzeta)), Rows = as.numeric(Rows),  Cols = as.numeric(parse_number(Cols)),
           FacGW = ifelse(val == 0 & Zetas == 1, "SB = 0 V GW = 1", ifelse(val == 1  & Zetas == 0, "SB = 1 V GW = 0",
                                                                           ifelse(val == 1 & Zetas == 1, "SB = 1 V GW = 1", "SB = 0 V GW = 0")))) %>% 
    filter(Rows<Cols) %>% count(FacGW)
  
  pdf(paste0("images/CompareZdistr_Gene", k, "GW.pdf"), width = 8, height = 6)
  tibble(SB = meanzeta[lower.tri(meanzeta)], GW = meanzeta_gw[lower.tri(meanzeta)]) %>% 
    pivot_longer(everything(), names_to = "type", values_to = "val") %>% 
    ggplot() +
    geom_density(aes(x = val, color = type, group = type), linewidth = 1) +
    scale_x_continuous(expand = c(0,0), limits = c(0, 1)) +
    scale_color_manual(values = c("tomato", "skyblue2")) +
    labs(x = expression(paste(hat(z)[jk])), color = "") +
    theme_bw() +
    theme(text = element_text(size = 18), panel.grid = element_blank(), legend.position = "top")
  dev.off()
  
  # Errore di stima sui missing 
  yest <- colMeans(missing_out)
  sqrt(mean((true_y-yest)^2))
  # plot(true_y, yest, pch = 20)
  # abline(0, 1, lwd = 1, col = "red2")
  
  # Se Confronti
  summary(crps_gw)
  summary(crps_sw)
  
  # Traceplots
  pdf(paste0("images/Trace_Gene",k,"_confronto.pdf"), width = 8, height = 6)
  tibble(GW = n_edge_gw, SB = n_edge_sw, iter = 1:nsim) %>% 
    pivot_longer(-iter, names_to = "model", values_to = "val") %>% 
    ggplot(aes(x = iter)) +
    geom_line(aes(y = val, color = model)) +
    scale_color_manual(values = c("tomato", "skyblue2")) +
    labs(x = "iteration", y = "n. of nonzero edges", color = "") +
    theme_bw() +
    theme(legend.position = "top", text = element_text(size = 16))
  dev.off()
  tibble(GW = n_edge_gw, SW = n_edge_sw, iter = 1:nsim) %>% apply(., 2, quantile, probs = c(0.025, 0.975))
  
  # Se Confronti
  summary(n_edge_gw/2)
  summary(n_edge_sw/2)          
  
}else if(modtype == "Doubs"){
  # Poisson: doubs ----------------------------------------------------------
  load("out/doubs_indipendent_speciesmcmc_out_lambda.Rdata") # cambiare l'RData
  
  # doubs <- readxl::read_excel("data/doubs_data.xlsx") %>% filter(Site != 8)
  # plot(doubs$x, doubs$y, pch = 20, type = "n")
  # text(doubs$x, doubs$y, labels = doubs$Site)
  
  # Lambda est
  pdf("images/EstimatedLambda_doubs.pdf", width = 8, height = 6)
  draw_lambda(meanlambda) + labs(fill= expression(paste(hat(Lambda))), x = "", y = "") + theme(legend.key.width = unit(1,"cm"))
  dev.off()
  
  # Zeta est
  pdf("images/EstimatedZeta_doubs.pdf", width = 8, height = 6)
  draw_zeta(round(meanzeta)) 
  dev.off()
  
  # sparsity
  1-sum(round(meanzeta)[lower.tri(meanzeta)]==1)/(k*(k-1)/2) # sb
  
  sapply(1:2000, \(i) sum(zmat_out[i, , ][lower.tri(zmat_out[i, , ])])) %>% summary
  sapply(1:2000, \(i) sum(zmat_out[i, , ][lower.tri(zmat_out[i, , ])])) %>% quantile(., probs = c(0.025, 0.975))
}else{
  stop("the argument must be either Gene or Doubs")
}

