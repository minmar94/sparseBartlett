#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Carico i pacchetti
require(tidyverse)
require(magrittr)

setwd("./simulations/")
fls_list <- list.files(path = "for_figures_normal/", pattern = "Normal", full.names = T)
modtype <- args[1]
fls_type <- fls_list[grepl(modtype, fls_list)]

class_index_all <- KL_all <- crps_all <- Nedges <- list()
if(modtype == "Banded"){
  for(ks in c(10, 25)){
    for(perm in c("0.0", "0.1")){
      fls_sub <- fls_type[grepl(paste0(ks, "_", perm), fls_type)]
      for(jj in 1:length(fls_sub)){
        load(fls_sub[jj])
        
        # Recovery of the sparsity
        class_index_all <- bind_rows(class_index_all, 
                                     bind_rows(
                                       class_index_gw, class_index_sw
                                     ) %>% mutate(k = ks, missing = perm, modtype = paste0(modtype,jj)))
        # KL-divergence
        KL_all <- bind_rows(KL_all, 
                            bind_rows(
                              tibble(GW = KLdiscr_gw, SB = KLdiscr_siw, k = ks,
                                     missing = perm, modtype = paste0(modtype,jj))
                            ))
        # Missing data
        crps_all <- bind_rows(crps_all, bind_rows(
          tibble(GW = rowMeans(crps_gw), SB = rowMeans(crps_sw), 
                 k = ks, missing = perm, modtype = paste0(modtype,jj))
        ))
        # N. of edges
        n_edges_GW <- sapply(1:2000, \(zz) sum(gw_replica[[20]]$zmat[,,zz][lower.tri(gw_replica[[20]]$zmat[,,zz])]))
        n_edges_SB <- sapply(1:2000, \(zz) sum(zmat_out[zz,,][lower.tri(zmat_out[zz,,])]))
        Nedges <- bind_rows(Nedges, tibble(
          SB = n_edges_SB
          ,
          GW = n_edges_GW) %>%
            mutate(k = ks, missing = perm, modtype = paste0(modtype,jj), iter = 1:2000,
                   Ntrue = sum(z_mat[lower.tri(z_mat)])))
      }
    }
  }
}else{
  for(ks in c(10, 25)){
    for(perm in c("0.0", "0.1")){
      fls_sub <- fls_type[grepl(paste0(ks, "_", perm), fls_type)]
      for(jj in 1:length(fls_sub)){
        load(fls_sub[jj])
        
        # Recovery of the sparsity
        class_index_all <- bind_rows(class_index_all, 
                                     bind_rows(
                                       class_index_gw, class_index_sw
                                     ) %>% mutate(k = ks, missing = perm, percsparse = c(0, 0.25, 0.5, 0.75)[jj]))
        # KL-divergence
        KL_all <- bind_rows(KL_all, 
                            bind_rows(
                              tibble(GW = KLdiscr_gw, SB = KLdiscr_siw, k = ks,
                                     missing = perm, percsparse = c(0, 0.25, 0.5, 0.75)[jj])
                            ))
        # Missing data
        crps_all <- bind_rows(crps_all, bind_rows(
          tibble(GW = rowMeans(crps_gw), SB = rowMeans(crps_sw), 
                 k = ks, missing = perm, percsparse = c(0, 0.25, 0.5, 0.75)[jj])
        ))
        # N. of edges
        n_edges_GW <- sapply(1:2000, \(zz) sum(gw_replica[[20]]$zmat[,,zz][lower.tri(gw_replica[[20]]$zmat[,,zz])]))
        n_edges_SB <- sapply(1:2000, \(zz) sum(zmat_out[zz,,][lower.tri(zmat_out[zz,,])]))
        Nedges <- bind_rows(Nedges, tibble(
          SB = n_edges_SB,
          GW = n_edges_GW) %>%
          mutate(k = ks, missing = perm, percsparse = c(0, 0.25, 0.5, 0.75)[jj], 
                 modtype = paste0(modtype,jj), iter = 1:2000,
                 Ntrue = sum(z_mat[lower.tri(z_mat)])))
      }
    }
  }
}


if(modtype == "Banded"){
  # Banded ------------------------------------------------------------------
  #### Recover the sparsity ####
  # PLOT
  pdf("images/Sim_normal_banded_k_miss_sparse_comparesensitivity.pdf", width = 12, height = 7)
  class_index_all %>% 
    pivot_longer(acc:spec, names_to = "classperf", values_to = "val") %>% 
    filter(classperf == "sens") %>% 
    mutate(model = ifelse(model == "gw", "GW", "SB"),
           classperf = ifelse(classperf == "acc", "accuracy", 
                              ifelse(classperf == "sens", "sensitivity", "specificity")),
           modtype = gsub("Banded", "band_", modtype)) %>% 
    ggplot(aes(x = factor(k), y = val, fill = model)) +
    geom_boxplot(position = "dodge") +
    facet_grid(rows = vars(missing), cols = vars(modtype), switch = "y") +
    scale_fill_manual(values = c("tomato", "skyblue2")) +
    labs(x = "p", fill = "", y = "sensitivity") +
    theme_bw() +
    theme(legend.position = "top", text = element_text(size = 16))
  dev.off()
  
  # TABLE
  # class_index_all %>% 
  #   filter(modtype == "Banded1") %>% 
  #   pivot_longer(acc:spec, names_to = "classperf", values_to = "val") %>% 
  #   group_by(k, missing, modtype, model, classperf) %>% 
  #   skimr::skim(val)
  
  #### KL-discrepancy per lambda ####
  # PLOT
  pdf("images/Sim_normal_banded_k_miss_sparse_compareKL.pdf", width = 12, height = 7)
  KL_all %>% 
    pivot_longer(GW:SB, names_to = "model", values_to = "val") %>% 
    mutate(modtype = gsub("Banded", "band_", modtype)) %>% 
    ggplot(aes(x = factor(k), y = log(val), fill = model)) +
    geom_boxplot(position = "dodge") +
    facet_grid(rows = vars(missing), cols = vars(modtype), switch = "y") +
    scale_fill_manual(values = c("tomato", "skyblue2")) +
    labs(x = "p", fill = "", y = "KL-discrepancy (log)") +
    theme_bw() +
    theme(legend.position = "top", text = element_text(size = 16))
  dev.off()
  # TABLE
  # KL_all %>% 
  #   pivot_longer(GW:SB, names_to = "model", values_to = "val") %>% 
  #   group_by(k, missing, modtype, model) %>% 
  #   skimr::skim(val)
  
  #### CRPS ####
  # PLOT
  # crps_all %>% mutate(modtype = gsub("Banded", "band_", modtype)) %>% 
  #   filter(missing == "0.1") %>% 
  #   pivot_longer(GW:SB, names_to = "model", values_to = "val") %>% 
  #   ggplot(aes(x = factor(k), y = val, fill = model)) +
  #   geom_boxplot(position = "dodge") +
  #   facet_grid(cols = vars(modtype), switch = "y") +
  #   scale_fill_manual(values = c("tomato", "skyblue2")) +
  #   labs(x = "k", fill = "", y = "CRPS") +
  #   theme_bw() +
  #   theme(legend.position = "top", text = element_text(size = 16))
  # TABLE
  crps_all %>% 
    filter(missing == "0.1") %>% 
    pivot_longer(GW:SB, names_to = "model", values_to = "val") %>% 
    group_by(k, missing, modtype, model) %>% 
    skimr::skim(val)
  
  
  ### Traceplots ###
  pdf("images/Sim_normal_banded_k_miss_sparse_compareNZero.pdf", width = 12, height = 7)
  Nedges %>% 
    # mutate(GW = (k^2-k)/2-GW, SB = (k^2-k)/2-SB) %>% 
    mutate(modtype = gsub("Banded", "band_", modtype)) %>% 
    filter(missing == "0.0") %>% 
    pivot_longer(SB:GW, names_to = "model", values_to = "val") %>% 
    ggplot(aes(x = iter, y = val, color = model, group = model)) +
    geom_line() +
    scale_color_manual(values = c("tomato", "skyblue2")) +
    labs(x = "iteration", y = "N. of nonzero edges", color = "") +
    facet_grid(rows = vars(k), cols = vars(modtype), switch = "y", scales = "free_y") +
    theme_bw() +
    theme(legend.position = "top", text = element_text(size = 16), axis.text.x = element_text(size = 12))
  dev.off()
  
  Nedges %>% 
    filter(missing == "0.0") %>% 
    pivot_longer(SB:GW, names_to = "model", values_to = "val") %>% 
    group_by(model, modtype, k) %>% 
    summarise(
      True = first(Ntrue), Mean = round(mean(val)), q1 = round(quantile(val, 0.025)), q2 = round(quantile(val, 0.975))
    ) %>% 
    ungroup() %>% 
    mutate(Est = paste0(Mean, " (", q1, ", ", q2, ")"), model = paste(model, modtype, sep="_")) %>% 
    select(model, k, Est) %>% pivot_wider(names_from = model, values_from = Est)
  
}else if(modtype == "Random"){
  # Random ------------------------------------------------------------------
  #### Recover the sparsity ####
  # PLOT
  pdf("images/Sim_normal_random_k_miss_sparse_comparesensitivity.pdf", width = 12, height = 7)
  class_index_all %>% 
    pivot_longer(acc:spec, names_to = "classperf", values_to = "val") %>% 
    filter(classperf == "sens", percsparse>0) %>% 
    mutate(model = ifelse(model == "gw", "GW", "SB"),
           classperf = ifelse(classperf == "acc", "accuracy", 
                              ifelse(classperf == "sens", "sensitivity", "specificity")),
           percsparse = paste0("rand_", percsparse)) %>% 
    ggplot(aes(x = factor(k), y = val, fill = model)) +
    geom_boxplot(position = "dodge") +
    facet_grid(rows = vars(missing), cols = vars(percsparse), switch = "y") +
    scale_fill_manual(values = c("tomato", "skyblue2")) +
    labs(x = "p", fill = "", y = "sensitivity") +
    theme_bw() +
    theme(legend.position = "top", text = element_text(size = 16))
  dev.off()
  # TABLE
  # class_index_all %>% 
  #   pivot_longer(acc:spec, names_to = "classperf", values_to = "val") %>% 
  #   filter(classperf == "acc") %>% 
  #   group_by(k, missing, model, percsparse) %>% 
  #   skimr::skim(val)
  
  #### KL-discrepancy per lambda ####
  # PLOT
  pdf("images/Sim_normal_random_k_miss_sparse_compareKL.pdf", width = 12, height = 7)
  KL_all %>% mutate(percsparse = paste0("rand_", percsparse)) %>% 
    pivot_longer(GW:SB, names_to = "model", values_to = "val") %>% 
    ggplot(aes(x = factor(k), y = log(val), fill = model)) +
    geom_boxplot(position = "dodge") +
    facet_grid(rows = vars(missing), cols = vars(percsparse), switch = "y") +
    scale_fill_manual(values = c("tomato", "skyblue2")) +
    labs(x = "p", fill = "", y = "KL-discrepancy (log)") +
    theme_bw() +
    theme(legend.position = "top", text = element_text(size = 16))
  dev.off()
  # TABLE
  # KL_all %>% 
  #   pivot_longer(GW:SB, names_to = "model", values_to = "val") %>% 
  #   filter(val<10) %>% 
  #   group_by(k, missing, percsparse, model) %>% 
  #   skimr::skim(val)
  
  #### CRPS ####
  # PLOT
  # pdf("images/Sim_normal_random_k_miss_sparse_comparecrps.pdf", width = 12, height = 7)
  # crps_all %>% 
  #   filter(missing == "0.1") %>% 
  #   pivot_longer(GW:SB, names_to = "model", values_to = "val") %>% 
  #   ggplot(aes(x = factor(k), y = val, fill = model)) +
  #   geom_boxplot(position = "dodge") +
  #   facet_grid(cols = vars(percsparse), switch = "y") +
  #   scale_fill_manual(values = c("tomato", "skyblue2")) +
  #   labs(x = "k", fill = "", y = "CRPS") +
  #   theme_bw() +
  #   theme(legend.position = "top", text = element_text(size = 16))
  # dev.off()
  # TABLE
  crps_all %>% 
    filter(missing == "0.1") %>% 
    pivot_longer(GW:SB, names_to = "model", values_to = "val") %>% 
    group_by(k, missing, percsparse, model) %>% 
    skimr::skim(val)
  
  #### TRACEPLOTS ####
  pdf("images/Sim_normal_random_k_miss_sparse_compareNZero.pdf", width = 12, height = 7)
  Nedges %>% 
    filter(missing == "0.0") %>% 
    mutate(percsparse = paste0("rand_", percsparse)) %>% 
    pivot_longer(SB:GW, names_to = "model", values_to = "val") %>% 
    ggplot(aes(x = iter, y = val, color = model, group = model)) +
    geom_line() +
    scale_color_manual(values = c("tomato", "skyblue2")) +
    labs(x = "iteration", y = "N. of nonzero edges", color = "") +
    facet_grid(rows = vars(k), cols = vars(percsparse), switch = "y", scales = "free_y") +
    theme_bw() +
    theme(legend.position = "top", text = element_text(size = 16), axis.text.x = element_text(size = 12))
  dev.off()
  
  Nedges %>% 
    filter(missing == "0.0") %>% 
    pivot_longer(SB:GW, names_to = "model", values_to = "val") %>% 
    group_by(model, percsparse, k) %>% 
    summarise(
      True = first(Ntrue), Mean = round(mean(val)), q1 = round(quantile(val, 0.025)), q2 = round(quantile(val, 0.975))
    ) %>% 
    ungroup() %>% 
    mutate(Est = paste0(Mean, " (", q1, ", ", q2, ")"), model = paste(model, percsparse, sep="_")) %>% 
    select(model, k, Est) %>% pivot_wider(names_from = model, values_from = Est)
  
}else{
  stop("The argument must be either Banded or Random")
}


