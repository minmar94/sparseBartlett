# Carico i pacchetti
require(tidyverse)
require(magrittr)

liktype <- "Poisson"
fls_list <- list.files(path = "Simulations/ConfrontiBuoni/", pattern = liktype, full.names = T)
modtype <- "Random" 
fls_type <- fls_list[grepl(modtype, fls_list)]

class_index_all <- KL_all <- crps_all <- list()
for(ks in c(10, 25)){
  for(perm in c("0.0", "0.1")){
    fls_sub <- fls_type[grepl(paste0(ks, " percmiss=", perm), fls_type)]
    for(jj in 1:length(fls_sub)){
      load(fls_sub[jj])
      class_index_all <- bind_rows(class_index_all, 
                                   class_index_sw %>% 
                                     mutate(k = ks, missing = perm, 
                                            modtype = paste0(modtype,jj)))
      
      KL_all <- bind_rows(KL_all, tibble(SB = KLdiscr_siw, k = ks,
                                 missing = perm, modtype = paste0(modtype,jj))
      )
      crps_all <- bind_rows(crps_all, 
        tibble(SB = rowMeans(crps_sw), 
               k = ks, missing = perm, modtype = paste0(modtype,jj)))
    }
  }
}




# gli oggetti di dimensione 2000 sono i risultati ad ogni iterazione dell'MCMC
# mentre il numero di repliche per ciascun scenario è pari a 20

# LA VERITA'
# lambda
# index_miss_mat
# y
# z_mat



# Random ------------------------------------------------------------------
#### Recover the sparsity ####
# PLOT
pdf("Images/Sim_poisson_random_k_miss_sparse_classperf.pdf", width = 12, height = 7)
class_index_all %>% 
  pivot_longer(acc:spec, names_to = "classperf", values_to = "val") %>% 
  mutate(classperf = ifelse(classperf == "acc", "accuracy", 
                            ifelse(classperf == "sens", "sensitivity", "specificity")),
         modtype = ifelse(modtype == "Random1", 0, 
                          ifelse(modtype == "Random2", 0.25, 
                                 ifelse(modtype == "Random3", 0.5, 0.75))),
         modtype = paste0("rand_", modtype)) %>% 
  filter(classperf == "sensitivity", modtype !="rand_0") %>% 
  ggplot(aes(x = factor(k), y = val, fill = I("skyblue2"))) +
  geom_boxplot(position = "dodge") +
  facet_grid(rows = vars(missing), cols = vars(modtype), switch = "y") +
  # scale_fill_manual(values = c( "tomato", "skyblue2", "darkgreen")) +
  labs(x = "p", y = "sensitivity", fill = "sparsity (%)") +
  theme_bw() +
  theme(legend.position = "top", text = element_text(size = 16))
dev.off()
# TABLE
class_index_all %>% 
  pivot_longer(acc:spec, names_to = "classperf", values_to = "val") %>% 
  group_by(k, missing, modtype, classperf) %>% 
  skimr::skim(val)

#### KL-discrepancy per lambda ####
# PLOT
pdf("Images/Sim_poisson_random_k_miss_sparse_KL.pdf", width = 12, height = 7)
KL_all %>% 
  mutate(
    modtype = ifelse(modtype == "Random1", 0, 
                     ifelse(modtype == "Random2", 0.25, 
                            ifelse(modtype == "Random3", 0.5, 0.75))),
    modtype = paste0("rand_", modtype)
  ) %>% 
  ggplot(aes(x = factor(k), y = log(SB), fill = I("skyblue2"))) +
  geom_boxplot(position = "dodge") +
  facet_grid(rows = vars(missing), cols = vars(modtype), switch = "y") +
  labs(x = "p", fill = "sparsity (%)", y = "KL-discrepancy (log)") +
  theme_bw() +
  theme(legend.position = "top", text = element_text(size = 16))
dev.off()
# TABLE
KL_all %>% 
  group_by(k, missing, modtype) %>% 
  skimr::skim(SB)

#### CRPS ####
# PLOT
jpeg("Images/Sim_poisson_random_k_miss_sparse_crps.jpeg", width = 800, height = 600, res = 150)
crps_all %>% 
  filter(missing == "0.1") %>% 
  ggplot(aes(x = factor(k), y = SB, fill = modtype)) +
  geom_boxplot(position = "dodge") +
  scale_fill_manual(values = c("tomato", "skyblue2", "darkgreen")) +
  labs(x = "k", fill = "", y = "CRPS") +
  theme_bw() +
  theme(legend.position = "top", text = element_text(size = 16))
dev.off()
# TABLE
crps_all %>% 
  filter(missing == "0.1") %>% 
  group_by(k, missing, modtype) %>% 
  skimr::skim(SB)



# Banded ------------------------------------------------------------------
#### Recover the sparsity ####
# PLOT
pdf("Images/Sim_poisson_banded_k_miss_sparse_classperf.pdf", width = 12, height = 7)
class_index_all %>% 
  pivot_longer(acc:spec, names_to = "classperf", values_to = "val") %>% 
  mutate(classperf = ifelse(classperf == "acc", "accuracy", 
                            ifelse(classperf == "sens", "sensitivity", "specificity")),
         modtype = gsub("Banded", "band_", modtype)) %>% 
  filter(classperf == "sensitivity") %>% 
  ggplot(aes(x = factor(k), y = val, fill = I("skyblue2"))) +
  geom_boxplot(position = "dodge") +
  facet_grid(rows = vars(missing), cols = vars(modtype), switch = "y") +
  # scale_fill_manual(values = c("tomato", "skyblue2", "darkgreen", "gold")) +
  labs(x = "p", fill = "", y = "sensitivity") +
  theme_bw() +
  theme(legend.position = "top", text = element_text(size = 16))
dev.off()
# TABLE
class_index_all %>% 
  pivot_longer(acc:spec, names_to = "classperf", values_to = "val") %>% 
  group_by(k, missing, model, modtype) %>% 
  skimr::skim(val)

#### KL-discrepancy per lambda ####
# PLOT
pdf("Images/Sim_poisson_banded_k_miss_sparse_KL.pdf", width = 12, height = 7)
KL_all %>% mutate(modtype = gsub("Banded", "band_", modtype)) %>% 
  # filter(SB<50) %>% # decommentare se Poisson
  ggplot(aes(x = factor(k), y = log(SB), fill = I("skyblue2"))) +
  geom_boxplot(position = "dodge") +
  facet_grid(rows = vars(missing), cols = vars(modtype), switch = "y") +
  # scale_fill_manual(values = c("tomato", "skyblue2", "darkgreen", "gold")) +
  labs(x = "p", fill = "", y = "KL-discrepancy (log)") +
  theme_bw() +
  theme(legend.position = "top", text = element_text(size = 16))
dev.off()
# TABLE
KL_all %>% 
  group_by(k, missing, modtype) %>% 
  skimr::skim(SB)

#### CRPS ####
# PLOT
crps_all %>% 
  mutate(modtype = ifelse(modtype == "Random1", "0.0", 
                          ifelse(modtype == "Random2", "0.25", 
                                 ifelse(modtype == "Random3", "0.5", "0.75")))) %>% 
  filter(missing == "0.1") %>% 
  ggplot(aes(x = factor(k), y = SB, fill = modtype)) +
  geom_boxplot(position = "dodge") +
  scale_fill_manual(values = c("tomato", "skyblue2", "darkgreen", "gold")) +
  labs(x = "k", fill = "", y = "CRPS") +
  theme_bw() +
  theme(legend.position = "top", text = element_text(size = 16))
dev.ff()
# TABLE
crps_all %>% 
  filter(missing == "0.1") %>% 
  group_by(k, missing, modtype) %>% 
  skimr::skim(SB)

