# Carico i pacchetti
require(tidyverse)
require(magrittr)

# Load RData
load("WSReal/GeneExp100_Confronti_GW.Rdata") # cambiare l'RData


# Quanta e quale sparsità?
corrplot::corrplot(meanzeta, is.corr = F, order = "hclust")
sum(round(meanzeta) == 0) / prod(dim(meanzeta))

# Se Confronti
sum(round(apply(mcmc_gwish$zmat, c(1, 2), mean)) == 0) / prod(dim(meanzeta))

# Lambda est
meanlambda %>%
  as.data.frame() %>%
  rownames_to_column(var = "Rows") %>%
  pivot_longer(-Rows, names_to = "Cols", values_to = "val") %>%
  ggplot(aes(Rows, Cols, fill = val)) +
  geom_raster() +
  scale_fill_distiller(palette = "RdBu")

# Errre di stima (no per Poisson)
yest <- colMeans(missing_out)
sqrt(mean((true_y - yest)^2))
plot(true_y, yest, pch = 20)
abline(0, 1, lwd = 1, col = "red2")

# Se Confronti
summary(crps_gw)
summary(crps_sw)

tibble(gw = n_edge_gw, sw = n_edge_sw) %>%
  ggplot(aes(x = 1:2000)) +
  geom_line(aes(y = gw)) +
  geom_line(aes(y = sw), color = I("red3")) +
  labs(x = "iter", y = "n. edges") +
  theme_bw()
