#### Load libraries ####

#Define R library path for external libraries (not needed on most computers, but necessary on CSC Puhti)
.libPaths(c("/projappl/project_2006897/project_rpackages_4.4.0", .libPaths()))
libpath <- .libPaths()[1]

#Load R packages used in this script
library("tidyverse") 
library("data.table")
library("janitor")
library("meta")
library("ggpubr")
library("svglite")
library("corrplot")
library("skiftiTools")


#### Data ####

options(scipen = 10000)

#Load data and filter out the 11,000 voxels where Generation R == 0; pivot longer for density plot

#Individual site effects
FA_cope <- fread("data/intercept/FA_cope_positive_FA_Skiftidata.txt", skip = 8) %>% 
  t() %>% as.data.frame() %>% 
  row_to_names(1) %>% 
  mutate_all(as.numeric) %>% 
  filter(genr != 0) %>% select(-preobe) %>% 
  pivot_longer(everything(), names_to = "site") %>% 
  mutate(site = factor(site, levels = c("fb1mo", "fb5y", "copsych", "genr", "nfbc")))
FA_sex_cope <- fread("data/intercept/FA_sex_cope_positive_FA_Skiftidata.txt", skip = 8) %>% 
  t() %>% as.data.frame() %>% 
  row_to_names(1) %>% 
  mutate_all(as.numeric) %>% 
  filter(genr != 0) %>% select(-preobe) %>%
  pivot_longer(everything(), names_to = "site") %>% 
  mutate(site = factor(site, levels = c("fb1mo", "fb5y", "copsych", "genr", "nfbc")))
MD_cope <- fread("data/intercept/MD_cope_positive_MD_Skiftidata.txt", skip = 8) %>% 
  t() %>% as.data.frame() %>% 
  row_to_names(1) %>% 
  mutate_all(as.numeric) %>% 
  filter(genr != 0) %>% select(-preobe) %>%
  pivot_longer(everything(), names_to = "site") %>% 
  mutate(site = factor(site, levels = c("fb1mo", "fb5y", "copsych", "genr", "nfbc")))
MD_sex_cope <- fread("data/intercept/MD_sex_cope_positive_MD_Skiftidata.txt", skip = 8) %>% 
  t() %>% as.data.frame() %>% 
  row_to_names(1) %>% 
  mutate_all(as.numeric) %>% 
  filter(genr != 0) %>% select(-preobe) %>%
  pivot_longer(everything(), names_to = "site") %>% 
  mutate(site = factor(site, levels = c("fb1mo", "fb5y", "copsych", "genr", "nfbc")))

#Individual site errors
FA_varcope <- fread("data/intercept/FA_varcope_positive_FA_Skiftidata.txt", skip = 8) %>% 
  t() %>% as.data.frame() %>% 
  row_to_names(1) %>% 
  mutate_all(as.numeric) %>% 
  filter(genr != 0) %>% select(-preobe) %>%
  pivot_longer(everything(), names_to = "site") %>% 
  mutate(site = factor(site, levels = c("fb1mo", "fb5y", "copsych", "genr", "nfbc"))) %>% 
  mutate(value = sqrt(value))
FA_sex_varcope <- fread("data/intercept/FA_sex_varcope_positive_FA_Skiftidata.txt", skip = 8) %>% 
  t() %>% as.data.frame() %>% 
  row_to_names(1) %>% 
  mutate_all(as.numeric) %>% 
  filter(genr != 0) %>% select(-preobe) %>%
  pivot_longer(everything(), names_to = "site") %>% 
  mutate(site = factor(site, levels = c("fb1mo", "fb5y", "copsych", "genr", "nfbc"))) %>% 
  mutate(value = sqrt(value))
MD_varcope <- fread("data/intercept/MD_varcope_positive_MD_Skiftidata.txt", skip = 8) %>% 
  t() %>% as.data.frame() %>% 
  row_to_names(1) %>% 
  mutate_all(as.numeric) %>% 
  filter(genr != 0) %>% select(-preobe) %>%
  pivot_longer(everything(), names_to = "site") %>% 
  mutate(site = factor(site, levels = c("fb1mo", "fb5y", "copsych", "genr", "nfbc"))) %>% 
  mutate(value = sqrt(value))
MD_sex_varcope <- fread("data/intercept/MD_sex_varcope_positive_MD_Skiftidata.txt", skip = 8) %>% 
  t() %>% as.data.frame() %>% 
  row_to_names(1) %>% 
  mutate_all(as.numeric) %>% 
  filter(genr != 0) %>% select(-preobe) %>%
  pivot_longer(everything(), names_to = "site") %>% 
  mutate(site = factor(site, levels = c("fb1mo", "fb5y", "copsych", "genr", "nfbc"))) %>% 
  mutate(value = sqrt(value))

#Meta-analysis results
FA_metagen_results <- readRDS("results/FA_positive_metagen_results_no_preobe.rds")
FA_sex_metagen_results <- readRDS("results/FA_sex_positive_metagen_results_no_preobe.rds")
MD_metagen_results <- readRDS("results/MD_positive_metagen_results_no_preobe.rds")
MD_sex_metagen_results <- readRDS("results/MD_sex_positive_metagen_results_no_preobe.rds")


#### Approximate standardised effect sizes ####

#SDx
SD_logCRP <- 0.8554
SD_logCRPsex <- 1.5485
SD_FAmean <- 0.06761415
SD_MDmean <- 0.00007255415

#SDy
SD_FA <- fread("data/FA_SD_FA_Skiftidata.txt", skip = 8) %>% t() %>% as.data.frame() %>% row_to_names(1) %>% 
  rownames_to_column("voxel") %>% rename("SDy" = "fb") %>% mutate(SDy = as.numeric(SDy))
SD_MD <- fread("data/MD_SD_MD_Skiftidata.txt", skip = 8) %>% t() %>% as.data.frame() %>% row_to_names(1) %>% 
  rownames_to_column("voxel") %>% rename("SDy" = "fb") %>% mutate(SDy = as.numeric(SDy))

#Create new variables
FA_metagen_results <- FA_metagen_results %>% 
  left_join(., SD_FA, by = "voxel") %>% 
  mutate(effect_standard = effect*SD_logCRP/SDy)

FA_sex_metagen_results <- FA_sex_metagen_results %>% 
  left_join(., SD_FA, by = "voxel") %>% 
  mutate(effect_standard = effect*SD_logCRPsex/SDy)

MD_metagen_results <- MD_metagen_results %>% 
  left_join(., SD_MD, by = "voxel") %>% 
  mutate(effect_standard = effect*SD_logCRP/SDy)

MD_sex_metagen_results <- MD_sex_metagen_results %>% 
  left_join(., SD_MD, by = "voxel") %>% 
  mutate(effect_standard = effect*SD_logCRPsex/SDy)

#Individual site effects
FA_function <- function(x, SDy_col) {x * SD_logCRP / SDy_col}

FA_sex_function <- function(x, SDy_col) {x * SD_logCRPsex / SDy_col}

FA_cope_standard <- fread("data/intercept/FA_cope_positive_FA_Skiftidata.txt", skip = 8) %>%
  t() %>% as.data.frame() %>% 
  row_to_names(1) %>% 
  mutate_all(as.numeric) %>% 
  filter(genr != 0) %>% 
  select(-preobe) %>%
  rownames_to_column("voxel") %>% 
  left_join(., SD_FA, by = "voxel") %>% 
  mutate(across(2:7, ~ FA_function(.x, SDy))) %>% 
  select(-voxel, -SDy) %>% 
  pivot_longer(everything(), names_to = "site") %>% 
  mutate(site = factor(site, levels = c("fb1mo", "fb5y", "copsych", "genr", "nfbc")))

FA_sex_cope_standard  <- fread("data/intercept/FA_sex_cope_positive_FA_Skiftidata.txt", skip = 8) %>% 
  t() %>% as.data.frame() %>% 
  row_to_names(1) %>% 
  mutate_all(as.numeric) %>% 
  filter(genr != 0) %>% 
  select(-preobe) %>%
  rownames_to_column("voxel") %>% 
  left_join(., SD_FA, by = "voxel") %>% 
  mutate(across(2:7, ~ FA_sex_function(.x, SDy))) %>% 
  select(-voxel, -SDy) %>% 
  pivot_longer(everything(), names_to = "site") %>% 
  mutate(site = factor(site, levels = c("fb1mo", "fb5y", "copsych", "genr", "nfbc")))

MD_cope_standard  <- fread("data/intercept/MD_cope_positive_MD_Skiftidata.txt", skip = 8) %>% 
  t() %>% as.data.frame() %>% 
  row_to_names(1) %>% 
  mutate_all(as.numeric) %>% 
  filter(genr != 0) %>% 
  select(-preobe) %>%
  rownames_to_column("voxel") %>% 
  left_join(., SD_MD, by = "voxel") %>% 
  mutate(across(2:7, ~ FA_function(.x, SDy))) %>% 
  select(-voxel, -SDy) %>% 
  pivot_longer(everything(), names_to = "site") %>% 
  mutate(site = factor(site, levels = c("fb1mo", "fb5y", "copsych", "genr", "nfbc")))

MD_sex_cope_standard  <- fread("data/intercept/MD_sex_cope_positive_MD_Skiftidata.txt", skip = 8) %>% 
  t() %>% as.data.frame() %>% 
  row_to_names(1) %>% 
  mutate_all(as.numeric) %>% 
  filter(genr != 0) %>% 
  select(-preobe) %>%
  rownames_to_column("voxel") %>% 
  left_join(., SD_MD, by = "voxel") %>% 
  mutate(across(2:7, ~ FA_sex_function(.x, SDy))) %>% 
  select(-voxel, -SDy) %>% 
  pivot_longer(everything(), names_to = "site") %>% 
  mutate(site = factor(site, levels = c("fb1mo", "fb5y", "copsych", "genr", "nfbc")))

#Individual site errors

FA_varcope_standard <- fread("data/intercept/FA_varcope_positive_FA_Skiftidata.txt", skip = 8) %>%
  t() %>% as.data.frame() %>% 
  row_to_names(1) %>% 
  mutate_all(as.numeric) %>% 
  filter(genr != 0) %>% 
  select(-preobe) %>%
  rownames_to_column("voxel") %>% 
  left_join(., SD_FA, by = "voxel") %>% 
  mutate(across(2:7, ~ FA_function(.x, SDy))) %>% 
  select(-voxel, -SDy) %>% 
  pivot_longer(everything(), names_to = "site") %>% 
  mutate(site = factor(site, levels = c("fb1mo", "fb5y", "copsych", "genr", "nfbc")))

FA_sex_varcope_standard  <- fread("data/intercept/FA_sex_varcope_positive_FA_Skiftidata.txt", skip = 8) %>% 
  t() %>% as.data.frame() %>% 
  row_to_names(1) %>% 
  mutate_all(as.numeric) %>% 
  filter(genr != 0) %>% 
  select(-preobe) %>%
  rownames_to_column("voxel") %>% 
  left_join(., SD_FA, by = "voxel") %>% 
  mutate(across(2:7, ~ FA_sex_function(.x, SDy))) %>% 
  select(-voxel, -SDy) %>% 
  pivot_longer(everything(), names_to = "site") %>% 
  mutate(site = factor(site, levels = c("fb1mo", "fb5y",  "copsych", "genr", "nfbc")))

MD_varcope_standard  <- fread("data/intercept/MD_varcope_positive_MD_Skiftidata.txt", skip = 8) %>% 
  t() %>% as.data.frame() %>% 
  row_to_names(1) %>% 
  mutate_all(as.numeric) %>% 
  filter(genr != 0) %>% 
  select(-preobe) %>%
  rownames_to_column("voxel") %>% 
  left_join(., SD_MD, by = "voxel") %>% 
  mutate(across(2:7, ~ FA_function(.x, SDy))) %>% 
  select(-voxel, -SDy) %>% 
  pivot_longer(everything(), names_to = "site") %>% 
  mutate(site = factor(site, levels = c("fb1mo", "fb5y", "copsych", "genr", "nfbc")))

MD_sex_varcope_standard  <- fread("data/intercept/MD_sex_varcope_positive_MD_Skiftidata.txt", skip = 8) %>% 
  t() %>% as.data.frame() %>% 
  row_to_names(1) %>% 
  mutate_all(as.numeric) %>% 
  filter(genr != 0) %>% 
  select(-preobe) %>%
  rownames_to_column("voxel") %>% 
  left_join(., SD_MD, by = "voxel") %>% 
  mutate(across(2:7, ~ FA_sex_function(.x, SDy))) %>% 
  select(-voxel, -SDy) %>% 
  pivot_longer(everything(), names_to = "site") %>% 
  mutate(site = factor(site, levels = c("fb1mo", "fb5y", "copsych", "genr", "nfbc")))

#### Plot distribution of effect sizes and errors in each site ####

#FA
FA_density <- ggplot(FA_cope_standard) +
  geom_density(aes(x=value, fill = site, colour = site), alpha = 0.3, linewidth = 1) +
  scale_fill_manual(values = c("#FF6467", "#7C86FF",  "#5EE9B5" , "#74D4FF", "#FE9A37")) +
  scale_colour_manual(values = c("#FF6467", "#7C86FF",  "#5EE9B5" , "#74D4FF", "#FE9A37")) +
  labs(x = "Covariate-adjusted coefficients", y = "Density") +
  xlim(-1, 1) +
  theme(panel.background = element_blank(), 
        legend.position = "none",
        plot.background = element_blank(),
        legend.background = element_blank(),
        legend.title = element_blank())

FA_error_density <- ggplot(FA_varcope) +
  geom_density(aes(x=value, fill = site, colour = site), alpha = 0.3, linewidth = 1) +
  scale_fill_manual(values = c("#FF6467", "#7C86FF",  "#5EE9B5" , "#74D4FF", "#FE9A37")) +
  scale_colour_manual(values = c("#FF6467", "#7C86FF",  "#5EE9B5" , "#74D4FF", "#FE9A37")) +
  labs(y = "Density", x = "Standard errors") +
  xlim(-0.0025, 0.035) +
  theme(panel.background = element_blank(), 
        plot.background = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank())

site_weights <- data.frame(weights = colMeans(FA_metagen_results[, 6:10])) %>% 
  rownames_to_column("site") %>% 
  mutate(site = factor(site, levels = c("fb1mo", "fb5y", "copsych", "genr", "nfbc"))) %>% 
  mutate(percent = (weights/(sum(weights))*100)) %>% 
  mutate(ymax = cumsum(percent)) %>% 
  mutate(ymin = c(0, head(ymax, n=-1)))

FA_weights <- ggplot(data = site_weights, aes(ymax=ymax, ymin=ymin, xmax =8, xmin=7, fill= site)) +
  geom_rect() +
  scale_fill_manual(values = c("#FF6467", "#7C86FF", "#5EE9B5" , "#74D4FF", "#FE9A37")) +
  coord_polar(theta="y", start = 6) + 
  xlim(c(6, 8)) + 
  theme_void() +
  labs(title = "Mean site weighting") +
  theme(plot.title = element_text(size = 10, hjust = 0.5, vjust = -5), legend.position = "none")

ggarrange(FA_density, FA_error_density, ncol = 2)
ggsave("plots/FA_each_site_no_preobe.svg", units = "px", width = 3000, height = 1000)

#MD
MD_density <- ggplot(MD_cope_standard) +
  geom_density(aes(x=value, fill = site, colour = site), alpha = 0.3, linewidth = 1) +
  scale_fill_manual(values = c("#FF6467", "#7C86FF", "#5EE9B5" , "#74D4FF", "#FE9A37")) +
  scale_colour_manual(values = c("#FF6467", "#7C86FF", "#5EE9B5" , "#74D4FF", "#FE9A37")) +
  labs(x = "Covariate-adjusted coefficients", y = "Density") +
  xlim(-1, 1) +
  theme(panel.background = element_blank(), 
        legend.position = "none",
        axis.title.y = element_blank(),
        plot.background = element_blank(),
        legend.background = element_blank(),
        legend.title = element_blank())

MD_error_density <- ggplot(MD_varcope) +
  geom_density(aes(x=value, fill = site, colour = site), alpha = 0.3, linewidth = 1) +
  scale_fill_manual(values = c("#FF6467", "#7C86FF",  "#5EE9B5" , "#74D4FF", "#FE9A37")) +
  scale_colour_manual(values = c("#FF6467", "#7C86FF",  "#5EE9B5" , "#74D4FF", "#FE9A37")) +
  labs(x = "Standard errors") +
  scale_x_continuous(limits = c(-0.000005, 0.00005), breaks = c(0, 0.00002, 0.00004)) +
  theme(panel.background = element_blank(), 
        plot.background = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank(),
        axis.title.y = element_blank())

site_weights <- data.frame(weights = colMeans(MD_metagen_results[, 6:10])) %>% 
  rownames_to_column("site") %>% 
  mutate(site = factor(site, levels = c("fb1mo", "fb5y", "copsych", "genr", "nfbc"))) %>% 
  mutate(percent = (weights/(sum(weights))*100)) %>% 
  mutate(ymax = cumsum(percent)) %>% 
  mutate(ymin = c(0, head(ymax, n=-1)))

MD_weights <- ggplot(data = site_weights, aes(ymax=ymax, ymin=ymin, xmax =8, xmin=7, fill= site)) +
  geom_rect() +
  scale_fill_manual(values = c("#FF6467", "#7C86FF", "#5EE9B5" , "#74D4FF", "#FE9A37")) +
  coord_polar(theta="y", start = 6) + 
  xlim(c(6, 8)) + 
  theme_void() +
  labs(title = "Mean site weighting") +
  theme(plot.title = element_text(size = 10, hjust = 0.5, vjust = -5), legend.position = "none")

ggarrange(MD_density, MD_error_density, ncol = 2)
ggsave("plots/MD_each_site_no_preobe.svg", units = "px", width = 3000, height = 1000)


#Save combined FA and MD
ggarrange(FA_density, MD_density, ncol = 2)
ggsave("plots/effect_sizes_each_site_no_preobe.svg", units = "px", width = 1500, height = 1000)

ggarrange(FA_error_density, MD_error_density, ncol = 2)
ggsave("plots/errors_each_site_no_preobe.svg", units = "px", width = 1500, height = 1000)

#FA x sex
FA_density <- ggplot(FA_sex_cope_standard) +
  geom_density(aes(x=value, fill = site, colour = site), alpha = 0.3, linewidth = 1.5) +
  scale_fill_manual(values = c("#FF6467", "#7C86FF", "#5EE9B5" , "#74D4FF", "#FE9A37")) +
  scale_colour_manual(values = c("#FF6467", "#7C86FF", "#5EE9B5" , "#74D4FF", "#FE9A37")) +
  labs(x = "Covariate-adjusted coefficients", y = "Density") +
  xlim(-3, 3) +
  theme(panel.background = element_blank(), 
        plot.background = element_blank(),
        legend.background = element_blank(),
        legend.title = element_blank())

FA_error_density <- ggplot(FA_sex_varcope) +
  geom_density(aes(x=value, fill = site, colour = site), alpha = 0.3, linewidth = 1.5) +
  scale_fill_manual(values = c("#FF6467", "#7C86FF",  "#5EE9B5" , "#74D4FF", "#FE9A37")) +
  scale_colour_manual(values = c("#FF6467", "#7C86FF",  "#5EE9B5" , "#74D4FF", "#FE9A37")) +
  labs(x = "Standard errors") +
  xlim(-0.003, 0.05) +
  theme(panel.background = element_blank(), 
        plot.background = element_blank(),
        legend.position = c(0.8, 0.65),
        legend.title = element_blank(),
        legend.background = element_blank(),
        axis.title.y = element_blank())

site_weights <- data.frame(weights = colMeans(FA_sex_metagen_results[, 8:12])) %>% 
  rownames_to_column("site") %>% 
  mutate(site = factor(site, levels = c("fb1mo", "fb5y", "copsych", "genr", "nfbc"))) %>% 
  mutate(percent = (weights/(sum(weights))*100)) %>% 
  mutate(ymax = cumsum(percent)) %>% 
  mutate(ymin = c(0, head(ymax, n=-1)))

FA_weights <- ggplot(data = site_weights, aes(ymax=ymax, ymin=ymin, xmax =8, xmin=7, fill= site)) +
  geom_rect() +
  scale_fill_manual(values = c("#FF6467", "#7C86FF", "#5EE9B5" , "#74D4FF", "#FE9A37")) +
  coord_polar(theta="y", start = 6) + 
  xlim(c(6, 8)) + 
  theme_void() +
  labs(title = "Mean site weighting") +
  theme(plot.title = element_text(size = 10, hjust = 0.5, vjust = -5), legend.position = "none")

ggarrange(FA_density, FA_weights, ncol = 2, widths = c(2,1))
ggsave("plots/FA_sex_each_site_no_preobe.svg", units = "px", width = 3000, height = 1000)

#MD x sex
MD_density <- ggplot(MD_sex_cope_standard) +
  geom_density(aes(x=value, fill = site, colour = site), alpha = 0.3, linewidth = 1.5) +
  scale_fill_manual(values = c("#FF6467", "#7C86FF",  "#5EE9B5" , "#74D4FF", "#FE9A37")) +
  scale_colour_manual(values = c("#FF6467", "#7C86FF", "#5EE9B5" , "#74D4FF", "#FE9A37")) +
  labs(x = "Covariate-adjusted coefficients", y = "Density") +
  xlim(-3, 3) +
  theme(panel.background = element_blank(), 
        plot.background = element_blank(),
        legend.background = element_blank(),
        legend.title = element_blank())

MD_error_density <- ggplot(MD_sex_varcope) +
  geom_density(aes(x=value, fill = site, colour = site), alpha = 0.3, linewidth = 1.5) +
  scale_fill_manual(values = c("#FF6467", "#7C86FF",  "#5EE9B5" , "#74D4FF", "#FE9A37")) +
  scale_colour_manual(values = c("#FF6467", "#7C86FF",  "#5EE9B5" , "#74D4FF", "#FE9A37")) +
  labs(x = "Standard errors") +
  xlim(-0.000005, 0.000075) +
  theme(panel.background = element_blank(), 
        plot.background = element_blank(),
        legend.position = c(0.8, 0.65),
        legend.title = element_blank(),
        legend.background = element_blank(),
        axis.title.y = element_blank())

site_weights <- data.frame(weights = colMeans(MD_sex_metagen_results[, 8:12])) %>% 
  rownames_to_column("site") %>% 
  mutate(site = factor(site, levels = c("fb1mo", "fb5y", "copsych", "genr", "nfbc"))) %>% 
  mutate(percent = (weights/(sum(weights))*100)) %>% 
  mutate(ymax = cumsum(percent)) %>% 
  mutate(ymin = c(0, head(ymax, n=-1)))

MD_weights <- ggplot(data = site_weights, aes(ymax=ymax, ymin=ymin, xmax =8, xmin=7, fill= site)) +
  geom_rect() +
  scale_fill_manual(values = c("#FF6467", "#7C86FF", "#5EE9B5" , "#74D4FF", "#FE9A37")) +
  coord_polar(theta="y", start = 6) + 
  xlim(c(6, 8)) + 
  theme_void() +
  labs(title = "Mean site weighting") +
  theme(plot.title = element_text(size = 10, hjust = 0.5, vjust = -5), legend.position = "none")

ggarrange(MD_density, MD_weights, ncol = 2, widths = c(2,1))
ggsave("plots/MD_sex_each_site_no_preobe.svg", units = "px", width = 3000, height = 1000)


#### Plot voxel-wise correlations of effect sizes in each site ####

FA_cope <- fread("data/intercept/FA_cope_positive_FA_Skiftidata.txt", skip = 8) %>% 
  t() %>% as.data.frame() %>% 
  row_to_names(1) %>% 
  mutate_all(as.numeric) %>% 
  filter(genr != 0) %>% select(-preobe)

svglite("plots/FA_cope_corrplot_no_preobe.svg", width=5, height=5, bg = "transparent")
corrplot(corr = cor(FA_cope, method = 'spearman'), 
         type = "lower", method = "color", number.cex = 0.2, tl.col = "black", bg = "transparent")
dev.off()

FA_sex_cope <- fread("data/intercept/FA_sex_cope_positive_FA_Skiftidata.txt", skip = 8) %>% 
  t() %>% as.data.frame() %>% 
  row_to_names(1) %>% 
  mutate_all(as.numeric) %>% 
  filter(genr != 0) %>% select(-preobe)

svglite("plots/FA_sex_cope_corrplot_no_preobe.svg", width=5, height=5, bg = "transparent")
corrplot(corr = cor(FA_sex_cope, method = 'spearman'), 
         type = "lower", method = "color", number.cex = 0.2, tl.col = "black", bg = "transparent")
dev.off()

MD_cope <- fread("data/intercept/MD_cope_positive_MD_Skiftidata.txt", skip = 8) %>% 
  t() %>% as.data.frame() %>% 
  row_to_names(1) %>% 
  mutate_all(as.numeric) %>% 
  filter(genr != 0) %>% select(-preobe)

svglite("plots/MD_cope_corrplot_no_preobe.svg", width=5, height=5, bg = "transparent")
corrplot(corr = cor(MD_cope, method = 'spearman'), 
         type = "lower", method = "color", number.cex = 0.2, tl.col = "black", bg = "transparent")
dev.off()

MD_sex_cope <- fread("data/intercept/MD_sex_cope_positive_MD_Skiftidata.txt", skip = 8) %>% 
  t() %>% as.data.frame() %>% 
  row_to_names(1) %>% 
  mutate_all(as.numeric) %>% 
  filter(genr != 0) %>% select(-preobe)

svglite("plots/MD_sex_cope_corrplot_no_preobe.svg", width=5, height=5, bg = "transparent")
corrplot(corr = cor(MD_sex_cope, method = 'spearman'), 
         type = "lower", method = "color", number.cex = 0.2, tl.col = "black", bg = "transparent")
dev.off()

#### Voxel-level p-value corrections ####
#Galwey
FA_galwey <- readRDS("results/FA_galwey.rds")
MD_galwey <- readRDS("results/MD_galwey.rds")

#FDR
FA_metagen_results <- FA_metagen_results %>% filter(genr != 0) %>% 
  arrange(pval) %>% 
  mutate(pval_rank = 1:nrow(.)) %>% 
  mutate(p_FDR = if_else(pval*nrow(.)/pval_rank > 1, 1.0, pval*nrow(.)/pval_rank)) %>% 
  mutate(meff = if_else(pval*FA_galwey < 1, pval*FA_galwey, 1.0)) %>% 
  arrange(voxel)

FA_sex_metagen_results <- FA_sex_metagen_results %>% filter(genr != 0) %>% 
  arrange(pval) %>% 
  mutate(pval_rank = 1:nrow(.)) %>% 
  mutate(p_FDR = if_else(pval*nrow(.)/pval_rank > 1, 1.0, pval*nrow(.)/pval_rank)) %>% 
  mutate(meff = if_else(pval*FA_galwey < 1, pval*FA_galwey, 1.0)) %>% 
  arrange(voxel)

MD_metagen_results <- MD_metagen_results %>% filter(genr != 0) %>% 
  arrange(pval) %>% 
  mutate(pval_rank = 1:nrow(.)) %>% 
  mutate(p_FDR = if_else(pval*nrow(.)/pval_rank > 1, 1.0, pval*nrow(.)/pval_rank)) %>% 
  mutate(meff = if_else(pval*MD_galwey < 1, pval*MD_galwey, 1.0)) %>% 
  arrange(voxel)

MD_sex_metagen_results <- MD_sex_metagen_results %>% filter(genr != 0) %>% 
  arrange(pval) %>% 
  mutate(pval_rank = 1:nrow(.)) %>% 
  mutate(p_FDR = if_else(pval*nrow(.)/pval_rank > 1, 1.0, pval*nrow(.)/pval_rank)) %>% 
  mutate(meff = if_else(pval*MD_galwey < 1, pval*MD_galwey, 1.0)) %>% 
  arrange(voxel)


#### Pooled effect volcano plots ####

#FA pooled
FA_metagen_results %>%  
  
  mutate(differences = case_when(
    pval < 0.05 & effect > 0 ~ "Nominal increase", 
    pval < 0.05 & effect< 0 ~ "Nominal decrease",
    TRUE ~ "Spurious")) %>% 
  
  ggplot(data = ., aes(x = effect_standard, y = -log10(pval), colour = differences, size = differences)) + 
  geom_point(position = "jitter") +
  xlab("Pooled effect") +
  ylab("-log10(p-value)") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", colour = "black", alpha = 0.75, size = 0.5) +
  geom_hline(yintercept = -log10(0.05/FA_galwey), linetype="dashed", colour = "black", alpha = 0.75, size = 0.5) +
  annotate("text", label = "p < 0.05", x=0, y=-log10(0.06), size = 4) +
  annotate("text", label = "padj < 0.05", x=0, y=-log10(0.06/FA_galwey), size = 4) +
  scale_colour_manual(values = c("#619CFF", "#F8766D","lightgrey")) +
  scale_size_manual(values = c(1, 1, 0.1), guide = "none") + 
  xlim(-0.5, 0.5) +
  ylim(0, 4.5) +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14), 
        legend.position = "none")

ggsave("plots/FA_volcano_no_preobe.svg", units = "px", width = 1500, height = 1500)
ggsave("plots/FA_volcano_no_preobe.png", units = "px", width = 1500, height = 1500)


#MD pooled
MD_metagen_results %>%  
  
  mutate(differences = case_when(
    pval < 0.05 & effect > 0 ~ "Nominal increase", 
    pval < 0.05 & effect < 0 ~ "Nominal decrease",
    TRUE ~ "Spurious")) %>% 
  
  ggplot(data = ., aes(x = effect_standard, y = -log10(pval), colour = differences, size = differences)) + 
  geom_point(position = "jitter") +
  xlab("Pooled effect") +
  ylab("-log10(p-value)") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", colour = "black", alpha = 0.75, size = 0.5) +
  geom_hline(yintercept = -log10(0.05/MD_galwey), linetype="dashed", colour = "black", alpha = 0.75, size = 0.5) +
  annotate("text", label = "p < 0.05", x=0, y=-log10(0.06), size = 4) +
  annotate("text", label = "padj < 0.05", x=0, y=-log10(0.06/MD_galwey), size = 4) +
  scale_colour_manual(values = c("#619CFF", "#F8766D","lightgrey")) +
  scale_size_manual(values = c(1, 1, 0.1), guide = "none") + 
  xlim(-0.5, 0.5) +
  ylim(0, 4.5) +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14), 
        legend.position = "none")

ggsave("plots/MD_volcano_no_preobe.svg", units = "px", width = 1500, height = 1500)
ggsave("plots/MD_volcano_no_preobe.png", units = "px", width = 1500, height = 1500)


#FA x sex pooled
FA_sex_metagen_results %>%  
  
  mutate(differences = case_when(
    pval < 0.05 & effect > 0 ~ "Nominal increase", 
    pval < 0.05 & effect < 0 ~ "Nominal decrease",
    TRUE ~ "Spurious")) %>% 
  
  ggplot(data = ., aes(x = effect_standard, y = -log10(pval), colour = differences, size = differences)) + 
  geom_point(position = "jitter") +
  xlab("Pooled effect") +
  ylab("-log10(p-value)") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", colour = "black", alpha = 0.75, size = 0.5) +
  geom_hline(yintercept = -log10(0.05/FA_galwey), linetype="dashed", colour = "black", alpha = 0.75, size = 0.5) +
  annotate("text", label = "p < 0.05", x=0, y=-log10(0.06), size = 4) +
  annotate("text", label = "padj < 0.05", x=0, y=-log10(0.06/FA_galwey), size = 4) +
  scale_colour_manual(values = c("#619CFF", "#F8766D","lightgrey")) +
  scale_size_manual(values = c(1, 1, 0.1), guide = "none") + 
  xlim(-1, 1) +
  ylim(0, 5) +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14), 
        legend.position = "none")

ggsave("plots/FA_sex_volcano_no_preobe.svg", units = "px", width = 1500, height = 1500)
ggsave("plots/FA_sex_volcano_no_preobe.png", units = "px", width = 1500, height = 1500)


#MD x sex pooled
MD_sex_metagen_results %>%  
  
  mutate(differences = case_when(
    pval < 0.05 & effect > 0 ~ "Nominal increase", 
    pval < 0.05 & effect < 0 ~ "Nominal decrease",
    TRUE ~ "Spurious")) %>% 
  
  ggplot(data = ., aes(x = effect_standard, y = -log10(pval), colour = differences, size = differences)) + 
  geom_point(position = "jitter") +
  xlab("Pooled effect") +
  ylab("-log10(p-value)") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", colour = "black", alpha = 0.75, size = 0.5) +
  geom_hline(yintercept = -log10(0.05/MD_galwey), linetype="dashed", colour = "black", alpha = 0.75, size = 0.5) +
  annotate("text", label = "p < 0.05", x=0, y=-log10(0.075), size = 4) +
  annotate("text", label = "padj < 0.05", x=0, y=-log10(0.075/MD_galwey), size = 4) +
  scale_colour_manual(values = c("#619CFF", "#F8766D","lightgrey")) +
  scale_size_manual(values = c(1, 1, 0.1), guide = "none") + 
  xlim(-1, 1) +
  ylim(0, 5) +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14), 
        legend.position = "none")

ggsave("plots/MD_sex_volcano_no_preobe.svg", units = "px", width = 1500, height = 1500)
ggsave("plots/MD_sex_volcano_no_preobe.png", units = "px", width = 1500, height = 1500)



#### Save skeletons ####

source("scripts/Skifti2Nifti.R")

tabular_results <- data.frame(FA_effect = FA_metagen_results$effect,
                              FA_effect_std = FA_metagen_results$effect_standard,
                              FA_pval = FA_metagen_results$pval,
                              FA_1minus_pval = 1-FA_metagen_results$pval,
                              FA_log_pval = -log10(FA_metagen_results$pval),
                              MD_effect = MD_metagen_results$effect,
                              MD_effect_std = MD_metagen_results$effect_standard,
                              MD_pval = MD_metagen_results$pval,
                              MD_1minus_pval = 1-MD_metagen_results$pval,
                              MD_log_pval = -log10(MD_metagen_results$pval))

Skifti_results <- list(reftype = "filename", refdata = "data/ENIGMA_TBSS_skeleton_mask.nii.gz", 
                       dim = c(3, 182, 218, 182, 1, 1, 1, 1), pixdim = c(-1, 1, 1, 1, 0, 0, 0, 0),
                       xform = matrix(c(-1, 0, 0, 90, 0, 1, 0, -126, 0, 0, 1, -72), nrow = 3, ncol = 4, byrow = TRUE),
                       data = t(tabular_results))

NIfTI_results <- skiftiTools::Skifti2Nifti(Skifti_results)


writeNifti(NIfTI_results[[1]], "results/FA_effect_no_preobe.nii.gz")
writeNifti(NIfTI_results[[2]], "results/FA_effect_standard_no_preobe.nii.gz")
writeNifti(NIfTI_results[[3]], "results/FA_pvalue_no_preobe.nii.gz")
writeNifti(NIfTI_results[[4]], "results/FA_1minus_pvalue_no_preobe.nii.gz")
writeNifti(NIfTI_results[[5]], "results/FA_minus_log10_pvalue_no_preobe.nii.gz")

writeNifti(NIfTI_results[[6]], "results/MD_effect_no_preobe.nii.gz")
writeNifti(NIfTI_results[[7]], "results/MD_effect_standard_no_preobe.nii.gz")
writeNifti(NIfTI_results[[8]], "results/MD_pvalue_no_preobe.nii.gz")
writeNifti(NIfTI_results[[9]], "results/MD_1minus_pvalue_no_preobe.nii.gz")
writeNifti(NIfTI_results[[10]], "results/MD_minus_log10_pvalue_no_preobe.nii.gz")


#For sex-interaction models

tabular_results <- data.frame(FA_effect_std = FA_sex_metagen_results$effect_standard,
                              FA_1minus_pval = 1-FA_sex_metagen_results$pval,
                              MD_effect_std = MD_sex_metagen_results$effect_standard,
                              MD_1minus_pval = 1-MD_sex_metagen_results$pval)

Skifti_results <- list(reftype = "filename", refdata = "data/ENIGMA_TBSS_skeleton_mask.nii.gz", 
                       dim = c(3, 182, 218, 182, 1, 1, 1, 1), pixdim = c(-1, 1, 1, 1, 0, 0, 0, 0),
                       xform = matrix(c(-1, 0, 0, 90, 0, 1, 0, -126, 0, 0, 1, -72), nrow = 3, ncol = 4, byrow = TRUE),
                       data = t(tabular_results))

NIfTI_results <- skiftiTools::Skifti2Nifti(Skifti_results)


writeNifti(NIfTI_results[[1]], "results/FA_sex_effect_standard_no_preobe.nii.gz")
writeNifti(NIfTI_results[[2]], "results/FA_sex_1minus_pvalue_no_preobe.nii.gz")
writeNifti(NIfTI_results[[3]], "results/MD_sex_effect_standard_no_preobe.nii.gz")
writeNifti(NIfTI_results[[4]], "results/MD_sex_1minus_pvalue_no_preobe.nii.gz")

#### Density plots for meta-regression ####

FA_metagen_results <- FA_metagen_results %>% 
  mutate(age_coef_standard = age_coef/SDy) %>% 
  mutate(measure = "FA")

MD_metagen_results <- MD_metagen_results %>% 
  mutate(age_coef_standard = age_coef/SDy) %>% 
  mutate(measure = "MD")

all_metagen_results <- rbind(FA_metagen_results, MD_metagen_results)

#Standardised coefficient
ggplot(all_metagen_results) +
  geom_density(aes(x=age_coef_standard, colour = measure, fill = measure), alpha = 0.3, linewidth = 0.75) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1, colour = "darkgrey") +
  scale_fill_manual(values = c("#BA4A00", "#619CFF")) +
  scale_colour_manual(values = c("#BA4A00", "#619CFF")) +
  labs(x = "Standardised β coefficient", y = "Density") +
  xlim(-0.05, 0.05) +
  theme(panel.background = element_blank(), 
        legend.position = "none",
        plot.background = element_blank(),
        legend.background = element_blank(),
        legend.title = element_blank())

ggsave("plots/meta_regression_effect_sizes.svg", units = "px", width = 1500, height = 1000)

#-lo10(pvalue)
ggplot(all_metagen_results) +
  geom_density(aes(x=-log10(age_pval), colour = measure, fill = measure), alpha = 0.3, linewidth = 0.75) +
  geom_vline(xintercept = -log10(0.05/MD_galwey), linetype = "dashed", size = 1, colour = "darkgrey") +
  
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", size = 1, colour = "darkgrey") +
  scale_fill_manual(values = c("#BA4A00", "#619CFF")) +
  scale_colour_manual(values = c("#BA4A00", "#619CFF")) +
  labs(x = "-log10(p-value)", y = "Density") +
  xlim(-0.075, 3.5) +
  theme(panel.background = element_blank(), 
        #legend.position = "none",
        plot.background = element_blank(),
        legend.background = element_blank(),
        legend.title = element_blank())

ggsave("plots/meta_regression_pvalues.svg", units = "px", width = 1500, height = 1000)
