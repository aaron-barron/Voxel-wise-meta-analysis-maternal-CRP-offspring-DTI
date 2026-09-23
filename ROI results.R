#### Libraries ####

#Define R library path for external libraries (not needed on most computers, but necessary on CSC Puhti)
.libPaths(c("/projappl/project_2006897/project_rpackages_4.4.0", .libPaths()))
libpath <- .libPaths()[1]

#Load R packages used in this script
library("tidyverse") 
library("data.table")
library("janitor")
library("svglite")


#### Data ####

#Standardise effect sizes
SD_logCRP <- 0.8554
SD_logCRPsex <- 1.5485
SD_FAmean <- 0.06761415
SD_MDmean <- 0.00007255415

FA_results <- readRDS("results/FA_ROI_results.rds") %>% mutate(effect_standard = effect*SD_logCRP/SD_FAmean) %>% 
  mutate(lower_standard = lower_effect*SD_logCRP/SD_FAmean) %>% mutate(upper_standard = upper_effect*SD_logCRP/SD_FAmean)

MD_results <- readRDS("results/MD_ROI_results.rds") %>% mutate(effect_standard = effect*SD_logCRP/SD_MDmean) %>% 
  mutate(lower_standard = lower_effect*SD_logCRP/SD_MDmean) %>% mutate(upper_standard = upper_effect*SD_logCRP/SD_MDmean)

FA_sex_results <- readRDS("results/FA_sex_ROI_results.rds") %>% mutate(effect_standard = effect*SD_logCRPsex/SD_FAmean) %>% 
  mutate(lower_standard = lower_effect*SD_logCRPsex/SD_FAmean) %>% mutate(upper_standard = upper_effect*SD_logCRPsex/SD_FAmean)

MD_sex_results <- readRDS("results/MD_sex_ROI_results.rds") %>% mutate(effect_standard = effect*SD_logCRPsex/SD_MDmean) %>% 
  mutate(lower_standard = lower_effect*SD_logCRPsex/SD_MDmean) %>% mutate(upper_standard = upper_effect*SD_logCRPsex/SD_MDmean)

FA_results_no_preobe <- readRDS("results/FA_ROI_no_preobe_results.rds") %>% mutate(effect_standard = effect*SD_logCRP/SD_FAmean) %>% 
  mutate(lower_standard = lower_effect*SD_logCRP/SD_FAmean) %>% mutate(upper_standard = upper_effect*SD_logCRP/SD_FAmean)

MD_results_no_preobe <- readRDS("results/MD_ROI_no_preobe_results.rds") %>% mutate(effect_standard = effect*SD_logCRP/SD_MDmean) %>% 
  mutate(lower_standard = lower_effect*SD_logCRP/SD_MDmean) %>% mutate(upper_standard = upper_effect*SD_logCRP/SD_MDmean)

FA_sex_results_no_preobe <- readRDS("results/FA_sex_ROI_no_preobe_results.rds") %>% mutate(effect_standard = effect*SD_logCRPsex/SD_FAmean) %>% 
  mutate(lower_standard = lower_effect*SD_logCRPsex/SD_FAmean) %>% mutate(upper_standard = upper_effect*SD_logCRPsex/SD_FAmean)

MD_sex_results_no_preobe <- readRDS("results/MD_sex_ROI_no_preobe_results.rds") %>% mutate(effect_standard = effect*SD_logCRPsex/SD_MDmean) %>% 
  mutate(lower_standard = lower_effect*SD_logCRP/SD_MDmean) %>% mutate(upper_standard = upper_effect*SD_logCRP/SD_MDmean)


#### Galwey meff estimation ####



#### Forest plots faceted per ROI ####

## FA ##
FA_results_forest <- FA_results %>% 
  rename(pooled_lower = lower_effect, pooled_upper = upper_effect) %>% 
  pivot_longer(cols = matches("_(lower|upper)$"), names_to = c("cohort", "bound"),
               names_pattern = "(.*)_(lower|upper)$", values_to = "value") %>% 
  pivot_wider(names_from = bound, values_from = value) %>%  
  mutate(mean = (lower + upper) / 2) %>% 
  mutate(lower = lower*SD_logCRP/SD_FAmean) %>% 
  mutate(upper = upper*SD_logCRP/SD_FAmean) %>% 
  mutate(mean = mean*SD_logCRP/SD_FAmean) %>% 
  mutate(cohort = factor(cohort, levels = c("fb1mo", "fb5y", "preobe", "copsych", "genr", "nfbc", "pooled")))

#Mean FA
ggplot(FA_results_forest %>% filter(ROI == "Mean") %>% filter(cohort != "pooled")) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", colour = "lightgrey") +
  geom_vline(xintercept = 6.5, size = 1, linetype = "dashed", colour = "lightgrey") +
  geom_point(aes(x = cohort, y = mean, colour = cohort), size = 2.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = cohort, colour = cohort), size = 1, width = 0.5, alpha = 0.75) +
  scale_colour_manual(values = c("#FF6467", "#7C86FF", "#E95FFC", "#5EE9B5" , "#74D4FF", "#FE9A37", "darkgrey")) +
  
  labs(y = "Standardised effect size") +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(size = 12, angle = 90),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        legend.position = "none")

ggsave("plots/rois/FA_forest_mean.svg", units = "px", width = 750, height = 750)


#FA in 48 ROIs, in 6x8 grid
ggplot(FA_results_forest %>% filter(ROI != "Mean")) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", colour = "lightgrey") +
  geom_vline(xintercept = 6.5, size = 1, linetype = "dashed", colour = "lightgrey") +
  geom_point(aes(x = cohort, y = mean, colour = cohort), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = cohort, colour = cohort), size = 1, width = 0.5, alpha = 0.75) +
  scale_colour_manual(values = c("#FF6467", "#7C86FF", "#E95FFC", "#5EE9B5" , "#74D4FF", "#FE9A37", "darkgrey")) +
  
  facet_wrap(~ROI, ncol = 6, nrow=8) +
  
  labs(y = "Standardised effect size") +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(size = 12, angle = 90),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        legend.position = "none")

ggsave("plots/rois/FA_forest_rois.svg", units = "px", width = 3000, height = 4000)

## MD ##
MD_results_forest <- MD_results %>% 
  rename(pooled_lower = lower_effect, pooled_upper = upper_effect) %>% 
  #mutate(preobe = preobe/10, preobe_lower = preobe_lower/10, preobe_upper = preobe_upper/10) %>% 
  pivot_longer(cols = matches("_(lower|upper)$"), names_to = c("cohort", "bound"),
               names_pattern = "(.*)_(lower|upper)$", values_to = "value") %>% 
  pivot_wider(names_from = bound, values_from = value) %>%  
  mutate(mean = (lower + upper) / 2) %>% 
  mutate(lower = lower*SD_logCRP/SD_MDmean) %>% 
  mutate(upper = upper*SD_logCRP/SD_MDmean) %>% 
  mutate(mean = mean*SD_logCRP/SD_MDmean) %>% 
  mutate(cohort = factor(cohort, levels = c("fb1mo", "fb5y", "preobe", "copsych", "genr", "nfbc", "pooled")))
  

#Mean MD
ggplot(MD_results_forest %>% filter(ROI == "Mean") %>% filter(cohort != "pooled")) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", colour = "lightgrey") +
  #geom_vline(xintercept = 6.5, size = 1, linetype = "dashed", colour = "lightgrey") +
  geom_point(aes(x = cohort, y = mean, colour = cohort), size = 2.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = cohort, colour = cohort), size = 1, width = 0.5, alpha = 0.75) +
  scale_colour_manual(values = c("#FF6467", "#7C86FF", "#E95FFC", "#5EE9B5" , "#74D4FF", "#FE9A37", "darkgrey")) +
  
  labs(y = "Standardised effect size") +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(size = 12, angle = 90),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        legend.position = "none")

ggsave("plots/rois/MD_forest_mean.svg", units = "px", width = 750, height = 750)


#MD in 48 ROIs, in 6x8 grid
ggplot(MD_results_forest %>% filter(ROI != "Mean")) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", colour = "lightgrey") +
  geom_vline(xintercept = 6.5, size = 1, linetype = "dashed", colour = "lightgrey") +
  geom_point(aes(x = cohort, y = mean, colour = cohort), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = cohort, colour = cohort), size = 1, width = 0.5, alpha = 0.75) +
  scale_colour_manual(values = c("#FF6467", "#7C86FF", "#E95FFC", "#5EE9B5" , "#74D4FF", "#FE9A37", "darkgrey")) +
  
  facet_wrap(~ROI, ncol = 6, nrow=8) +
  
  labs(y = "Standardised effect size") +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(size = 12, angle = 90),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        legend.position = "none")

ggsave("plots/rois/MD_forest_rois.svg", units = "px", width = 3000, height = 4000)


## FA*sex ##
FA_sex_results_forest <- FA_sex_results %>% 
  rename(pooled_lower = lower_effect, pooled_upper = upper_effect) %>% 
  pivot_longer(cols = matches("_(lower|upper)$"), names_to = c("cohort", "bound"),
               names_pattern = "(.*)_(lower|upper)$", values_to = "value") %>% 
  pivot_wider(names_from = bound, values_from = value) %>%  
  mutate(mean = (lower + upper) / 2) %>% 
  mutate(lower = lower*SD_logCRPsex/SD_FAmean) %>% 
  mutate(upper = upper*SD_logCRPsex/SD_FAmean) %>% 
  mutate(mean = mean*SD_logCRPsex/SD_FAmean) %>% 
  mutate(cohort = factor(cohort, levels = c("fb1mo", "fb5y", "preobe", "copsych", "genr", "nfbc", "pooled")))

#Mean FA_sex
ggplot(FA_sex_results_forest %>% filter(ROI == "Mean")) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", colour = "lightgrey") +
  geom_vline(xintercept = 6.5, size = 1, linetype = "dashed", colour = "lightgrey") +
  geom_point(aes(x = cohort, y = mean, colour = cohort), size = 5) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = cohort, colour = cohort), size = 1, width = 0.5, alpha = 0.75) +
  scale_colour_manual(values = c("#FF6467", "#7C86FF", "#E95FFC", "#5EE9B5" , "#74D4FF", "#FE9A37", "darkgrey")) +
  
  labs(y = "Standardised effect size") +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(size = 12, angle = 90),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        legend.position = "none")

ggsave("plots/rois/FA_sex_forest_mean.svg", units = "px", width = 1000, height = 1000)

## MD*sex ##
MD_sex_results_forest <- MD_sex_results %>% 
  rename(pooled_lower = lower_effect, pooled_upper = upper_effect) %>% 
  pivot_longer(cols = matches("_(lower|upper)$"), names_to = c("cohort", "bound"),
               names_pattern = "(.*)_(lower|upper)$", values_to = "value") %>% 
  pivot_wider(names_from = bound, values_from = value) %>%  
  mutate(mean = (lower + upper) / 2) %>% 
  mutate(lower = lower*SD_logCRPsex/SD_MDmean) %>% 
  mutate(upper = upper*SD_logCRPsex/SD_MDmean) %>% 
  mutate(mean = mean*SD_logCRPsex/SD_MDmean) %>% 
  mutate(cohort = factor(cohort, levels = c("fb1mo", "fb5y", "preobe", "copsych", "genr", "nfbc", "pooled")))

#Mean MD_sex
ggplot(MD_sex_results_forest %>% filter(ROI == "Mean")) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", colour = "lightgrey") +
  geom_vline(xintercept = 6.5, size = 1, linetype = "dashed", colour = "lightgrey") +
  geom_point(aes(x = cohort, y = mean, colour = cohort), size = 5) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = cohort, colour = cohort), size = 1, width = 0.5, alpha = 0.75) +
  scale_colour_manual(values = c("#FF6467", "#7C86FF", "#E95FFC", "#5EE9B5" , "#74D4FF", "#FE9A37", "darkgrey")) +
  
  labs(y = "Standardised effect size") +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(size = 12, angle = 90),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        legend.position = "none")

ggsave("plots/rois/MD_sex_forest_mean.svg", units = "px", width = 1000, height = 1000)


#MD_sex in 48 ROIs, in 6x8 grid
ggplot(MD_sex_results_forest %>% filter(ROI != "Mean")) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", colour = "lightgrey") +
  geom_vline(xintercept = 6.5, size = 1, linetype = "dashed", colour = "lightgrey") +
  geom_point(aes(x = cohort, y = mean, colour = cohort), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = cohort, colour = cohort), size = 1, width = 0.5, alpha = 0.75) +
  scale_colour_manual(values = c("#FF6467", "#7C86FF", "#E95FFC", "#5EE9B5" , "#74D4FF", "#FE9A37", "darkgrey")) +
  
  facet_wrap(~ROI, ncol = 6, nrow=8) +
  
  labs(y = "Standardised effect size") +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(size = 12, angle = 90),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        legend.position = "none")

ggsave("plots/rois/MD_sex_forest_rois.svg", units = "px", width = 3000, height = 4000)



#FA_sex in 48 ROIs, in 6x8 grid
ggplot(FA_sex_results_forest %>% filter(ROI != "Mean")) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", colour = "lightgrey") +
  geom_vline(xintercept = 6.5, size = 1, linetype = "dashed", colour = "lightgrey") +
  geom_point(aes(x = cohort, y = mean, colour = cohort), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = cohort, colour = cohort), size = 1, width = 0.5, alpha = 0.75) +
  scale_colour_manual(values = c("#FF6467", "#7C86FF", "#E95FFC", "#5EE9B5" , "#74D4FF", "#FE9A37", "darkgrey")) +
  
  facet_wrap(~ROI, ncol = 6, nrow=8) +
  
  labs(y = "Standardised effect size") +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(size = 12, angle = 90),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        legend.position = "none")

ggsave("plots/rois/FA_sex_forest_rois.svg", units = "px", width = 3000, height = 4000)


#### Forest plots faceted per ROI - no PREOBE ####

## FA ##
FA_results_no_preobe_forest <- FA_results_no_preobe %>% 
  rename(pooled_lower = lower_effect, pooled_upper = upper_effect) %>% 
  pivot_longer(cols = matches("_(lower|upper)$"), names_to = c("cohort", "bound"),
               names_pattern = "(.*)_(lower|upper)$", values_to = "value") %>% 
  pivot_wider(names_from = bound, values_from = value) %>%  
  mutate(mean = (lower + upper) / 2) %>% 
  mutate(lower = lower*SD_logCRP/SD_FAmean) %>% 
  mutate(upper = upper*SD_logCRP/SD_FAmean) %>% 
  mutate(mean = mean*SD_logCRP/SD_FAmean) %>% 
  mutate(cohort = factor(cohort, levels = c("fb1mo", "fb5y", "copsych", "genr", "nfbc", "pooled")))

#Mean FA
ggplot(FA_results_no_preobe_forest %>% filter(ROI == "Mean")) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", colour = "lightgrey") +
  geom_vline(xintercept = 5.5, size = 1, linetype = "dashed", colour = "lightgrey") +
  geom_point(aes(x = cohort, y = mean, colour = cohort), size = 5) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = cohort, colour = cohort), size = 1, width = 0.5, alpha = 0.75) +
  scale_colour_manual(values = c("#FF6467", "#7C86FF",  "#5EE9B5" , "#74D4FF", "#FE9A37", "darkgrey")) +
  
  labs(y = "Standardised effect size") +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(size = 12, angle = 90),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        legend.position = "none")

ggsave("plots/rois/FA_forest_mean_no_preobe.svg", units = "px", width = 1000, height = 1000)


#FA in 48 ROIs, in 6x8 grid
ggplot(FA_results_no_preobe_forest %>% filter(ROI != "Mean")) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", colour = "lightgrey") +
  geom_vline(xintercept = 5.5, size = 1, linetype = "dashed", colour = "lightgrey") +
  geom_point(aes(x = cohort, y = mean, colour = cohort), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = cohort, colour = cohort), size = 1, width = 0.5, alpha = 0.75) +
  scale_colour_manual(values = c("#FF6467", "#7C86FF",  "#5EE9B5" , "#74D4FF", "#FE9A37", "darkgrey")) +
  
  facet_wrap(~ROI, ncol = 6, nrow=8) +
  
  labs(y = "Standardised effect size") +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(size = 12, angle = 90),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        legend.position = "none")

ggsave("plots/rois/FA_forest_rois_no_preobe.svg", units = "px", width = 3000, height = 4000)

## MD ##
MD_results_no_preobe_forest <- MD_results_no_preobe %>% 
  rename(pooled_lower = lower_effect, pooled_upper = upper_effect) %>% 
  pivot_longer(cols = matches("_(lower|upper)$"), names_to = c("cohort", "bound"),
               names_pattern = "(.*)_(lower|upper)$", values_to = "value") %>% 
  pivot_wider(names_from = bound, values_from = value) %>%  
  mutate(mean = (lower + upper) / 2) %>% 
  mutate(lower = lower*SD_logCRP/SD_MDmean) %>% 
  mutate(upper = upper*SD_logCRP/SD_MDmean) %>% 
  mutate(mean = mean*SD_logCRP/SD_MDmean) %>% 
  mutate(cohort = factor(cohort, levels = c("fb1mo", "fb5y", "copsych", "genr", "nfbc", "pooled")))

#Mean MD
ggplot(MD_results_no_preobe_forest %>% filter(ROI == "Mean")) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", colour = "lightgrey") +
  geom_vline(xintercept = 5.5, size = 1, linetype = "dashed", colour = "lightgrey") +
  geom_point(aes(x = cohort, y = mean, colour = cohort), size = 5) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = cohort, colour = cohort), size = 1, width = 0.5, alpha = 0.75) +
  scale_colour_manual(values = c("#FF6467", "#7C86FF",  "#5EE9B5" , "#74D4FF", "#FE9A37", "darkgrey")) +
  
  labs(y = "Standardised effect size") +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(size = 12, angle = 90),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        legend.position = "none")

ggsave("plots/rois/MD_forest_mean_no_preobe.svg", units = "px", width = 1000, height = 1000)


#MD in 48 ROIs, in 6x8 grid
ggplot(MD_results_no_preobe_forest %>% filter(ROI != "Mean")) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", colour = "lightgrey") +
  geom_vline(xintercept = 5.5, size = 1, linetype = "dashed", colour = "lightgrey") +
  geom_point(aes(x = cohort, y = mean, colour = cohort), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = cohort, colour = cohort), size = 1, width = 0.5, alpha = 0.75) +
  scale_colour_manual(values = c("#FF6467", "#7C86FF",  "#5EE9B5" , "#74D4FF", "#FE9A37", "darkgrey")) +
  
  facet_wrap(~ROI, ncol = 6, nrow=8) +
  
  labs(y = "Standardised effect size") +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(size = 12, angle = 90),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        legend.position = "none")

ggsave("plots/rois/MD_forest_rois_no_preobe.svg", units = "px", width = 3000, height = 4000)


## FA*sex ##
FA_sex_results_no_preobe_forest <- FA_sex_results_no_preobe %>% 
  rename(pooled_lower = lower_effect, pooled_upper = upper_effect) %>% 
  pivot_longer(cols = matches("_(lower|upper)$"), names_to = c("cohort", "bound"),
               names_pattern = "(.*)_(lower|upper)$", values_to = "value") %>% 
  pivot_wider(names_from = bound, values_from = value) %>%  
  mutate(mean = (lower + upper) / 2) %>% 
  mutate(lower = lower*SD_logCRPsex/SD_FAmean) %>% 
  mutate(upper = upper*SD_logCRPsex/SD_FAmean) %>% 
  mutate(mean = mean*SD_logCRPsex/SD_FAmean) %>% 
  mutate(cohort = factor(cohort, levels = c("fb1mo", "fb5y", "copsych", "genr", "nfbc", "pooled")))

#Mean FA_sex
ggplot(FA_sex_results_no_preobe_forest %>% filter(ROI == "Mean")) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", colour = "lightgrey") +
  geom_vline(xintercept = 5.5, size = 1, linetype = "dashed", colour = "lightgrey") +
  geom_point(aes(x = cohort, y = mean, colour = cohort), size = 5) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = cohort, colour = cohort), size = 1, width = 0.5, alpha = 0.75) +
  scale_colour_manual(values = c("#FF6467", "#7C86FF",  "#5EE9B5" , "#74D4FF", "#FE9A37", "darkgrey")) +
  
  labs(y = "Standardised effect size") +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(size = 12, angle = 90),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        legend.position = "none")

ggsave("plots/rois/FA_sex_forest_mean_no_preobe.svg", units = "px", width = 1000, height = 1000)

#FA_sex in 48 ROIs, in 6x8 grid
ggplot(FA_sex_results_no_preobe_forest %>% filter(ROI != "Mean")) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", colour = "lightgrey") +
  geom_vline(xintercept = 5.5, size = 1, linetype = "dashed", colour = "lightgrey") +
  geom_point(aes(x = cohort, y = mean, colour = cohort), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = cohort, colour = cohort), size = 1, width = 0.5, alpha = 0.75) +
  scale_colour_manual(values = c("#FF6467", "#7C86FF",  "#5EE9B5" , "#74D4FF", "#FE9A37", "darkgrey")) +
  
  facet_wrap(~ROI, ncol = 6, nrow=8) +
  
  labs(y = "Standardised effect size") +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(size = 12, angle = 90),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        legend.position = "none")

ggsave("plots/rois/FA_sex_forest_rois_no_preobe.svg", units = "px", width = 3000, height = 4000)

## MD*sex ##
MD_sex_results_no_preobe_forest <- MD_sex_results_no_preobe %>% 
  rename(pooled_lower = lower_effect, pooled_upper = upper_effect) %>% 
  pivot_longer(cols = matches("_(lower|upper)$"), names_to = c("cohort", "bound"),
               names_pattern = "(.*)_(lower|upper)$", values_to = "value") %>% 
  pivot_wider(names_from = bound, values_from = value) %>%  
  mutate(mean = (lower + upper) / 2) %>% 
  mutate(lower = lower*SD_logCRPsex/SD_MDmean) %>% 
  mutate(upper = upper*SD_logCRPsex/SD_MDmean) %>% 
  mutate(mean = mean*SD_logCRPsex/SD_MDmean) %>% 
  mutate(cohort = factor(cohort, levels = c("fb1mo", "fb5y", "copsych", "genr", "nfbc", "pooled")))

#Mean MD_sex
ggplot(MD_sex_results_no_preobe_forest %>% filter(ROI == "Mean")) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", colour = "lightgrey") +
  geom_vline(xintercept = 5.5, size = 1, linetype = "dashed", colour = "lightgrey") +
  geom_point(aes(x = cohort, y = mean, colour = cohort), size = 5) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = cohort, colour = cohort), size = 1, width = 0.5, alpha = 0.75) +
  scale_colour_manual(values = c("#FF6467", "#7C86FF",  "#5EE9B5" , "#74D4FF", "#FE9A37", "darkgrey")) +
  
  labs(y = "Standardised effect size") +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(size = 12, angle = 90),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        legend.position = "none")

ggsave("plots/rois/MD_sex_forest_mean_no_preobe.svg", units = "px", width = 1000, height = 1000)


#MD_sex in 48 ROIs, in 6x8 grid
ggplot(MD_sex_results_no_preobe_forest %>% filter(ROI != "Mean")) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", colour = "lightgrey") +
  geom_vline(xintercept = 5.5, size = 1, linetype = "dashed", colour = "lightgrey") +
  geom_point(aes(x = cohort, y = mean, colour = cohort), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = cohort, colour = cohort), size = 1, width = 0.5, alpha = 0.75) +
  scale_colour_manual(values = c("#FF6467", "#7C86FF",  "#5EE9B5" , "#74D4FF", "#FE9A37", "darkgrey")) +
  
  facet_wrap(~ROI, ncol = 6, nrow=8) +
  
  labs(y = "Standardised effect size") +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(size = 12, angle = 90),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        legend.position = "none")

ggsave("plots/rois/MD_sex_forest_rois_no_preobe.svg", units = "px", width = 3000, height = 4000)







#### Ordered bar plots ####

#FA effect size
ggplot(data = FA_results) + 
  geom_errorbar(aes(ymin = lower_standard, ymax = upper_standard, x = reorder(ROI, effect_standard), colour = effect_standard), 
                size = 1, width = 0.5, alpha = 0.75) +
  geom_bar(aes(x=reorder(ROI, effect_standard), y=effect_standard, fill = effect_standard),
           stat = "identity", colour = "black", alpha = 0.95) + 
  
  scale_fill_gradient2(midpoint = 0, low = "#1F618D", mid = "white", high = "#922B21") +
  scale_colour_gradient2(midpoint = 0, low = "#1F618D", mid = "white", high = "#922B21") +
  
  labs(y = "Standardised effect size") + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, size = 10), 
        legend.position = "none", 
        plot.background = element_blank(), 
        panel.background = element_blank())

ggsave("plots/rois/FA_ordered_barplot.svg", units = "px", width = 3000, height = 1000)

#MD effect size
ggplot(data = MD_results) + 
  geom_errorbar(aes(ymin = lower_standard, ymax = upper_standard, x = reorder(ROI, effect_standard), colour = effect_standard), 
                size = 1, width = 0.5, alpha = 0.75) +
  geom_bar(aes(x=reorder(ROI, effect_standard), y=effect_standard, fill = effect_standard),
           stat = "identity", colour = "black", alpha = 0.95) + 
  
  scale_fill_gradient2(midpoint = 0, low = "#1F618D", mid = "white", high = "#922B21") +
  scale_colour_gradient2(midpoint = 0, low = "#1F618D", mid = "white", high = "#922B21") +
  
  labs(y = "Standardised effect size") + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, size = 10), 
        legend.position = "none", 
        plot.background = element_blank(), 
        panel.background = element_blank())

ggsave("plots/rois/MD_ordered_barplot.svg", units = "px", width = 3000, height = 1000)


#FA sex effect size
ggplot(data = FA_sex_results) + 
  geom_errorbar(aes(ymin = lower_standard, ymax = upper_standard, x = reorder(ROI, effect_standard), colour = effect_standard), 
                size = 1, width = 0.5, alpha = 0.75) +
  geom_bar(aes(x=reorder(ROI, effect_standard), y=effect_standard, fill = effect_standard),
           stat = "identity", colour = "black", alpha = 0.95) + 
  
  scale_fill_gradient2(midpoint = 0, low = "#1F618D", mid = "white", high = "#922B21") +
  scale_colour_gradient2(midpoint = 0, low = "#1F618D", mid = "white", high = "#922B21") +
  
  labs(y = "Standardised effect size") + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, size = 10), 
        legend.position = "none", 
        plot.background = element_blank(), 
        panel.background = element_blank())

ggsave("plots/rois/FA_sex_ordered_barplot.svg", units = "px", width = 3000, height = 1000)

#MD sex effect size
ggplot(data = MD_sex_results) + 
  geom_errorbar(aes(ymin = lower_standard, ymax = upper_standard, x = reorder(ROI, effect_standard), colour = effect_standard), 
                size = 1, width = 0.5, alpha = 0.75) +
  geom_bar(aes(x=reorder(ROI, effect_standard), y=effect_standard, fill = effect_standard),
           stat = "identity", colour = "black", alpha = 0.95) + 
  
  scale_fill_gradient2(midpoint = 0, low = "#1F618D", mid = "white", high = "#922B21") +
  scale_colour_gradient2(midpoint = 0, low = "#1F618D", mid = "white", high = "#922B21") +
  
  labs(y = "Standardised effect size") + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, size = 10), 
        legend.position = "none", 
        plot.background = element_blank(), 
        panel.background = element_blank())

ggsave("plots/rois/MD_sex_ordered_barplot.svg", units = "px", width = 3000, height = 1000)


#### Ordered bar plots - no PREOBE ####

#FA effect size
ggplot(data = FA_results_no_preobe) + 
  geom_errorbar(aes(ymin = lower_standard, ymax = upper_standard, x = reorder(ROI, effect_standard), colour = effect_standard), 
                size = 1, width = 0.5, alpha = 0.75) +
  geom_bar(aes(x=reorder(ROI, effect_standard), y=effect_standard, fill = effect_standard),
           stat = "identity", colour = "black", alpha = 0.95) + 
  
  scale_fill_gradient2(midpoint = 0, low = "#1F618D", mid = "white", high = "#922B21") +
  scale_colour_gradient2(midpoint = 0, low = "#1F618D", mid = "white", high = "#922B21") +
  
  labs(y = "Standardised effect size") + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, size = 10), 
        legend.position = "none", 
        plot.background = element_blank(), 
        panel.background = element_blank())

ggsave("plots/rois/FA_ordered_barplot_no_preobe.svg", units = "px", width = 3000, height = 1000)

#MD effect size
ggplot(data = MD_results_no_preobe) + 
  geom_errorbar(aes(ymin = lower_standard, ymax = upper_standard, x = reorder(ROI, effect_standard), colour = effect_standard), 
                size = 1, width = 0.5, alpha = 0.75) +
  geom_bar(aes(x=reorder(ROI, effect_standard), y=effect_standard, fill = effect_standard),
           stat = "identity", colour = "black", alpha = 0.95) + 
  
  scale_fill_gradient2(midpoint = 0, low = "#1F618D", mid = "white", high = "#922B21") +
  scale_colour_gradient2(midpoint = 0, low = "#1F618D", mid = "white", high = "#922B21") +
  
  labs(y = "Standardised effect size") + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, size = 10), 
        legend.position = "none", 
        plot.background = element_blank(), 
        panel.background = element_blank())

ggsave("plots/rois/MD_ordered_barplot_no_preobe.svg", units = "px", width = 3000, height = 1000)


#FA sex effect size
ggplot(data = FA_sex_results_no_preobe) + 
  geom_errorbar(aes(ymin = lower_standard, ymax = upper_standard, x = reorder(ROI, effect_standard), colour = effect_standard), 
                size = 1, width = 0.5, alpha = 0.75) +
  geom_bar(aes(x=reorder(ROI, effect_standard), y=effect_standard, fill = effect_standard),
           stat = "identity", colour = "black", alpha = 0.95) + 
  
  scale_fill_gradient2(midpoint = 0, low = "#1F618D", mid = "white", high = "#922B21") +
  scale_colour_gradient2(midpoint = 0, low = "#1F618D", mid = "white", high = "#922B21") +
  
  labs(y = "Standardised effect size") + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, size = 10), 
        legend.position = "none", 
        plot.background = element_blank(), 
        panel.background = element_blank())

ggsave("plots/rois/FA_sex_ordered_barplot_no_preobe.svg", units = "px", width = 3000, height = 1000)

#MD sex effect size
ggplot(data = MD_sex_results_no_preobe) + 
  geom_errorbar(aes(ymin = lower_standard, ymax = upper_standard, x = reorder(ROI, effect_standard), colour = effect_standard), 
                size = 1, width = 0.5, alpha = 0.75) +
  geom_bar(aes(x=reorder(ROI, effect_standard), y=effect_standard, fill = effect_standard),
           stat = "identity", colour = "black", alpha = 0.95) + 
  
  scale_fill_gradient2(midpoint = 0, low = "#1F618D", mid = "white", high = "#922B21") +
  scale_colour_gradient2(midpoint = 0, low = "#1F618D", mid = "white", high = "#922B21") +
  
  labs(y = "Standardised effect size") + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, size = 10), 
        legend.position = "none", 
        plot.background = element_blank(), 
        panel.background = element_blank())

ggsave("plots/rois/MD_sex_ordered_barplot_no_preobe.svg", units = "px", width = 3000, height = 1000)

