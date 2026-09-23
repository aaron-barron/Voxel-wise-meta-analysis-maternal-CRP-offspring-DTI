#### Libraries ####

#Define R library path for external libraries (not needed on most computers, but necessary on CSC Puhti)
.libPaths(c("/projappl/project_2006897/project_rpackages_4.4.0", .libPaths()))
libpath <- .libPaths()[1]

#Load R packages used in this script
library("tidyverse") 
library("data.table")
library("janitor")
library("svglite")
library("meta")

sites <- c("fb1mo", "fb5y", "preobe", "copsych", "genr", "nfbc")

#### Age for meta-regression ####

ages <- c(0.07, 5.4, 10, 9.9, 26.5)

#### FA ####

#Effect sizes
FA_cope <- read.csv("data/ROIs/FA_cope_ROIs.csv", header = F)[-3, ] %>% mutate(Mean = rowMeans(.))

colnames(FA_cope) <- c("MCP", "PCT", "GCC", "BCC", "SCC", "FX",
                       "CST_R", "CST_L", "ML_R", "ML_L", "ICP_R", "ICP_L",
                       "SCP_R", "SCP_L", "CP_R", "CP_L", "ALIC_R", "ALIC_L",
                       "PLIC_R", "PLIC_L", "RLIC_R", "RLIC_L", "ACR_R", "ACR_L",
                       "SCR_R", "SCR_L", "PCR_R", "PCR_L", "PTR_R", "PTR_L",
                       "SS_R", "SS_L", "EC_R", "EC_L", "CGC_R", "CGC_L",
                       "CGH_R", "CGH_L", "FX_ST_R", "FX_ST_L", "SLF_R", "SLF_L",
                       "SFOF_R", "SFOF_L", "IFOF_R", "IFOF_L", "UNC_R", "UNC_L", "Mean")

#Errors
FA_varcope <- read.csv("data/ROIs/FA_varcope_ROIs.csv", header = F)[-3, ] %>% mutate(Mean = rowMeans(.))

colnames(FA_varcope) <- c("MCP", "PCT", "GCC", "BCC", "SCC", "FX",
                          "CST_R", "CST_L", "ML_R", "ML_L", "ICP_R", "ICP_L",
                          "SCP_R", "SCP_L", "CP_R", "CP_L", "ALIC_R", "ALIC_L",
                          "PLIC_R", "PLIC_L", "RLIC_R", "RLIC_L", "ACR_R", "ACR_L",
                          "SCR_R", "SCR_L", "PCR_R", "PCR_L", "PTR_R", "PTR_L",
                          "SS_R", "SS_L", "EC_R", "EC_L", "CGC_R", "CGC_L",
                          "CGH_R", "CGH_L", "FX_ST_R", "FX_ST_L", "SLF_R", "SLF_L",
                          "SFOF_R", "SFOF_L", "IFOF_R", "IFOF_L", "UNC_R", "UNC_L", "Mean")

FA_SE <- sqrt(FA_varcope)
rm(FA_varcope)


#Meta-analysis
#Loop meta-analysis model over voxels in parallel
FA_results <- data.frame()

for(i in colnames(FA_cope)){
  print(i)
  
  #Define random-effects meta-analysis for voxel i
  model <- metagen(TE = FA_cope[[i]], seTE = FA_SE[[i]], common = FALSE, random = TRUE,
                   control = list(maxiter = 10000, stepadj = 0.5))
  
  #Meta regression with age
  regression <- metareg(model, ~ ages)
  
  #Extract voxel ID
  ROI <- i
  
  #Extract model results
  effect <- model$TE.random
  lower_effect <- model$lower.random
  upper_effect  <- model$upper.random
  error <- model$seTE.random
  tstat <- model$statistic.random
  pval <- model$pval.random
  
  #extract weights from each site
  weights <- model$w.random
  weights <- if(length(weights) < 6) c(weights, rep(NA, 6 - length(weights))) else weights[1:6]
  
  fb1mo <- weights[1]
  fb5y <- weights[2]
  copsych <- weights[3]
  genr <- weights[4]
  nfbc <- weights[5]
  
  #lower estimates
  lower <- model$lower
  lower <- if(length(lower) < 6) c(lower, rep(NA, 6 - length(lower))) else lower[1:6]
  
  fb1mo_lower <- lower[1]
  fb5y_lower <- lower[2]
  copsych_lower <- lower[3]
  genr_lower <- lower[4]
  nfbc_lower <- lower[5]
  
  #upper estimates
  upper <- model$upper
  upper <- if(length(upper) < 6) c(upper, rep(NA, 6 - length(upper))) else upper[1:6]
  
  fb1mo_upper <- upper[1]
  fb5y_upper <- upper[2]
  copsych_upper <- upper[3]
  genr_upper <- upper[4]
  nfbc_upper <- upper[5]
  
  #extract results from meta-regression against age
  age_coef<- regression$beta[2]
  age_pval  <- regression$pval[2]
  
  #Combine results
  i_results <- data.frame(ROI, effect, error, tstat, pval, fb1mo, fb5y, copsych, genr, nfbc, lower_effect, upper_effect,
                          fb1mo_lower, fb5y_lower, copsych_lower, genr_lower, nfbc_lower,
                          fb1mo_upper, fb5y_upper, copsych_upper, genr_upper, nfbc_upper, age_coef, age_pval)
  FA_results <- rbind(FA_results, i_results)
  
}

#Write results
write_rds(FA_results, "results/FA_ROI_no_preobe_results.rds")


#### FA_sex ####

#Effect sizes
FA_sex_cope <- read.csv("data/ROIs/FA_sex_cope_ROIs.csv", header = F)[-3, ] %>% mutate(Mean = rowMeans(.))

colnames(FA_sex_cope) <- c("MCP", "PCT", "GCC", "BCC", "SCC", "FX",
                           "CST_R", "CST_L", "ML_R", "ML_L", "ICP_R", "ICP_L",
                           "SCP_R", "SCP_L", "CP_R", "CP_L", "ALIC_R", "ALIC_L",
                           "PLIC_R", "PLIC_L", "RLIC_R", "RLIC_L", "ACR_R", "ACR_L",
                           "SCR_R", "SCR_L", "PCR_R", "PCR_L", "PTR_R", "PTR_L",
                           "SS_R", "SS_L", "EC_R", "EC_L", "CGC_R", "CGC_L",
                           "CGH_R", "CGH_L", "FX_ST_R", "FX_ST_L", "SLF_R", "SLF_L",
                           "SFOF_R", "SFOF_L", "IFOF_R", "IFOF_L", "UNC_R", "UNC_L", "Mean")

#Errors
FA_sex_varcope <- read.csv("data/ROIs/FA_sex_varcope_ROIs.csv", header = F)[-3, ] %>% mutate(Mean = rowMeans(.))

colnames(FA_sex_varcope) <- c("MCP", "PCT", "GCC", "BCC", "SCC", "FX",
                              "CST_R", "CST_L", "ML_R", "ML_L", "ICP_R", "ICP_L",
                              "SCP_R", "SCP_L", "CP_R", "CP_L", "ALIC_R", "ALIC_L",
                              "PLIC_R", "PLIC_L", "RLIC_R", "RLIC_L", "ACR_R", "ACR_L",
                              "SCR_R", "SCR_L", "PCR_R", "PCR_L", "PTR_R", "PTR_L",
                              "SS_R", "SS_L", "EC_R", "EC_L", "CGC_R", "CGC_L",
                              "CGH_R", "CGH_L", "FX_ST_R", "FX_ST_L", "SLF_R", "SLF_L",
                              "SFOF_R", "SFOF_L", "IFOF_R", "IFOF_L", "UNC_R", "UNC_L", "Mean")

FA_sex_SE <- sqrt(FA_sex_varcope)
rm(FA_sex_varcope)


#Meta-analysis
#Loop meta-analysis model over voxels in parallel
FA_sex_results <- data.frame()

for(i in colnames(FA_sex_cope)){
  print(i)
  
  #Define random-effects meta-analysis for voxel i
  model <- metagen(TE = FA_sex_cope[[i]], seTE = FA_sex_SE[[i]], common = FALSE, random = TRUE,
                   control = list(maxiter = 10000, stepadj = 0.5))
  
  #Meta regression with age
  regression <- metareg(model, ~ ages)
  
  #Extract voxel ID
  ROI <- i
  
  #Extract model results
  effect <- model$TE.random
  lower_effect <- model$lower.random
  upper_effect  <- model$upper.random
  error <- model$seTE.random
  tstat <- model$statistic.random
  pval <- model$pval.random
  
  #extract weights from each site
  weights <- model$w.random
  weights <- if(length(weights) < 6) c(weights, rep(NA, 6 - length(weights))) else weights[1:6]
  
  fb1mo <- weights[1]
  fb5y <- weights[2]
  copsych <- weights[3]
  genr <- weights[4]
  nfbc <- weights[5]
  
  #lower estimates
  lower <- model$lower
  lower <- if(length(lower) < 6) c(lower, rep(NA, 6 - length(lower))) else lower[1:6]
  
  fb1mo_lower <- lower[1]
  fb5y_lower <- lower[2]
  copsych_lower <- lower[3]
  genr_lower <- lower[4]
  nfbc_lower <- lower[5]
  
  #upper estimates
  upper <- model$upper
  upper <- if(length(upper) < 6) c(upper, rep(NA, 6 - length(upper))) else upper[1:6]
  
  fb1mo_upper <- upper[1]
  fb5y_upper <- upper[2]
  copsych_upper <- upper[3]
  genr_upper <- upper[4]
  nfbc_upper <- upper[5]
  
  #extract results from meta-regression against age
  age_coef<- regression$beta[2]
  age_pval  <- regression$pval[2]
  
  #Combine results
  i_results <- data.frame(ROI, effect, error, tstat, pval, fb1mo, fb5y, copsych, genr, nfbc, lower_effect, upper_effect,
                          fb1mo_lower, fb5y_lower, copsych_lower, genr_lower, nfbc_lower,
                          fb1mo_upper, fb5y_upper, copsych_upper, genr_upper, nfbc_upper, age_coef, age_pval)
  FA_sex_results <- rbind(FA_sex_results, i_results)
  
}

#Write results
write_rds(FA_sex_results, "results/FA_sex_ROI_no_preobe_results.rds")


#### MD ####

#Effect sizes
MD_cope <- read.csv("data/ROIs/MD_cope_ROIs.csv", header = F)[-3, ] %>% mutate(Mean = rowMeans(.))

colnames(MD_cope) <- c("MCP", "PCT", "GCC", "BCC", "SCC", "FX",
                       "CST_R", "CST_L", "ML_R", "ML_L", "ICP_R", "ICP_L",
                       "SCP_R", "SCP_L", "CP_R", "CP_L", "ALIC_R", "ALIC_L",
                       "PLIC_R", "PLIC_L", "RLIC_R", "RLIC_L", "ACR_R", "ACR_L",
                       "SCR_R", "SCR_L", "PCR_R", "PCR_L", "PTR_R", "PTR_L",
                       "SS_R", "SS_L", "EC_R", "EC_L", "CGC_R", "CGC_L",
                       "CGH_R", "CGH_L", "FX_ST_R", "FX_ST_L", "SLF_R", "SLF_L",
                       "SFOF_R", "SFOF_L", "IFOF_R", "IFOF_L", "UNC_R", "UNC_L", "Mean")

#Errors
MD_varcope <- read.csv("data/ROIs/MD_varcope_ROIs.csv", header = F)[-3, ] %>% mutate(Mean = rowMeans(.))

colnames(MD_varcope) <- c("MCP", "PCT", "GCC", "BCC", "SCC", "FX",
                          "CST_R", "CST_L", "ML_R", "ML_L", "ICP_R", "ICP_L",
                          "SCP_R", "SCP_L", "CP_R", "CP_L", "ALIC_R", "ALIC_L",
                          "PLIC_R", "PLIC_L", "RLIC_R", "RLIC_L", "ACR_R", "ACR_L",
                          "SCR_R", "SCR_L", "PCR_R", "PCR_L", "PTR_R", "PTR_L",
                          "SS_R", "SS_L", "EC_R", "EC_L", "CGC_R", "CGC_L",
                          "CGH_R", "CGH_L", "FX_ST_R", "FX_ST_L", "SLF_R", "SLF_L",
                          "SFOF_R", "SFOF_L", "IFOF_R", "IFOF_L", "UNC_R", "UNC_L", "Mean")

MD_SE <- sqrt(MD_varcope)
rm(MD_varcope)


#Meta-analysis
#Loop meta-analysis model over voxels in parallel
MD_results <- data.frame()

for(i in colnames(MD_cope)){
  print(i)
  
  #Define random-effects meta-analysis for voxel i
  model <- metagen(TE = MD_cope[[i]], seTE = MD_SE[[i]], common = FALSE, random = TRUE,
                   control = list(maxiter = 10000, stepadj = 0.5))
  
  #Meta regression with age
  regression <- metareg(model, ~ ages)
  
  #Extract voxel ID
  ROI <- i
  
  #Extract model results
  effect <- model$TE.random
  lower_effect <- model$lower.random
  upper_effect  <- model$upper.random
  error <- model$seTE.random
  tstat <- model$statistic.random
  pval <- model$pval.random
  
  #extract weights from each site
  weights <- model$w.random
  weights <- if(length(weights) < 6) c(weights, rep(NA, 6 - length(weights))) else weights[1:6]
  
  fb1mo <- weights[1]
  fb5y <- weights[2]
  copsych <- weights[3]
  genr <- weights[4]
  nfbc <- weights[5]
  
  #lower estimates
  lower <- model$lower
  lower <- if(length(lower) < 6) c(lower, rep(NA, 6 - length(lower))) else lower[1:6]
  
  fb1mo_lower <- lower[1]
  fb5y_lower <- lower[2]
  copsych_lower <- lower[3]
  genr_lower <- lower[4]
  nfbc_lower <- lower[5]
  
  #upper estimates
  upper <- model$upper
  upper <- if(length(upper) < 6) c(upper, rep(NA, 6 - length(upper))) else upper[1:6]
  
  fb1mo_upper <- upper[1]
  fb5y_upper <- upper[2]
  copsych_upper <- upper[3]
  genr_upper <- upper[4]
  nfbc_upper <- upper[5]
  
  #extract results from meta-regression against age
  age_coef<- regression$beta[2]
  age_pval  <- regression$pval[2]
  
  #Combine results
  i_results <- data.frame(ROI, effect, error, tstat, pval, fb1mo, fb5y, copsych, genr, nfbc, lower_effect, upper_effect,
                          fb1mo_lower, fb5y_lower, copsych_lower, genr_lower, nfbc_lower,
                          fb1mo_upper, fb5y_upper, copsych_upper, genr_upper, nfbc_upper, age_coef, age_pval)
  MD_results <- rbind(MD_results, i_results)
  
}

#Write results
write_rds(MD_results, "results/MD_ROI_no_preobe_results.rds")


#### MD_sex ####

#Effect sizes
MD_sex_cope <- read.csv("data/ROIs/MD_sex_cope_ROIs.csv", header = F)[-3, ] %>% mutate(Mean = rowMeans(.))

colnames(MD_sex_cope) <- c("MCP", "PCT", "GCC", "BCC", "SCC", "FX",
                           "CST_R", "CST_L", "ML_R", "ML_L", "ICP_R", "ICP_L",
                           "SCP_R", "SCP_L", "CP_R", "CP_L", "ALIC_R", "ALIC_L",
                           "PLIC_R", "PLIC_L", "RLIC_R", "RLIC_L", "ACR_R", "ACR_L",
                           "SCR_R", "SCR_L", "PCR_R", "PCR_L", "PTR_R", "PTR_L",
                           "SS_R", "SS_L", "EC_R", "EC_L", "CGC_R", "CGC_L",
                           "CGH_R", "CGH_L", "FX_ST_R", "FX_ST_L", "SLF_R", "SLF_L",
                           "SFOF_R", "SFOF_L", "IFOF_R", "IFOF_L", "UNC_R", "UNC_L", "Mean")

#Errors
MD_sex_varcope <- read.csv("data/ROIs/MD_sex_varcope_ROIs.csv", header = F)[-3, ] %>% mutate(Mean = rowMeans(.))

colnames(MD_sex_varcope) <- c("MCP", "PCT", "GCC", "BCC", "SCC", "FX",
                              "CST_R", "CST_L", "ML_R", "ML_L", "ICP_R", "ICP_L",
                              "SCP_R", "SCP_L", "CP_R", "CP_L", "ALIC_R", "ALIC_L",
                              "PLIC_R", "PLIC_L", "RLIC_R", "RLIC_L", "ACR_R", "ACR_L",
                              "SCR_R", "SCR_L", "PCR_R", "PCR_L", "PTR_R", "PTR_L",
                              "SS_R", "SS_L", "EC_R", "EC_L", "CGC_R", "CGC_L",
                              "CGH_R", "CGH_L", "FX_ST_R", "FX_ST_L", "SLF_R", "SLF_L",
                              "SFOF_R", "SFOF_L", "IFOF_R", "IFOF_L", "UNC_R", "UNC_L", "Mean")

MD_sex_SE <- sqrt(MD_sex_varcope)
rm(MD_sex_varcope)


#Meta-analysis
#Loop meta-analysis model over voxels in parallel
MD_sex_results <- data.frame()

for(i in colnames(MD_sex_cope)){
  print(i)
  
  #Define random-effects meta-analysis for voxel i
  model <- metagen(TE = MD_sex_cope[[i]], seTE = MD_sex_SE[[i]], common = FALSE, random = TRUE,
                   control = list(maxiter = 10000, stepadj = 0.5))
  
  #Meta regression with age
  regression <- metareg(model, ~ ages)
  
  #Extract voxel ID
  ROI <- i
  
  #Extract model results
  effect <- model$TE.random
  lower_effect <- model$lower.random
  upper_effect  <- model$upper.random
  error <- model$seTE.random
  tstat <- model$statistic.random
  pval <- model$pval.random
  
  #extract weights from each site
  weights <- model$w.random
  weights <- if(length(weights) < 6) c(weights, rep(NA, 6 - length(weights))) else weights[1:6]
  
  fb1mo <- weights[1]
  fb5y <- weights[2]
  copsych <- weights[3]
  genr <- weights[4]
  nfbc <- weights[5]
  
  #lower estimates
  lower <- model$lower
  lower <- if(length(lower) < 6) c(lower, rep(NA, 6 - length(lower))) else lower[1:6]
  
  fb1mo_lower <- lower[1]
  fb5y_lower <- lower[2]
  copsych_lower <- lower[3]
  genr_lower <- lower[4]
  nfbc_lower <- lower[5]
  
  #upper estimates
  upper <- model$upper
  upper <- if(length(upper) < 6) c(upper, rep(NA, 6 - length(upper))) else upper[1:6]
  
  fb1mo_upper <- upper[1]
  fb5y_upper <- upper[2]
  copsych_upper <- upper[3]
  genr_upper <- upper[4]
  nfbc_upper <- upper[5]
  
  #extract results from meta-regression against age
  age_coef<- regression$beta[2]
  age_pval  <- regression$pval[2]
  
  #Combine results
  i_results <- data.frame(ROI, effect, error, tstat, pval, fb1mo, fb5y, copsych, genr, nfbc, lower_effect, upper_effect,
                          fb1mo_lower, fb5y_lower, copsych_lower, genr_lower, nfbc_lower,
                          fb1mo_upper, fb5y_upper, copsych_upper, genr_upper, nfbc_upper, age_coef, age_pval)
  MD_sex_results <- rbind(MD_sex_results, i_results)
  
}

#Write results
write_rds(MD_sex_results, "results/MD_sex_ROI_no_preobe_results.rds")
