#### Load libraries ####

#Define R library path for external libraries (not needed on most computers, but necessary on CSC Puhti)
.libPaths(c("/projappl/project_2006897/project_rpackages_4.4.0", .libPaths()))
libpath <- .libPaths()[1]

#Load R packages used in this script
library("tidyverse") 
library("data.table")
library("janitor")
library("meta")
library("doMPI")

#### Age for meta-regression ####

ages <- c(0.07, 5.4, 6.4, 10, 9.9, 26.5)

#### Set up parallel environment ####
cluster <- startMPIcluster(count=39)
registerDoMPI(cluster)
cat("Number of workers: ", getDoParWorkers(), "\n")

print(paste("Script started at", Sys.time()))

exportDoMPI(cluster, "ages")

#### FA Meta analysis for positive association ####

#Load data
FA_cope_positive <- fread("data/intercept/FA_cope_positive_FA_Skiftidata.txt", skip = 8) %>% as.data.frame()
FA_varcope_positive <- fread("data/intercept/FA_varcope_positive_FA_Skiftidata.txt", skip = 8) %>% as.data.frame()

#Extract site name
site <- FA_cope_positive[, 1]
FA_cope_positive <- FA_cope_positive[, -1]
FA_varcope_positive <- FA_varcope_positive[, -1]

#Estimate standard errors of effect size
FA_SE_positive <- sqrt(FA_varcope_positive)
rm(FA_varcope_positive)

#Split data into chunks for efficient parallelisation
chunks <- split(1:ncol(FA_cope_positive), cut(1:ncol(FA_cope_positive), 39))

#Time the loop
start <- Sys.time()

#Loop meta-analysis model over voxels in parallel
meta_results <- foreach(k = chunks, .combine = "rbind", .export = c("ages")) %dopar% {
  
  k_results <- data.frame()
  
  #Define a for-loop within each of the workers
  for(i in seq_along(k)){
  
  #Column index
  column_index <- k[i]
    
  #Define random-effects meta-analysis for voxel i
  model <- metagen(TE = FA_cope_positive[, column_index], seTE = FA_SE_positive[, column_index], common = FALSE, random = TRUE,
                   ccontrol = list(maxiter = 10000, stepadj = 0.5), method.tau = "DL")
  
  #Meta-regression against age
  regression <- metareg(model, ~ ages)
  
  #Extract voxel ID
  voxel <- colnames(FA_cope_positive)[column_index]
  
  #Extract model results
  effect <- model$TE.random
  lower <- model$lower.random
  upper <- model$upper.random
  error <- model$seTE.random
  tstat <- model$statistic.random
  pval <- model$pval.random
  
  #extract weights from each site
  weights <- model$w.random
  weights <- if(length(weights) < 6) c(weights, rep(NA, 6 - length(weights))) else weights[1:6]
  
  fb1mo <- weights[1]
  fb5y <- weights[2]
  preobe <- weights[3]
  copsych <- weights[4]
  genr <- weights[5]
  nfbc <- weights[6]
  
  #extract results from meta-regression against age
  age_coef<- regression$beta[2]
  age_pval  <- regression$pval[2]
  
  #Combine results
  i_results <- data.frame(voxel, effect, error, tstat, pval, fb1mo, fb5y, preobe, copsych, genr, nfbc, age_coef, age_pval)
  k_results <- rbind(k_results, i_results)
  
  }
  
  k_results
  
}

end <- Sys.time()
print(paste(c("Time taken =", end - start)))

#Write results
write_rds(meta_results, "results/FA_positive_metagen_results.rds")

#Clear environment
rm(FA_cope_positive)
rm(FA_SE_positive)
rm(meta_results)


#### FA sex Meta analysis for positive association ####

#Load data
FA_cope_positive <- fread("data/intercept/FA_sex_cope_positive_FA_Skiftidata.txt", skip = 8) %>% as.data.frame()
FA_varcope_positive <- fread("data/intercept/FA_sex_varcope_positive_FA_Skiftidata.txt", skip = 8) %>% as.data.frame()

#Extract site name
site <- FA_cope_positive[, 1]
FA_cope_positive <- FA_cope_positive[, -1]
FA_varcope_positive <- FA_varcope_positive[, -1]

#Estimate standard errors of effect size
FA_SE_positive <- sqrt(FA_varcope_positive)
rm(FA_varcope_positive)

#Split data into chunks for efficient parallelisation
chunks <- split(1:ncol(FA_cope_positive), cut(1:ncol(FA_cope_positive), 39))

#Time the loop
start <- Sys.time()

#Loop meta-analysis model over voxels in parallel
meta_results <- foreach(k = chunks, .combine = "rbind", .export = c("ages")) %dopar% {
  
  k_results <- data.frame()
  
  #Define a for-loop within each of the workers
  for(i in seq_along(k)){
    
    #Column index
    column_index <- k[i]
    
    #Define random-effects meta-analysis for voxel i
    model <- metagen(TE = FA_cope_positive[, column_index], seTE = FA_SE_positive[, column_index], common = FALSE, random = TRUE,
                     ccontrol = list(maxiter = 10000, stepadj = 0.5), method.tau = "DL")
    
    #Meta-regression against age
    regression <- metareg(model, ~ ages)
    
    #Extract voxel ID
    voxel <- colnames(FA_cope_positive)[column_index]
    
    #Extract model results
    effect <- model$TE.random
    lower <- model$lower.random
    upper <- model$upper.random
    error <- model$seTE.random
    tstat <- model$statistic.random
    pval <- model$pval.random
    
    #extract weights from each site
    weights <- model$w.random
    weights <- if(length(weights) < 6) c(weights, rep(NA, 6 - length(weights))) else weights[1:6]
    
    fb1mo <- weights[1]
    fb5y <- weights[2]
    preobe <- weights[3]
    copsych <- weights[4]
    genr <- weights[5]
    nfbc <- weights[6]
    
    #extract results from meta-regression against age
    age_coef<- regression$beta[2]
    age_pval  <- regression$pval[2]
    
    #Combine results
    i_results <- data.frame(voxel, effect, error, tstat, pval, fb1mo, fb5y, preobe, copsych, genr, nfbc, age_coef, age_pval)
    k_results <- rbind(k_results, i_results)
    
  }
  
  k_results
  
}

end <- Sys.time()
print(paste(c("Time taken =", end - start)))

#Write results
write_rds(meta_results, "results/FA_sex_positive_metagen_results.rds")

#Clear environment
rm(FA_cope_positive)
rm(FA_SE_positive)
rm(meta_results)


#### MD Meta analysis for positive association ####

#Load data
MD_cope_positive <- fread("data/intercept/MD_cope_positive_MD_Skiftidata.txt", skip = 8) %>% as.data.frame()
MD_varcope_positive <- fread("data/intercept/MD_varcope_positive_MD_Skiftidata.txt", skip = 8) %>% as.data.frame()

#Extract site name
site <- MD_cope_positive[, 1]
MD_cope_positive <- MD_cope_positive[, -1]
MD_varcope_positive <- MD_varcope_positive[, -1]

#Estimate standard errors of effect size
MD_SE_positive <- sqrt(MD_varcope_positive)
rm(MD_varcope_positive)

#Split data into chunks for efficient parallelisation
chunks <- split(1:ncol(MD_cope_positive), cut(1:ncol(MD_cope_positive), 39))

#Time the loop
start <- Sys.time()

#Loop meta-analysis model over voxels in parallel
meta_results <- foreach(k = chunks, .combine = "rbind", .export = c("ages")) %dopar% {
  
  k_results <- data.frame()
  
  #Define a for-loop within each of the workers
  for(i in seq_along(k)){
    
    #Column index
    column_index <- k[i]
    
    #Define random-effects meta-analysis for voxel i
    model <- metagen(TE = MD_cope_positive[, column_index], seTE = MD_SE_positive[, column_index], common = FALSE, random = TRUE,
                     ccontrol = list(maxiter = 10000, stepadj = 0.5), method.tau = "DL")
    
    #Meta-regression against age
    regression <- metareg(model, ~ ages)
    
    #Extract voxel ID
    voxel <- colnames(MD_cope_positive)[column_index]
    
    #Extract model results
    effect <- model$TE.random
    lower <- model$lower.random
    upper <- model$upper.random
    error <- model$seTE.random
    tstat <- model$statistic.random
    pval <- model$pval.random
    
    #extract weights from each site
    weights <- model$w.random
    weights <- if(length(weights) < 6) c(weights, rep(NA, 6 - length(weights))) else weights[1:6]
    
    fb1mo <- weights[1]
    fb5y <- weights[2]
    preobe <- weights[3]
    copsych <- weights[4]
    genr <- weights[5]
    nfbc <- weights[6]
    
    #extract results from meta-regression against age
    age_coef<- regression$beta[2]
    age_pval  <- regression$pval[2]
    
    #Combine results
    i_results <- data.frame(voxel, effect, error, tstat, pval, fb1mo, fb5y, preobe, copsych, genr, nfbc, age_coef, age_pval)
    k_results <- rbind(k_results, i_results)
    
  }
  
  k_results
  
}

end <- Sys.time()
print(paste(c("Time taken =", end - start)))

#Write results
write_rds(meta_results, "results/MD_positive_metagen_results.rds")

#Clear environment
rm(MD_cope_positive)
rm(MD_SE_positive)
rm(meta_results)


#### MD sex Meta analysis for positive association ####

#Load data
MD_cope_positive <- fread("data/intercept/MD_sex_cope_positive_MD_Skiftidata.txt", skip = 8) %>% as.data.frame()
MD_varcope_positive <- fread("data/intercept/MD_sex_varcope_positive_MD_Skiftidata.txt", skip = 8) %>% as.data.frame()

#Extract site name
site <- MD_cope_positive[, 1]
MD_cope_positive <- MD_cope_positive[, -1]
MD_varcope_positive <- MD_varcope_positive[, -1]

#Estimate standard errors of effect size
MD_SE_positive <- sqrt(MD_varcope_positive)
rm(MD_varcope_positive)

#Split data into chunks for efficient parallelisation
chunks <- split(1:ncol(MD_cope_positive), cut(1:ncol(MD_cope_positive), 39))

#Time the loop
start <- Sys.time()

#Loop meta-analysis model over voxels in parallel
meta_results <- foreach(k = chunks, .combine = "rbind", .export = c("ages")) %dopar% {
  
  k_results <- data.frame()
  
  #Define a for-loop within each of the workers
  for(i in seq_along(k)){
    
    #Column index
    column_index <- k[i]
    
    #Define random-effects meta-analysis for voxel i
    model <- metagen(TE = MD_cope_positive[, column_index], seTE = MD_SE_positive[, column_index], common = FALSE, random = TRUE,
                     ccontrol = list(maxiter = 10000, stepadj = 0.5), method.tau = "DL")
    
    #Meta-regression against age
    regression <- metareg(model, ~ ages)
    
    #Extract voxel ID
    voxel <- colnames(MD_cope_positive)[column_index]
    
    #Extract model results
    effect <- model$TE.random
    lower <- model$lower.random
    upper <- model$upper.random
    error <- model$seTE.random
    tstat <- model$statistic.random
    pval <- model$pval.random
    
    #extract weights from each site
    weights <- model$w.random
    weights <- if(length(weights) < 6) c(weights, rep(NA, 6 - length(weights))) else weights[1:6]
    
    fb1mo <- weights[1]
    fb5y <- weights[2]
    preobe <- weights[3]
    copsych <- weights[4]
    genr <- weights[5]
    nfbc <- weights[6]
    
    #extract results from meta-regression against age
    age_coef<- regression$beta[2]
    age_pval  <- regression$pval[2]
    
    #Combine results
    i_results <- data.frame(voxel, effect, error, tstat, pval, fb1mo, fb5y, preobe, copsych, genr, nfbc, age_coef, age_pval)
    k_results <- rbind(k_results, i_results)
    
  }
  
  k_results
  
}

end <- Sys.time()
print(paste(c("Time taken =", end - start)))

#Write results
write_rds(meta_results, "results/MD_sex_positive_metagen_results.rds")

#Clear environment
rm(MD_cope_positive)
rm(MD_SE_positive)
rm(meta_results)


#### Close cluster ####
closeCluster(cluster)
mpi.quit()