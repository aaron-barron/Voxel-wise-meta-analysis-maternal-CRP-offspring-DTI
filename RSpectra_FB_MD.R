#### Load libraries ####

#Define R library path for external libraries (not needed on most computers, but necessary on CSC Puhti)
.libPaths(c("/projappl/project_2006897/project_rpackages_4.4.0", .libPaths()))
libpath <- .libPaths()[1]

#Load R packages used in this script
library("tidyverse") 
library("data.table")
library("poolr")
library("RSpectra")


#Collect voxel names
all_voxels <- readRDS("results/MD_positive_metagen_results.rds")[, 1]
common_voxels <- (readRDS("results/MD_positive_metagen_results.rds") %>% filter(genr != 0))[, 1]

#### MD ####


#Load correlation matrix
MD_cormat <- read_rds("data/FB_5Y_MD_cormat_result.rds")

#Filter only to the 106,220 voxels that are common to all sites
MD_cormat <- as.matrix(MD_cormat)
rownames(MD_cormat) <- all_voxels
colnames(MD_cormat) <- all_voxels

print(paste("Number of voxels in total =", dim(MD_cormat)))
MD_cormat <- MD_cormat[common_voxels, common_voxels]
print(paste("Number of voxels used =", dim(MD_cormat)))

diag(MD_cormat) <- 1
gc()

#Eigendecomposition of cormat
start_time <- Sys.time()
eigenvalues <- eigs_sym(MD_cormat, k = 100, opts = list(retvec = FALSE))$values
print("eigenvalues successfully estimated")
end_time <- Sys.time()
eigen_time <- end_time - start_time
paste("Time taken for full eigendecomposition =", eigen_time)
write_rds(eigenvalues, "results/MD_FB_eigenvalues.rds")


rm(MD_cormat)
gc()

#Estimate meff from eigenvalues 
MD_galwey <- meff(eigen = eigenvalues, method = "galwey")
print(paste("meff for mean MD =", MD_galwey))
write_rds(MD_galwey, "results/MD_mean_galwey.rds")

#Clear memory
rm(MD_cormat)
gc()
Sys.sleep(30)
