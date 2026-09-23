#The correlation matrix has been  estimated in another parallel script called Correlations_parallel_doMPI.R
#The output is opened in this script and used to project a correlation plot
library("corrplot")
library("tidyverse")

#FA----------------------------------------------------------------------------

#Create image device
png("plots/correlations/dHCP_FA_image_all_voxels.png", units = "px",
    width=10000, height=10000, bg = "transparent", res = 2500)
par(mar = c(0, 0, 0, 0))

#Source skifti file
start_time <- Sys.time()
FA_cormat <- readRDS("results/FA_cormat_result.rds") %>% as.matrix()
end_time <- Sys.time()
print("Time taken to load correlation matrix rds file:")
print(end_time - start_time)

#Plot
start_time <- Sys.time()
image(FA_cormat, Rowv = NA, Colv = NA, col = COL2('RdBu', 200), 
      zlim = c(-1,1), axes = FALSE)

end_time <- Sys.time()
print("Time taken to project correlation plot:")
print(end_time - start_time)

#Write image to disk
start_time <- Sys.time()
dev.off()
end_time <- Sys.time()
print("Time taken to save plot as PNG:")
print(end_time - start_time)

#Remove FA correlation matrix from memory
rm(FA_cormat)