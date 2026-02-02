library(terra)

# Set input wd
shortcut_wd <- "/lustre1/g/sbs_bonebrake/Eugene/SDMin/shortcut"#SDM output result file destination
setwd(shortcut_wd) #Ensure saved to right location

# Import envi var
current.climate <- terra::rast("current_climate.tif")

# Extract data and run PCA
matrix <- terra::values(current.climate,na.rm=TRUE)
pca <- prcomp(matrix, center = TRUE, scale = TRUE)
PCAsum <- summary(pca)
summary(pca)

# Extract PCs
PCenv <- predict(current.climate, pca)

# Save PCs
terra::writeRaster(PCenv,"current_env_PC.tif",overwrite=FALSE)

# Save and load PCA summary (proportion of variance explained)
saveRDS(PCAsum, "ALLenv_PCAsummary.RData")
readRDS("ALLenv_PCAsummary.RData")
