library(dplyr)
library(terra)


#select truncation
Xmini <- 69 #left x-axis, originally 64
Xmaxi <- 161.6 #right x-axis, originally 166
Ymini <- -10 #lower y-axis, originally -13
Ymaxi <- 36 #upper y-axis, originally 37


wd <- "/lustre1/g/sbs_bonebrake/Eugene/SDMin"
setwd(wd)

#Import species data
data <- read.csv("Tropical_Asian_Butterflies_June2025.csv")

colnames(data) <- gsub("decimalLat","decimalLatitude",colnames(data))
colnames(data) <- gsub("decimalLon","decimalLongitude",colnames(data))

data$decimalLatitude <- as.numeric(data$decimalLatitude)
data$decimalLongitude <- as.numeric(data$decimalLongitude)

#Truncate data
data <- subset(data,decimalLatitude>=Ymini)
data <- subset(data,decimalLatitude<=Ymaxi)
data <- subset(data,decimalLongitude>=Xmini)
data <- subset(data,decimalLongitude<=Xmaxi)

data$year <- as.numeric(data$year)
data <- data[! c(1:nrow(data)) %in% which(data$year<1970),]
data[data$lohman_final_genus_species == "Polygonia c-album","lohman_final_genus_species"] <- "Polygonia c album"
data[data$lohman_final_genus_species == "Polygonia c-aureum","lohman_final_genus_species"] <- "Polygonia c aureum"

Count <- data %>% count(lohman_final_genus_species) #Count number of records, Feeds data into count(), saves result in count
Count <- subset(Count, n > 9) #Threshold was 25, changed to 10 for more species

Cleaned.data <- data[which(data$lohman_final_genus_species %in% Count$lohman_final_genus_species),] #Returns  position/index of the value which satisfies given condition

Species.list <- sort(unique(Cleaned.data$lohman_final_genus_species)) #Alphabetically sort species name

Cleaned.data <- Cleaned.data[,c('lohman_final_genus_species',"decimalLongitude","decimalLatitude")]


###Prepare rasters
Global1 <- terra::vect("world.shp")

#Truncate map
Global <- terra::crop(Global1,ext(Xmini,Xmaxi,Ymini,Ymaxi)) #Crop raster file to define study area // extent(left x-axis, right x-axis, lower y-axis, upper y-axis)

Global <- terra::project(Global,"epsg:6933") 
r.rast <- terra::rast(Global, resolution=c(10000,10000))
r <- terra::rasterize(Global, r.rast)

#Environmental rasters
bioclimsource <- "worldclim" # "worldclim" OR "Current"
current_wd <- paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMin/",bioclimsource)

#Select useful bioclimatic variables
file.current.list <-list.files(current_wd)
bio.seq <- grepl("01|04|05|06|12|13|14|15",file.current.list) 
file.current.list <- file.current.list[bio.seq]
current.climate <- c(paste0(current_wd,"/",file.current.list)) # Stacks all previously selected variables together (saved tgt)

current.climate <- terra::rast(current.climate)
#current.climate <- terra::project(current.climate,"epsg:4326")
current.climate <- terra::crop(current.climate,ext(Xmini,Xmaxi,Ymini,Ymaxi))
current.climate <- terra::project(current.climate,"epsg:6933") 
current.climate <- terra::resample(current.climate, r, method="average")

NDVImean <- terra::rast(paste0(current_wd,"/NDVImean.tif"))
names(NDVImean) <- "NDVImean"
NDVImean <- terra::project(NDVImean,"epsg:4326")
NDVImean <- terra::crop(NDVImean,ext(Xmini,Xmaxi,Ymini,Ymaxi))
NDVImean <- terra::project(NDVImean,"epsg:6933")
NDVImean <- terra::resample(NDVImean, r, method="average")

ph <- terra::rast(paste0(current_wd,"/Soilgrids_phh2o.tif"))
names(ph) <- "Soilgrids_soil_pH"
ph <- terra::project(ph,"epsg:4326")
ph <- terra::crop(ph,ext(Xmini,Xmaxi,Ymini,Ymaxi))
ph <- terra::project(ph,"epsg:6933")
ph <- terra::resample(ph, r, method="average")

SOC <- terra::rast(paste0(current_wd,"/Soilgrids_soc.tif"))
names(SOC) <- "Soilgrids_soil_organic_carbon"
SOC <- terra::project(SOC,"epsg:4326")
SOC <- terra::crop(SOC,ext(Xmini,Xmaxi,Ymini,Ymaxi))
SOC <- terra::project(SOC,"epsg:6933")
SOC <- terra::resample(SOC, r, method="average")

NITRO <- terra::rast(paste0(current_wd,"/Soilgrids_nitrogen.tif"))
names(NITRO) <- "Soilgrids_total_nitrogen"
NITRO <- terra::project(NITRO,"epsg:4326")
NITRO <- terra::crop(NITRO,ext(Xmini,Xmaxi,Ymini,Ymaxi))
NITRO <- terra::project(NITRO,"epsg:6933")
NITRO <- terra::resample(NITRO, r, method="average")

c.height <- terra::rast(paste0(current_wd,"/Cheight.tif"))
names(c.height) <- "Canopy_Height"
c.height <- terra::project(c.height,"epsg:4326")
c.height <- terra::crop(c.height,ext(Xmini,Xmaxi,Ymini,Ymaxi))
c.height <- terra::project(c.height,"epsg:6933")
c.height <- terra::resample(c.height, r, method="average")

current.climate <- c(current.climate, NDVImean, c.height, ph, SOC, NITRO)

#########################################################################
#########################################################################
cec <- terra::rast(paste0(current_wd,"/Soilgrids_cec.tif"))
names(cec) <- "Soilgrids_cation_exchange_capacity"
cec <- terra::project(cec,"epsg:4326")
cec <- terra::crop(cec,ext(Xmini,Xmaxi,Ymini,Ymaxi))
cec <- terra::project(cec,"epsg:6933")
cec <- terra::resample(cec, r, method="average")

cfvo <- terra::rast(paste0(current_wd,"/Soilgrids_cfvo.tif"))
names(cfvo) <- "Soilgrids_coarse_fragments"
cfvo <- terra::project(cfvo,"epsg:4326")
cfvo <- terra::crop(cfvo,ext(Xmini,Xmaxi,Ymini,Ymaxi))
cfvo <- terra::project(cfvo,"epsg:6933")
cfvo <- terra::resample(cfvo, r, method="average")

sand <- terra::rast(paste0(current_wd,"/Soilgrids_sand.tif"))
names(sand) <- "Soilgrids_sand"
sand <- terra::project(sand,"epsg:4326")
sand <- terra::crop(sand,ext(Xmini,Xmaxi,Ymini,Ymaxi))
sand <- terra::project(sand,"epsg:6933")
sand <- terra::resample(sand, r, method="average")

silt <- terra::rast(paste0(current_wd,"/Soilgrids_silt.tif"))
names(silt) <- "Soilgrids_soil_silt"
silt <- terra::project(silt,"epsg:4326")
silt <- terra::crop(silt,ext(Xmini,Xmaxi,Ymini,Ymaxi))
silt <- terra::project(silt,"epsg:6933")
silt <- terra::resample(silt, r, method="average")
#########################################################################
#########################################################################



###Save files
shortcut_wd <- "/lustre1/g/sbs_bonebrake/Eugene/SDMin/shortcut"#SDM output result file destination
setwd(shortcut_wd) #Ensure saved to right location

terra::writeRaster(current.climate, "current_climate.tif", overwrite=FALSE)

terra::writeRaster(NDVImean, "NDVImean.tif", overwrite=FALSE)
terra::writeRaster(c.height, "c_height.tif", overwrite=FALSE)
terra::writeRaster(ph, "ph.tif", overwrite=FALSE)
terra::writeRaster(SOC, "SOC.tif", overwrite=FALSE)
terra::writeRaster(NITRO, "NITRO.tif", overwrite=FALSE)

#########################################################################
#########################################################################
terra::writeRaster(cec, "CEC.tif", overwrite=FALSE)
terra::writeRaster(cfvo, "CFVO.tif", overwrite=FALSE)
terra::writeRaster(sand, "SAND.tif", overwrite=FALSE)
terra::writeRaster(silt, "SILT.tif", overwrite=FALSE)
#########################################################################
#########################################################################

#terra::writeRaster(dens.ras, "dens_ras.tif", overwrite=FALSE)

write.csv(Cleaned.data, "Cleaned_published_data.csv", row.names = FALSE)

