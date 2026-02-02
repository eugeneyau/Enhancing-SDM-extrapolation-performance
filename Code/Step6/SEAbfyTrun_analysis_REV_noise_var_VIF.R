library(terra)
library(usdm)


# Define study area
Xmini <- 69 #left x-axis
Xmaxi <- 161.6 #right x-axis
Ymini <- 0 #lower y-axis
Ymaxi <- 36 #upper y-axis


# Import world map
Global <- terra::vect("/lustre1/g/sbs_bonebrake/Eugene/SDMin/world.shp")
Global <- terra::crop(Global,ext(Xmini,Xmaxi,Ymini,Ymaxi)) #Crop raster file to define study area // extent(left x-axis, right x-axis, lower y-axis, upper y-axis)
Global <- terra::project(Global,"epsg:6933") 


# Import and stack envi var
shortcut_wd <- "/lustre1/g/sbs_bonebrake/Eugene/SDMin/shortcut"#SDM output result file destination
setwd(shortcut_wd) #Ensure saved to right location

current.climate <- terra::rast("current_climate.tif")

cec <- terra::rast("CEC.tif")
cfvo <- terra::rast("CFVO.tif")
sand <- terra::rast("SAND.tif")
silt <- terra::rast("SILT.tif")
current.climate <- c(current.climate,cec,cfvo,sand,silt)

current.climate <- terra::crop(current.climate,terra::ext(Global))


# Calculate VIF
usdm::vif(current.climate, size=500000)


