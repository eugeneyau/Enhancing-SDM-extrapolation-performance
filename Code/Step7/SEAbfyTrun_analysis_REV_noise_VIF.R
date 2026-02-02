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






### Results ###
> usdm::vif(current.climate, size=50000)
                         Variables        VIF
1                     WorldClim_bio01 232.121249
2                     WorldClim_bio04  20.434554
3                     WorldClim_bio05  72.308180
4                     WorldClim_bio06 135.618886
5                     WorldClim_bio12  27.128650
6                     WorldClim_bio13  14.893679
7                     WorldClim_bio14   6.490471
8                     WorldClim_bio15   6.732263
9                            NDVImean   3.678139
10                      Canopy_Height   4.922268
11                  Soilgrids_soil_pH   6.450514
12      Soilgrids_soil_organic_carbon   2.993907
13           Soilgrids_total_nitrogen   2.555104
14 Soilgrids_cation_exchange_capacity   2.053696
15         Soilgrids_coarse_fragments   3.491813
16                     Soilgrids_sand   4.166125
17                Soilgrids_soil_silt   4.338276

> usdm::vif(current.climate, size=100000)
                         Variables        VIF
1                     WorldClim_bio01 233.234160
2                     WorldClim_bio04  20.121835
3                     WorldClim_bio05  73.070657
4                     WorldClim_bio06 134.742902
5                     WorldClim_bio12  26.573440
6                     WorldClim_bio13  14.594325
7                     WorldClim_bio14   6.262962
8                     WorldClim_bio15   6.651498
9                            NDVImean   3.718006
10                      Canopy_Height   4.958766
11                  Soilgrids_soil_pH   6.531421
12      Soilgrids_soil_organic_carbon   2.974842
13           Soilgrids_total_nitrogen   2.571805
14 Soilgrids_cation_exchange_capacity   2.047438
15         Soilgrids_coarse_fragments   3.507846
16                     Soilgrids_sand   4.198958
17                Soilgrids_soil_silt   4.402932


> usdm::vif(current.climate, size=200000)
                         Variables        VIF
1                     WorldClim_bio01 232.780850
2                     WorldClim_bio04  20.164935
3                     WorldClim_bio05  73.002528
4                     WorldClim_bio06 134.582966
5                     WorldClim_bio12  26.553036
6                     WorldClim_bio13  14.605866
7                     WorldClim_bio14   6.315668
8                     WorldClim_bio15   6.642753
9                            NDVImean   3.721439
10                      Canopy_Height   4.963390
11                  Soilgrids_soil_pH   6.501751
12      Soilgrids_soil_organic_carbon   2.994850
13           Soilgrids_total_nitrogen   2.571406
14 Soilgrids_cation_exchange_capacity   2.044012
15         Soilgrids_coarse_fragments   3.496890
16                     Soilgrids_sand   4.214937
17                Soilgrids_soil_silt   4.417435
Warning message:
  [spatSample] fewer values returned than requested 


> usdm::vif(current.climate, size=500000)
                         Variables        VIF
1                     WorldClim_bio01 232.780850
2                     WorldClim_bio04  20.164935
3                     WorldClim_bio05  73.002528
4                     WorldClim_bio06 134.582966
5                     WorldClim_bio12  26.553036
6                     WorldClim_bio13  14.605866
7                     WorldClim_bio14   6.315668
8                     WorldClim_bio15   6.642753
9                            NDVImean   3.721439
10                      Canopy_Height   4.963390
11                  Soilgrids_soil_pH   6.501751
12      Soilgrids_soil_organic_carbon   2.994850
13           Soilgrids_total_nitrogen   2.571406
14 Soilgrids_cation_exchange_capacity   2.044012
15         Soilgrids_coarse_fragments   3.496890
16                     Soilgrids_sand   4.214937
17                Soilgrids_soil_silt   4.417435
Warning messages:
  1: [spatSample] requested sample size is larger than the number of cells 
2: [spatSample] more non-NA cells requested than available 






