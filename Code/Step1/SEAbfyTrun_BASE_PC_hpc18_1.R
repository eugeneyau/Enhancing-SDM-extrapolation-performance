
library(biomod2)
library(maptools)
library(BBmisc)
library(dplyr)
library(terra)

########################## Settings ##########################
m <- 1:370
#m <- 371:740
#m <- 741:1110
#m <- 1111:1481

PAnb=4

#Define study area
Xmini <- 69 #left x-axis
Xmaxi <- 161.6 #right x-axis
Ymini <- -10 #lower y-axis
Ymaxi <- 36 #upper y-axis

#Model run name
datetime <- "biomod424_BASE_PC"
##############################################################

#Set directory to save results
model.files.dir <- paste0("results_", datetime)
dir.create(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDM_BASE/",model.files.dir)) 
species.list.dir <- paste0("species_list_", datetime)
dir.create(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDM_BASE/",species.list.dir)) 

###Import input data
wd <- "/lustre1/g/sbs_bonebrake/Eugene/SDMin"
setwd(wd)

#Full map
Global <- terra::vect("world.shp")
Global <- terra::crop(Global,ext(Xmini,Xmaxi,Ymini,Ymaxi)) #Crop raster file to define study area // extent(left x-axis, right x-axis, lower y-axis, upper y-axis)
Global <- terra::project(Global,"epsg:6933") 
r.rast <- terra::rast(Global, resolution=c(10000,10000))
r <- terra::rasterize(Global, r.rast)

#Set shortcut wd
shortcut_wd <- "/lustre1/g/sbs_bonebrake/Eugene/SDMin/shortcut"#SDM output result file destination
setwd(shortcut_wd) #Ensure saved to right location

#Import species data
Cleaned.data <- read.csv("Cleaned_published_data_20.csv")

Species.list <- unique(Cleaned.data$lohman_final_genus_species) #Enters species name

#Cleaning Species.list
Species.list[Species.list == "Polygonia c-album"] <- "Polygonia c album" #R cannot handle hyphen, have to change to underscore 
Species.list[Species.list == "Polygonia c-aureum"] <- "Polygonia c aureum"
Species.list <- Species.list[ !Species.list == 'not present'] # Remove the "not present" value in Species.list
Species.list <- sort(Species.list) # Alphabetically sort list

###Get pending species list
splist <- list.files(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDM_BASE/",species.list.dir))
splist <- gsub("\\.", " ", splist)
finished <- sub("_Finally_Finished.txt", "", splist)
insufficient.rec <- sub("_Insufficient_Records_Failed.txt", "", splist)
all.done <- c(finished,insufficient.rec)
bool <- Species.list %in% all.done
pending_sp_list <- Species.list[!bool]

Species.list <- sort(pending_sp_list)

#m <- 1:ceiling(length(Species.list)/4)
#m <- floor(length(Species.list)/4):ceiling(length(Species.list)/2)
#m <- floor(length(Species.list)/2):ceiling(length(Species.list)*3/4)
#m <- floor(length(Species.list)*3/4):length(Species.list)
###

#Import envi var
#readRDS("ALLenv_PCAsummary.RData") # PC1-7 explains >95% total variance
current.climate <- terra::rast("current_env_PC.tif")[[1:7]]

######

save_wd <- paste0("/lustre1/g/sbs_bonebrake/Eugene/SDM_BASE/",model.files.dir) #SDM output result file destination
setwd(save_wd) #Ensure saved to right location

ModelVirtualSP <- function(m) {
  try({
    message("Species", Species.list[[m]]) #Let me know which SDM the program is running
    error.msg <- NULL
    SDM_data <- subset(Cleaned.data,lohman_final_genus_species == Species.list[m]) #(Subsets/draws out data of single species for analysis)
    SDM_data <- SDM_data[!is.na(SDM_data$decimalLatitude),] #Remove negative/NA data
    myRespXY <- SDM_data[,c("decimalLongitude","decimalLatitude")] #Draws out latitude and longitude for each records
    Projected.co <- terra::vect(myRespXY, geom=c("decimalLongitude","decimalLatitude"), crs="epsg:4326") #Converts to point file and give CRS, for it to be identified as GIS stuff
    T.Projected.co <- terra::project(Projected.co,"epsg:6933")
    T.Projected.co.raster <- terra::rasterize(T.Projected.co, r, vals=1) #Use base raster to make raster which has data only at occurrence point
    T.Projected.co.raster[T.Projected.co.raster>=1] <- 1 #Limit the max. number of data points to 1 in every 10km x 10km grid square 
    
    #Create a file "df" to record 10km x 10km coordinates of grid squares with >=1 data points
    df <- terra::crds(T.Projected.co.raster, na.rm=TRUE) 
    
    if (nrow(df) > 24) { #Check if we have more than 25 EFFECTIVE records (10km x 10km grid squares)
      
      myBiomodOption <- BIOMOD_ModelingOptions() #Default
      
      vectXY <- terra::vect(df, crs="epsg:6933")
      
      #Assign value 1 to presence
      values(vectXY) <- data.frame(rec=rep(1, nrow(df)))
      
      #Format data
      myBiomodData <- BIOMOD_FormatingData(resp.var = vectXY, 
                                           expl.var = current.climate, 
                                           resp.name = Species.list[m], 
                                           PA.nb.rep = PAnb, 
                                           PA.nb.absences = nrow(df), 
                                           PA.strategy = "random", 
                                           filter.raster = TRUE)  
      
      #Fit SDMs
      myBiomodModelOut <- BIOMOD_Modeling( 
        bm.format = myBiomodData, 
        models = c("ANN","CTA","GLM","GBM","MARS","MAXNET","RF","XGBOOST"), 
        bm.options = myBiomodOption, 
        CV.strategy = "random",
        CV.nb.rep = 4,
        CV.perc = 0.9, 
        prevalence = 0.5, #Default setting on model weighting
        metric.eval = c("TSS"), #Model evaluation methods
        scale.models = TRUE, #models predictions scaled with a binomial GLM or not *****Might think about it*****
        CV.do.full.models = FALSE, 
        modeling.id = paste("FirstModeling",sep=""), #Name of model (Species name_Firstmodeling)
        var.import = 5, #number of permutations to be done for each variable to estimate variable importance
        seed.val = 66) 
      
      #csv file for performance of each models built
      write.csv(get_evaluations(myBiomodModelOut),
                file.path(paste(gsub(" |_",".",Species.list[m]),"/",Species.list[m],"_formal_models_evaluation.csv", sep="")))
      #csv file for importance of each variables
      write.csv(get_variables_importance(myBiomodModelOut),
                file.path(paste(gsub(" |_",".",Species.list[m]),"/",Species.list[m],"_formal_models_variables_importance.csv", sep="")))
      
      #Merging the 15 models into one model
      myBiomodEM <- try({BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                                 models.chosen = "all",
                                                 em.by = "all", 
                                                 metric.select = c("TSS"),
                                                 metric.select.thresh = NULL, 
                                                 em.algo = c('EMca','EMmean','EMmedian'))}) 
      
      #If successfully merged single models: Each of the single models individually project
      myBiomodProj <- BIOMOD_Projection(
        bm.mod = myBiomodModelOut,
        new.env = current.climate,
        proj.name = "Current",
        models.chosen = "all",
        metric.binary = "TSS",
        compress = "xz",
        build.clamping.mask = T,
        output.format = ".grd")
      
      #Then get ensemble projection
      if (length(myBiomodEM@em.models_kept) > 1) {
        
        myBiomodEF <- BIOMOD_EnsembleForecasting(
          bm.em = myBiomodEM,
          bm.proj = myBiomodProj,
          models.chosen = 'all',
          metric.binary = "TSS")
        
        write.table("Finished", file.path(paste("/lustre1/g/sbs_bonebrake/Eugene/SDM_BASE/",species.list.dir,"/",Species.list[m],"_Finally_Finished.txt", sep="")))
        
      } else if (length(myBiomodEM@em.models_kept) == 1) {
        
        myBiomodEF <- try({BIOMOD_EnsembleForecasting(
          bm.em = myBiomodEM,
          bm.proj = myBiomodProj,
          models.chosen = 'all',
          metric.binary = "TSS")})
        
        write.table("Finished", file.path(paste("/lustre1/g/sbs_bonebrake/Eugene/SDM_BASE/",species.list.dir,"/",Species.list[m],"_Single_Model_Finished.txt", sep="")))
        
      } else { write.table("Failed", file.path(paste("/lustre1/g/sbs_bonebrake/Eugene/SDM_BASE/",species.list.dir,"/",Species.list[m],"_EM_Failed.txt", sep="")))}
      
      write.table("Finished", file.path(paste(gsub(" |_",".",Species.list[m]),"/",Species.list[m],"_Finished.txt", sep="")))
      
    } else {
      write.table("Failed", file.path(paste("/lustre1/g/sbs_bonebrake/Eugene/SDM_BASE/",species.list.dir,"/",Species.list[m],"_Insufficient_Records_Failed.txt", sep="")))
    }
  })}

library(parallel)
mclapply(m, ModelVirtualSP, mc.cores=18)