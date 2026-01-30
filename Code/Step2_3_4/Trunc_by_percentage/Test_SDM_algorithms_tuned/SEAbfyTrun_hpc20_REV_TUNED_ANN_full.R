
library(biomod2)
library(maptools)
library(BBmisc)
library(dplyr)
library(terra)

########################## Settings ##########################
#m <- 1:984

PAnb=5
CVnb=5

#Define study area
Xmini <- 69 #left x-axis
Xmaxi <- 161.6 #right x-axis
Ymini <- 0 #lower y-axis
Ymaxi <- 36 #upper y-axis

#Select truncation
TrunAmount <- 0 #0/0.1/0.2/0.3, portion of data to truncate
Trunarea <- "SOUTH" # SOUTH/NORTH

#Model run name
SDMalgo <- "ANN" #"ANN","CTA","GLM","GBM","MARS","RF" // No tuning for "MAXNET","XGBOOST"
runid <- "bmod424_trunc_rev_tuned"

#Virtual species location
vspdir <- "biomod424_BASE_REV"
projection <- "/proj_Current/proj_Current_" #Helps R navigate in your file storing SDM results
EMmethod <- "EMmean"
##############################################################

#Set directory to save results
TrunPrec <- TrunAmount*100
model.files.dir <- paste0("results_",runid,"_",SDMalgo,"_",Trunarea,"_trun",TrunPrec,"perc")
dir.create(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/",model.files.dir)) 
species.list.dir <- paste0("species_list_",runid,"_",SDMalgo,"_",Trunarea,"_trun",TrunPrec,"perc")
dir.create(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/",species.list.dir)) 

###Import input data
wd <- "/lustre1/g/sbs_bonebrake/Eugene/SDMin"
setwd(wd)

#Full map
Global <- terra::vect("world.shp")
Global <- terra::crop(Global,ext(Xmini,Xmaxi,Ymini,Ymaxi)) #Crop raster file to define study area // extent(left x-axis, right x-axis, lower y-axis, upper y-axis)
Global <- terra::project(Global,"epsg:6933") 

###Get species list and sample size
sampledf <- read.csv(paste0(vspdir,".csv"))
Species.list <- sort(sampledf$Species)

#Get pending species list
splist <- list.files(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/",species.list.dir))
finished <- sub("_Finally_Finished.txt", "", splist)
insufficient.rec <- sub("_Insufficient_Records_Failed.txt", "", splist)
all.done <- c(finished,insufficient.rec)
bool <- Species.list %in% all.done
pending_sp_list <- Species.list[!bool]

Species.list <- sort(pending_sp_list)

m <- 1:length(Species.list)

###

#Import envi var
shortcut_wd <- "/lustre1/g/sbs_bonebrake/Eugene/SDMin/shortcut"#SDM output result file destination
setwd(shortcut_wd) #Ensure saved to right location

current.climate <- terra::rast("current_climate.tif")
current.climate <- terra::crop(current.climate,ext(Global))

######

save_wd <- paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/",model.files.dir) #SDM output result file destination
setwd(save_wd) #Ensure saved to right locationsetwd(save_wd) #Ensure saved to right location


MyBiomodSF <- function(m) {
  try({
    #Get virtual distribution
    rangemap <- terra::rast(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDM_BASE/results_",vspdir,"/",Species.list[m],projection,Species.list[m],"_ensemble_TSSbin.tif"))
    rangemap <- rangemap[EMmethod]
    rangemap <- terra::crop(rangemap,ext(Global))
    Obs_df <- terra::as.data.frame(rangemap, xy=TRUE)
    colnames(Obs_df) <- c("x","y","Suitability")
    Obs_df <- subset(Obs_df, Suitability == 1)
    
    #Sample data
    samplesize <- subset(sampledf,Species == Species.list[m])[,2]
    maxsamplesize <- floor(nrow(Obs_df)/5)
    if (maxsamplesize < samplesize) {samplesize <- maxsamplesize}
    
    df <- Obs_df[sample(1:nrow(Obs_df), samplesize, replace=FALSE),]
    df <- df[,-3]
    
    #Get vector of sampled data
    vectXY <- terra::vect(df,geom=c("x","y"),crs="epsg:6933")
    
    if (nrow(df) > 20) { #GBM works with at least 21 records
      
      values(vectXY) <- data.frame(rec=rep(1, nrow(df)))
      
      #Format data
      myBiomodData <- BIOMOD_FormatingData(resp.var = vectXY, 
                                           expl.var = current.climate, 
                                           resp.name = Species.list[m], 
                                           PA.nb.rep = PAnb, 
                                           PA.nb.absences = nrow(df), 
                                           PA.strategy = "random", 
                                           filter.raster = TRUE)  
      
      #Model tuning
      bm.tuning <- BIOMOD_Tuning(
        models = SDMalgo,
        metric.eval = "TSS",
        bm.format = myBiomodData, 
        ME.env = current.climate)
      myBiomodOption <- bm.tuning$models.options
      
      #SDM fitting
      myBiomodModelOut <- BIOMOD_Modeling( 
        bm.format = myBiomodData, 
        models = SDMalgo,
        bm.options = myBiomodOption, 
        CV.strategy = "random",
        CV.nb.rep = CVnb,
        CV.perc = 0.9, 
        prevalence = 0.5, #Default setting on model weighting
        metric.eval = "TSS", #Model evaluation methods
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
      
      
      ########## No TSS filtering ##########
      
      #Ensemble modelling
      myBiomodEMNA <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                              models.chosen = "all",
                                              em.by = "all", 
                                              metric.select = c("TSS"),
                                              metric.select.thresh = NULL,
                                              em.algo = c('EMmean', 'EMca')) 
      
      #Single models individually estimates habitat suitability (Projection)
      myBiomodProjNA <- BIOMOD_Projection(
        bm.mod = myBiomodModelOut,
        new.env = current.climate,
        proj.name = "TSSNA",
        models.chosen = "all",
        metric.binary = "TSS",
        compress = "xz",
        build.clamping.mask = T, # Very important
        output.format = ".grd")
      
      myBiomodEF <- BIOMOD_EnsembleForecasting(
        bm.em = myBiomodEMNA,
        bm.proj = myBiomodProjNA,
        models.chosen = 'all',
        metric.binary = "TSS")
      
      write.table("Finished", file.path(paste(gsub(" |_",".",Species.list[m]),"/",Species.list[m],"_Finished.txt", sep="")))
      write.table("Finished", file.path(paste("/lustre1/g/sbs_bonebrake/Eugene/SDMout/",species.list.dir,"/",Species.list[m],"_Finally_Finished.txt", sep="")))
      
    } else {
      write.table("Failed", 
                  file.path(paste("/lustre1/g/sbs_bonebrake/Eugene/SDMout/",species.list.dir,"/",Species.list[m],"_Insufficient_Records_Failed.txt", sep="")))
    }
  })}

library(parallel)
mclapply(m, MyBiomodSF, mc.cores=20)