library(terra)

filename <- "biomod424_BASE_REV"
#filename <- "biomod424_BASE_PC"
wd <- "/lustre1/g/sbs_bonebrake/Eugene/SDM_BASE"


### Get list of modelled species
projection <- "/proj_Current/proj_Current_" #Helps R navigate in your file storing SDM results, need to change if your file structure is different
sp.list <- list.files(paste0(wd,"/results_",filename)) #Extract file names of all files in SDM results folder, should be a list of species names

poss0 <- paste0(wd,"/results_",filename,"/",sp.list,projection,sp.list,"_ensemble_TSSbin.tif") #Generate file path of SDM-produced prediction (.tif raster file) for all species
tf <- file.exists(poss0) #Check which species has SDM output (.tif file)
sp.list <- sp.list[tf] #Generate list of species with SDM output (.tif file), we only want to work on them


# Create df for documenting
columns <- c("Species","Samplesize","Correction") 
sampledf <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(sampledf) <- columns

for (S in 1:length(sp.list)) {
  message("Sampling: ",sp.list[[S]])
  rangemap <- terra::rast(paste0(wd,"/results_",filename,"/",sp.list[[S]],projection,sp.list[[S]],"_ensemble_TSSbin.tif"))
  landarea <- nrow(as.data.frame(rangemap,na.rm=TRUE))
  rangesize <- sum(as.data.frame(rangemap,na.rm=TRUE)[,1])
  maxsamplesize <- floor(rangesize/5)
  range.perc <- rangesize/landarea
  Correction_tf <- "FALSE"
  
  samplesize <- sample(100:250, 1)
  
  if (range.perc < 0.03) {samplesize <- sample(50:100, 1)}
  
  if (range.perc > 0.15) {samplesize <- sample(250:300, 1)}
  
  if (maxsamplesize < samplesize) {
    samplesize <- maxsamplesize
    Correction_tf <- "TRUE"}
  
  sampledf <<- rbind(sampledf, data.frame("Species"=sp.list[[S]], 
                                          "Samplesize"=samplesize,
                                          "Correction"=Correction_tf))
}

write.csv(sampledf,paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMin/",filename,".csv"),row.names = FALSE) #save as csv


#plot(density(sampledf$Samplesize))
