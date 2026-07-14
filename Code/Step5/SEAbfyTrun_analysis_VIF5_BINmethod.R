library(terra)
options(scipen = 999)

########################## Settings ##########################
Analysis <- "BINMETH"

filename.base <- "biomod424_BASE_VIF5" #Base distribution
filename <- "bmod424_trunc_VIF5_SOUTH_trun20perc"
#filename <- "bmod424_trunc_VIF5_SOUTH_trun0perc"
#filename <- "bmod424_trunc_VIF5_NORTH_trun20perc"

VirtualEM <- "EMmean"
ProjEM <- "EMmean"

BASE_METRIC <- "TSS"
METRICS <- c("TSS","ROC","POD","POFD")

BASEprojection <- "/proj_Current/proj_Current_"
TSSprojection <- "proj_TSSNA"
projection <- paste0("/",TSSprojection,"/",TSSprojection,"_")

save.filename <- paste0(Analysis,"_VIF5_TruncS20_",ProjEM) # Name the df output file
##############################################################

species.list.out <- list.files(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename))
basepath <- paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename,"/")

########## Common species list ##########

# Dynamically check file existence for all metrics
valid_lists <- list()
for (m in METRICS) {
  poss <- paste0(basepath, species.list.out, projection, species.list.out, "_ensemble_", m, "bin.tif")
  valid_lists[[m]] <- species.list.out[file.exists(poss)]
}

common.list <- Reduce(intersect, valid_lists)

### Clean memory
rm(valid_lists)
gc(full = TRUE)

########## FPOS FNEG analysis ##########

### Create df for documentation
columns <- "Species"
for (m in METRICS) {
  columns <- c(columns, 
               paste0(m,".FP"), paste0(m,".FP.perc"), 
               paste0(m,".FN"), paste0(m,".FN.perc"), 
               paste0(m,".F"),  paste0(m,".F.perc"))
}

diff.perc.df <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(diff.perc.df) <- columns

### Do the analysis

# Automated calculation of error for each species
Accuracy_Fun <- function(i) {
  
  # Load baseline map once per species
  truerange <- terra::rast(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDM_BASE/results_",filename.base,"/",common.list[[i]],BASEprojection,common.list[[i]],"_ensemble_",BASE_METRIC,"bin.tif"))
  truerange <- truerange[VirtualEM]
  
  # Initialize species data row
  sp.data <- list(common.list[[i]])
  
  # Dynamically process each metric map
  for (m in METRICS) {
    map_current <- terra::rast(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename,"/",common.list[[i]],projection,common.list[[i]],"_ensemble_",m,"bin.tif"))[ProjEM]
    
    truerange_cropped <- terra::crop(truerange, terra::ext(map_current))
    
    # Get range
    truerange_size <- as.numeric(terra::global(truerange_cropped, "sum", na.rm = TRUE))
    
    # False Positives
    FPOS <- sum(values(truerange_cropped < map_current), na.rm=TRUE)
    FPOS.perc <- FPOS/truerange_size*100
    
    # False Negatives
    FNEG <- sum(values(truerange_cropped > map_current), na.rm=TRUE)
    FNEG.perc <- FNEG/truerange_size*100
    
    # Total wrong predictions
    FSUM <- FPOS+FNEG
    FSUM.perc <- FPOS.perc+FNEG.perc
    
    # Append to current species list
    sp.data <- c(sp.data, FPOS, FPOS.perc, FNEG, FNEG.perc, FSUM, FSUM.perc)
  }
  
  diff.perc.df[nrow(diff.perc.df)+ 1,] <<- sp.data
}

for(idx in 1:length(common.list)) {
  Accuracy_Fun(idx)}

### Save output df
saveRDS(diff.perc.df,paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",save.filename,".RData"))

### Load output df
save.filename <- "BINMETH_VIF5_TruncS20_EMca"
save.filename <- "BINMETH_VIF5_TruncS20_EMmean"
diff.perc.df <- readRDS(paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",save.filename,".RData"))


########## Visualise results: EM ##########

### False positive number
plot(density(diff.perc.df$TSS.FP,n=10000),xlim = c(-1500,10000),ylim = c(0,0.00038),col="green")
lines(density(diff.perc.df$ROC.FP,n=10000),col="cyan")
lines(density(diff.perc.df$POD.FP,n=10000),col="blue")
lines(density(diff.perc.df$POFD.FP,n=10000),col="purple")

### False negative number
plot(density(diff.perc.df$TSS.FN,n=10000),xlim = c(-2000,10000),ylim = c(0,0.0006),col="green")
lines(density(diff.perc.df$ROC.FN,n=10000),col="cyan")
lines(density(diff.perc.df$POD.FN,n=10000),col="blue")
lines(density(diff.perc.df$POFD.FN,n=10000),col="purple")

### Overall inaccuracy number
plot(density(diff.perc.df$TSS.F,n=10000),xlim = c(-1000,25000),ylim = c(0,0.00015),col="green")
lines(density(diff.perc.df$ROC.F,n=10000),col="cyan")
lines(density(diff.perc.df$POD.F,n=10000),col="blue")
lines(density(diff.perc.df$POFD.F,n=10000),col="purple")



########## Stastistically test results ##########

ks.test(diff.perc.df$ROC.F,
        diff.perc.df$TSS.F,
        alternative = "two.sided",exact = NULL, 
        simulate.p.value=FALSE)
save.filename

ks.test(diff.perc.df$TSS.F, #x
        diff.perc.df$POFD.F, #y
        alternative = "greater",exact = NULL, #Alternative: X>Y
        simulate.p.value=FALSE)
save.filename

ks.test(diff.perc.df$POFD.F, #x
        diff.perc.df$POD.F, #y
        alternative = "greater",exact = NULL, #Alternative: X>Y
        simulate.p.value=FALSE)
save.filename


########## Quick stats ##########

mean(diff.perc.df$EMca.F)
# Trun0 = 3988.066
# Trun20 = 5276.335
# 32.3031
mean(diff.perc.df$EMmean.F)
# Trun0 = 3754.568
# Trun20 = 5910.977
# 57.43428
mean(diff.perc.df$EMwmean.F)
# Trun0 = 3741.895
# Trun20 = 5859.224
# 56.5844
mean(diff.perc.df$EMmedian.F)
# Trun0 = 4110.209
# Trun20 = 7509.839
# 82.71185

# (Trun20-Trun0)/Trun0*100
# mean(c(32.3031,57.43428,56.5844,82.71185))