
library(terra)
options(scipen = 999)

# Define the list of algorithms
algo <- sort(c("ANN","CTA","GLM","MARS","MAXNET","RF","XGBOOST","GBM","tuned_ANN","tuned_CTA","tuned_GLM","tuned_GBM","tuned_MARS","tuned_RF"))
#algo <- sort(c("ANN","CTA","GLM","MARS","MAXNET","RF","XGBOOST","GBM","tuned_ANN","tuned_CTA","tuned_GLM","tuned_GBM","tuned_0interact_MARS","tuned_RF"))

# Specify filenames
RunID <- "VIF5"
fullarea <- "SOUTH"
trunarea <- "SOUTH"
Virtualsp_dir <- "biomod424_BASE_VIF5" #Virtual distribution
All_algo_full_dir <- "bmod424_trunc_VIF5_SOUTH_trun0perc" #All algo FULL
All_algo_trun_dir <- "bmod424_trunc_VIF5_SOUTH_trun20perc" #All algo TRUN
Best_algo_full_dir <- "NA" #Best algo FULL
Best_algo_trun_dir <- "bmod424_trunc_VIF5_bestalgo_SOUTH_trun20perc" #Best algo TRUN

# List filenames
filenames <- c(
  Virtualsp_dir,
  All_algo_full_dir,
  All_algo_trun_dir,
  sort(c(paste0("bmod424_trunc_",RunID,"_", algo, "_",fullarea,"_trun0perc"),
         paste0("bmod424_trunc_",RunID,"_", algo, "_",trunarea,"_trun20perc")))
)

if (!Best_algo_full_dir=="NA") {filenames <- c(filenames, Best_algo_full_dir)}
if (!Best_algo_trun_dir=="NA") {filenames <- c(filenames, Best_algo_trun_dir)}

# Settings
VirtualEM <- "EMmean"
EM <- "EMca"
AlgoEM <- "EMca"
TSSprojection <- "proj_TSSNA"

save.filename <- paste0("Algo_VIF5_",EM,"_",TSSprojection,"_trun",trunarea) # Name the df output file

projection <- paste0("/",TSSprojection,"/",TSSprojection,"_")
BASEprojection <- "/proj_Current/proj_Current_"
base_dir <- "/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_"

# Function to find common species across all directories
get_common_list <- function(f) {
  basepath <- paste0(base_dir,filenames[f])
  species.list.out <- list.files(basepath)
  poss <- paste0(basepath,"/",species.list.out,projection,species.list.out,"_ensemble_TSSbin.tif")
  tf <- file.exists(poss)
  splist <- species.list.out[tf]
  if (f==2) {common.list <<- splist}
  common.list <<- intersect(common.list, splist) 
}

lapply(2:length(filenames),get_common_list)

# Create data frame for documentation
multialgo <- c("All.f.F", "All.f.F.perc", "All.f.FN", "All.f.FN.perc", "All.f.FP", "All.f.FP.perc",
               "All.t.F", "All.t.F.perc", "All.t.FN", "All.t.FN.perc", "All.t.FP", "All.t.FP.perc")

if (!Best_algo_full_dir=="NA") {multialgo <- c(multialgo, "Best.f.F", "Best.f.F.perc", "Best.f.FN", "Best.f.FN.perc", "Best.f.FP", "Best.f.FP.perc")}
if (!Best_algo_trun_dir=="NA") {multialgo <- c(multialgo, "Best.t.F", "Best.t.F.perc", "Best.t.FN", "Best.t.FN.perc", "Best.t.FP", "Best.t.FP.perc")}

columns <- c("Species", multialgo,
             sort(c(paste0(algo, ".f.F"), paste0(algo, ".f.F.perc"),
                    paste0(algo, ".f.FN"), paste0(algo, ".f.FN.perc"),
                    paste0(algo, ".f.FP"), paste0(algo, ".f.FP.perc"),
                    paste0(algo, ".t.F"), paste0(algo, ".t.F.perc"),
                    paste0(algo, ".t.FN"), paste0(algo, ".t.FN.perc"),
                    paste0(algo, ".t.FP"), paste0(algo, ".t.FP.perc")))
)  

diff.perc.df <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(diff.perc.df) <- columns

# Function to calculate accuracy metrics for a given species and algorithm
Accuracy_Fun <- function(i) {
  
  # Load maps
  map_all_f <- terra::rast(paste0(base_dir,All_algo_full_dir,"/",common.list[[i]],projection,common.list[[i]],"_ensemble_TSSbin.tif"))
  map_all_f <- map_all_f[EM]
  map_all_t <- terra::rast(paste0(base_dir,All_algo_trun_dir,"/",common.list[[i]],projection,common.list[[i]],"_ensemble_TSSbin.tif"))
  map_all_t <- map_all_t[EM]
  
  truerange <- terra::rast(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDM_BASE/results_",Virtualsp_dir,"/",common.list[[i]],BASEprojection,common.list[[i]],"_ensemble_TSSbin.tif"))
  truerange <- truerange[VirtualEM]
  truerange <- terra::crop(truerange, terra::ext(map_all_f))
  
  truerange_size <- subset(terra::freq(truerange), value == 1)[, 3]
  
  # False Positives
  FPOS_f <- sum(values(truerange<map_all_f), na.rm=TRUE)
  FPOS.perc_f <- FPOS_f/truerange_size*100
  FPOS_t <- sum(values(truerange<map_all_t), na.rm=TRUE)
  FPOS.perc_t <- FPOS_t/truerange_size*100
  
  # False Negatives
  FNEG_f <- sum(values(truerange>map_all_f), na.rm=TRUE)
  FNEG.perc_f <- FNEG_f/truerange_size*100
  FNEG_t <- sum(values(truerange>map_all_t), na.rm=TRUE)
  FNEG.perc_t <- FNEG_t/truerange_size*100
  
  # Total wrong predictions
  FSUM_f <- FPOS_f+FNEG_f
  FSUM.perc_f <- FPOS.perc_f+FNEG.perc_f
  FSUM_t <- FPOS_t+FNEG_t
  FSUM.perc_t <- FPOS.perc_t+FNEG.perc_t
  
  sp.data <- list(common.list[[i]],
                  FSUM_f,FSUM.perc_f,FNEG_f,FNEG.perc_f,FPOS_f,FPOS.perc_f,
                  FSUM_t,FSUM.perc_t,FNEG_t,FNEG.perc_t,FPOS_t,FPOS.perc_t)
  
  
  if (!Best_algo_full_dir=="NA") {
    # Load maps
    map_f <- terra::rast(paste0(base_dir,Best_algo_full_dir,"/",common.list[[i]],projection,common.list[[i]],"_ensemble_TSSbin.tif"))
    map_f <- map_f[EM]
    
    # False Positives
    FPOS_f <- sum(values(truerange<map_f), na.rm=TRUE)
    FPOS.perc_f <- FPOS_f/truerange_size*100
    
    # False Negatives
    FNEG_f <- sum(values(truerange>map_f), na.rm=TRUE)
    FNEG.perc_f <- FNEG_f/truerange_size*100
    
    # Total wrong predictions
    FSUM_f <- FPOS_f+FNEG_f
    FSUM.perc_f <- FPOS.perc_f+FNEG.perc_f
    
    sp.data <- c(sp.data, list(FSUM_f,FSUM.perc_f,FNEG_f,FNEG.perc_f,FPOS_f,FPOS.perc_f))}
  
  
  if (!Best_algo_trun_dir=="NA") {
    # Load maps
    map_t <- terra::rast(paste0(base_dir,Best_algo_trun_dir,"/",common.list[[i]],projection,common.list[[i]],"_ensemble_TSSbin.tif"))
    map_t <- map_t[EM]
    
    # False Positives
    FPOS_t <- sum(values(truerange<map_t), na.rm=TRUE)
    FPOS.perc_t <- FPOS_t/truerange_size*100
    
    # False Negatives
    FNEG_t <- sum(values(truerange>map_t), na.rm=TRUE)
    FNEG.perc_t <- FNEG_t/truerange_size*100
    
    # Total wrong predictions
    FSUM_t <- FPOS_t+FNEG_t
    FSUM.perc_t <- FPOS.perc_t+FNEG.perc_t
    
    sp.data <- c(sp.data, list(FSUM_t,FSUM.perc_t,FNEG_t,FNEG.perc_t,FPOS_t,FPOS.perc_t))}
  
  
  for (k in 1:length(algo)) {
    
    # Get filenames
    Algo_full_dir <- paste0(base_dir,"bmod424_trunc_",RunID,"_",algo[k],"_",fullarea,"_trun0perc")
    Algo_trun_dir <- paste0(base_dir,"bmod424_trunc_",RunID,"_",algo[k],"_",trunarea,"_trun20perc")
    
    # Load maps
    map_f <- terra::rast(paste0(Algo_full_dir,"/",common.list[[i]],projection,common.list[[i]],"_ensemble_TSSbin.tif"))
    map_f <- map_f[AlgoEM]
    map_t <- terra::rast(paste0(Algo_trun_dir,"/",common.list[[i]],projection,common.list[[i]],"_ensemble_TSSbin.tif"))
    map_t <- map_t[AlgoEM]
    
    # False Positives
    FPOS_f <- sum(values(truerange<map_f), na.rm=TRUE)
    FPOS.perc_f <- FPOS_f/truerange_size*100
    FPOS_t <- sum(values(truerange<map_t), na.rm=TRUE)
    FPOS.perc_t <- FPOS_t/truerange_size*100
    
    # False Negatives
    FNEG_f <- sum(values(truerange>map_f), na.rm=TRUE)
    FNEG.perc_f <- FNEG_f/truerange_size*100
    FNEG_t <- sum(values(truerange>map_t), na.rm=TRUE)
    FNEG.perc_t <- FNEG_t/truerange_size*100
    
    # Total wrong predictions
    FSUM_f <- FPOS_f+FNEG_f
    FSUM.perc_f <- FPOS.perc_f+FNEG.perc_f
    FSUM_t <- FPOS_t+FNEG_t
    FSUM.perc_t <- FPOS.perc_t+FNEG.perc_t
    
    sp.data <- c(sp.data, list(FSUM_f,FSUM.perc_f,FNEG_f,FNEG.perc_f,FPOS_f,FPOS.perc_f,
                               FSUM_t,FSUM.perc_t,FNEG_t,FNEG.perc_t,FPOS_t,FPOS.perc_t))
    
  }
  
  diff.perc.df[nrow(diff.perc.df)+1,] <<- sp.data
  
}

i=1:length(common.list)

lapply(i,Accuracy_Fun)


# Save the results
saveRDS(diff.perc.df,paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",save.filename,".RData"))



########## Visualise results: Algo ##########

### Load output df
diff.perc.df <- readRDS(paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",save.filename,".RData"))

### False positive percentage
plot(density(diff.perc.df$All.f.FP.perc,n=10000),xlim = c(-10,150))
lines(density(diff.perc.df$All.t.FP.perc,n=10000),lty=3)
lines(density(diff.perc.df$GLM.f.FP.perc,n=10000),col="#332288")
lines(density(diff.perc.df$GLM.t.FP.perc,n=10000),col="#332288",lty=3)
lines(density(diff.perc.df$CTA.f.FP.perc,n=10000),col="#88CCEE")
lines(density(diff.perc.df$CTA.t.FP.perc,n=10000),col="#88CCEE",lty=3)
lines(density(diff.perc.df$MARS.f.FP.perc,n=10000),col="#44AA99")
lines(density(diff.perc.df$MARS.t.FP.perc,n=10000),col="#44AA99",lty=3)
lines(density(diff.perc.df$MAXNET.f.FP.perc,n=10000),col="#999933")
lines(density(diff.perc.df$MAXNET.t.FP.perc,n=10000),col="#999933",lty=3)
lines(density(diff.perc.df$XGBOOST.f.FP.perc,n=10000),col="#DDCC77")
lines(density(diff.perc.df$XGBOOST.t.FP.perc,n=10000),col="#DDCC77",lty=3)
lines(density(diff.perc.df$RF.f.FP.perc,n=10000),col="#CC6677")
lines(density(diff.perc.df$RF.t.FP.perc,n=10000),col="#CC6677",lty=3)
lines(density(diff.perc.df$ANN.f.FP.perc,n=10000),col="#882255")
lines(density(diff.perc.df$ANN.t.FP.perc,n=10000),col="#882255",lty=3)
lines(density(diff.perc.df$GBM.f.FP.perc,n=10000),col="#AA4499")
lines(density(diff.perc.df$GBM.t.FP.perc,n=10000),col="#AA4499",lty=3)

plot(density(diff.perc.df$All.t.FP.perc,n=10000),xlim = c(-10,150),lty=3)
lines(density(diff.perc.df$GLM.t.FP.perc,n=10000),col="#332288",lty=3)
lines(density(diff.perc.df$CTA.t.FP.perc,n=10000),col="#88CCEE",lty=3)
lines(density(diff.perc.df$MARS.t.FP.perc,n=10000),col="#44AA99",lty=3)
lines(density(diff.perc.df$MAXNET.t.FP.perc,n=10000),col="#999933",lty=3)
lines(density(diff.perc.df$XGBOOST.t.FP.perc,n=10000),col="#DDCC77",lty=3)
lines(density(diff.perc.df$RF.t.FP.perc,n=10000),col="#CC6677",lty=3)
lines(density(diff.perc.df$ANN.t.FP.perc,n=10000),col="#882255",lty=3)
lines(density(diff.perc.df$GBM.t.FP.perc,n=10000),col="#AA4499",lty=3)

plot(density(diff.perc.df$GLM.t.FP.perc,n=10000),col="#332288",xlim=c(-10,170),ylim=c(0,0.020))
lines(density(diff.perc.df$tuned_GLM.t.FP.perc,n=10000),col="#332288",lty=2)
lines(density(diff.perc.df$CTA.t.FP.perc,n=10000),col="#88CCEE")
lines(density(diff.perc.df$tuned_CTA.t.FP.perc,n=10000),col="#88CCEE",lty=2)
lines(density(diff.perc.df$MARS.t.FP.perc,n=10000),col="#44AA99")
lines(density(diff.perc.df$tuned_MARS.t.FP.perc,n=10000),col="#44AA99",lty=2)
lines(density(diff.perc.df$RF.t.FP.perc,n=10000),col="#CC6677")
lines(density(diff.perc.df$tuned_RF.t.FP.perc,n=10000),col="#CC6677",lty=2)
lines(density(diff.perc.df$ANN.t.FP.perc,n=10000),col="#882255")
lines(density(diff.perc.df$tuned_ANN.t.FP.perc,n=10000),col="#882255",lty=2)
lines(density(diff.perc.df$GBM.t.FP.perc,n=10000),col="#AA4499")
lines(density(diff.perc.df$tuned_GBM.t.FP.perc,n=10000),col="#AA4499",lty=2)


### False positive number
plot(density(diff.perc.df$All.f.FP,n=10000),xlim = c(-1000,15000))
lines(density(diff.perc.df$All.t.FP,n=10000),lty=3)
lines(density(diff.perc.df$GLM.f.FP,n=10000),col="#332288")
lines(density(diff.perc.df$GLM.t.FP,n=10000),col="#332288",lty=3)
lines(density(diff.perc.df$CTA.f.FP,n=10000),col="#88CCEE")
lines(density(diff.perc.df$CTA.t.FP,n=10000),col="#88CCEE",lty=3)
lines(density(diff.perc.df$MARS.f.FP,n=10000),col="#44AA99")
lines(density(diff.perc.df$MARS.t.FP,n=10000),col="#44AA99",lty=3)
lines(density(diff.perc.df$MAXNET.f.FP,n=10000),col="#999933")
lines(density(diff.perc.df$MAXNET.t.FP,n=10000),col="#999933",lty=3)
lines(density(diff.perc.df$XGBOOST.f.FP,n=10000),col="#DDCC77")
lines(density(diff.perc.df$XGBOOST.t.FP,n=10000),col="#DDCC77",lty=3)
lines(density(diff.perc.df$RF.f.FP,n=10000),col="#CC6677")
lines(density(diff.perc.df$RF.t.FP,n=10000),col="#CC6677",lty=3)
lines(density(diff.perc.df$ANN.f.FP,n=10000),col="#882255")
lines(density(diff.perc.df$ANN.t.FP,n=10000),col="#882255",lty=3)
lines(density(diff.perc.df$GBM.f.FP,n=10000),col="#AA4499")
lines(density(diff.perc.df$GBM.t.FP,n=10000),col="#AA4499",lty=3)

plot(density(diff.perc.df$All.t.FP,n=10000),xlim = c(-1000,18000),ylim = c(0,0.00028),lty=3)
lines(density(diff.perc.df$GLM.t.FP,n=10000),col="#332288",lty=3)
lines(density(diff.perc.df$CTA.t.FP,n=10000),col="#88CCEE",lty=3)
lines(density(diff.perc.df$MARS.t.FP,n=10000),col="#44AA99",lty=3)
lines(density(diff.perc.df$MAXNET.t.FP,n=10000),col="#999933",lty=3)
lines(density(diff.perc.df$XGBOOST.t.FP,n=10000),col="#DDCC77",lty=3)
lines(density(diff.perc.df$RF.t.FP,n=10000),col="#CC6677",lty=3)
lines(density(diff.perc.df$ANN.t.FP,n=10000),col="#882255",lty=3)
lines(density(diff.perc.df$GBM.t.FP,n=10000),col="#AA4499",lty=3)

plot(density(diff.perc.df$GLM.t.FP,n=10000),col="#332288",xlim = c(-1000,15000),ylim = c(0,0.00037))
lines(density(diff.perc.df$CTA.t.FP,n=10000),col="#88CCEE")
lines(density(diff.perc.df$MARS.t.FP,n=10000),col="#44AA99")
lines(density(diff.perc.df$RF.t.FP,n=10000),col="#CC6677")
lines(density(diff.perc.df$ANN.t.FP,n=10000),col="#882255")
lines(density(diff.perc.df$GBM.t.FP,n=10000),col="#AA4499")
lines(density(diff.perc.df$tuned_GLM.t.FP,n=10000),col="#332288",lty=2)
lines(density(diff.perc.df$tuned_CTA.t.FP,n=10000),col="#88CCEE",lty=2)
lines(density(diff.perc.df$tuned_MARS.t.FP,n=10000),col="#44AA99",lty=2)
lines(density(diff.perc.df$tuned_RF.t.FP,n=10000),col="#CC6677",lty=2)
lines(density(diff.perc.df$tuned_ANN.t.FP,n=10000),col="#882255",lty=2)
lines(density(diff.perc.df$tuned_GBM.t.FP,n=10000),col="#AA4499",lty=2)


### False negative percentage
plot(density(diff.perc.df$All.f.FN.perc,n=10000),xlim = c(-5,45),ylim = c(0,0.12),lty=3)
lines(density(diff.perc.df$All.t.FN.perc,n=10000),lty=3)
lines(density(diff.perc.df$GLM.f.FN.perc,n=10000),col="#332288")
lines(density(diff.perc.df$GLM.t.FN.perc,n=10000),col="#332288",lty=3)
lines(density(diff.perc.df$CTA.f.FN.perc,n=10000),col="#88CCEE")
lines(density(diff.perc.df$CTA.t.FN.perc,n=10000),col="#88CCEE",lty=3)
lines(density(diff.perc.df$MARS.f.FN.perc,n=10000),col="#44AA99")
lines(density(diff.perc.df$MARS.t.FN.perc,n=10000),col="#44AA99",lty=3)
lines(density(diff.perc.df$MAXNET.f.FN.perc,n=10000),col="#999933")
lines(density(diff.perc.df$MAXNET.t.FN.perc,n=10000),col="#999933",lty=3)
lines(density(diff.perc.df$XGBOOST.f.FN.perc,n=10000),col="#DDCC77")
lines(density(diff.perc.df$XGBOOST.t.FN.perc,n=10000),col="#DDCC77",lty=3)
lines(density(diff.perc.df$RF.f.FN.perc,n=10000),col="#CC6677")
lines(density(diff.perc.df$RF.t.FN.perc,n=10000),col="#CC6677",lty=3)
lines(density(diff.perc.df$ANN.f.FN.perc,n=10000),col="#882255")
lines(density(diff.perc.df$ANN.t.FN.perc,n=10000),col="#882255",lty=3)
lines(density(diff.perc.df$GBM.f.FN.perc,n=10000),col="#AA4499")
lines(density(diff.perc.df$GBM.t.FN.perc,n=10000),col="#AA4499",lty=3)

plot(density(diff.perc.df$All.t.FN.perc,n=10000),xlim = c(-5,45),ylim = c(0,0.11),lty=3)
lines(density(diff.perc.df$GLM.t.FN.perc,n=10000),col="#332288",lty=3)
lines(density(diff.perc.df$CTA.t.FN.perc,n=10000),col="#88CCEE",lty=3)
lines(density(diff.perc.df$MARS.t.FN.perc,n=10000),col="#44AA99",lty=3)
lines(density(diff.perc.df$MAXNET.t.FN.perc,n=10000),col="#999933",lty=3)
lines(density(diff.perc.df$XGBOOST.t.FN.perc,n=10000),col="#DDCC77",lty=3)
lines(density(diff.perc.df$RF.t.FN.perc,n=10000),col="#CC6677",lty=3)
lines(density(diff.perc.df$ANN.t.FN.perc,n=10000),col="#882255",lty=3)
lines(density(diff.perc.df$GBM.t.FN.perc,n=10000),col="#AA4499",lty=3)

plot(density(diff.perc.df$GLM.t.FN.perc,n=10000),col="#332288",xlim = c(-5,45),ylim = c(0,0.11))
lines(density(diff.perc.df$tuned_GLM.t.FN.perc,n=10000),col="#332288",lty=2)
lines(density(diff.perc.df$CTA.t.FN.perc,n=10000),col="#88CCEE")
lines(density(diff.perc.df$tuned_CTA.t.FN.perc,n=10000),col="#88CCEE",lty=2)
lines(density(diff.perc.df$MARS.t.FN.perc,n=10000),col="#44AA99")
lines(density(diff.perc.df$tuned_MARS.t.FN.perc,n=10000),col="#44AA99",lty=2)
lines(density(diff.perc.df$RF.t.FN.perc,n=10000),col="#CC6677")
lines(density(diff.perc.df$tuned_RF.t.FN.perc,n=10000),col="#CC6677",lty=2)
lines(density(diff.perc.df$ANN.t.FN.perc,n=10000),col="#882255")
lines(density(diff.perc.df$tuned_ANN.t.FN.perc,n=10000),col="#882255",lty=2)
lines(density(diff.perc.df$GBM.t.FN.perc,n=10000),col="#AA4499")
lines(density(diff.perc.df$tuned_GBM.t.FN.perc,n=10000),col="#AA4499",lty=2)


### False negative number
plot(density(diff.perc.df$All.f.FN,n=10000),xlim = c(-800,3500),ylim = c(0,0.0011))
lines(density(diff.perc.df$All.t.FN,n=10000),lty=3)
lines(density(diff.perc.df$GLM.f.FN,n=10000),col="#332288")
lines(density(diff.perc.df$GLM.t.FN,n=10000),col="#332288",lty=3)
lines(density(diff.perc.df$CTA.f.FN,n=10000),col="#88CCEE")
lines(density(diff.perc.df$CTA.t.FN,n=10000),col="#88CCEE",lty=3)
lines(density(diff.perc.df$MARS.f.FN,n=10000),col="#44AA99")
lines(density(diff.perc.df$MARS.t.FN,n=10000),col="#44AA99",lty=3)
lines(density(diff.perc.df$MAXNET.f.FN,n=10000),col="#999933")
lines(density(diff.perc.df$MAXNET.t.FN,n=10000),col="#999933",lty=3)
lines(density(diff.perc.df$XGBOOST.f.FN,n=10000),col="#DDCC77")
lines(density(diff.perc.df$XGBOOST.t.FN,n=10000),col="#DDCC77",lty=3)
lines(density(diff.perc.df$RF.f.FN,n=10000),col="#CC6677")
lines(density(diff.perc.df$RF.t.FN,n=10000),col="#CC6677",lty=3)
lines(density(diff.perc.df$ANN.f.FN,n=10000),col="#882255")
lines(density(diff.perc.df$ANN.t.FN,n=10000),col="#882255",lty=3)
lines(density(diff.perc.df$GBM.f.FN,n=10000),col="#AA4499")
lines(density(diff.perc.df$GBM.t.FN,n=10000),col="#AA4499",lty=3)

plot(density(diff.perc.df$All.t.FN,n=10000),xlim = c(-800,3500),ylim = c(0,0.0010),lty=3)
lines(density(diff.perc.df$GLM.t.FN,n=10000),col="#332288",lty=3)
lines(density(diff.perc.df$CTA.t.FN,n=10000),col="#88CCEE",lty=3)
lines(density(diff.perc.df$MARS.t.FN,n=10000),col="#44AA99",lty=3)
lines(density(diff.perc.df$MAXNET.t.FN,n=10000),col="#999933",lty=3)
lines(density(diff.perc.df$XGBOOST.t.FN,n=10000),col="#DDCC77",lty=3)
lines(density(diff.perc.df$RF.t.FN,n=10000),col="#CC6677",lty=3)
lines(density(diff.perc.df$ANN.t.FN,n=10000),col="#882255",lty=3)
lines(density(diff.perc.df$GBM.t.FN,n=10000),col="#AA4499",lty=3)

plot(density(diff.perc.df$GLM.t.FN,n=10000),col="#332288",xlim = c(-800,3500),ylim = c(0,0.0010))
lines(density(diff.perc.df$tuned_GLM.t.FN,n=10000),col="#332288",lty=2)
lines(density(diff.perc.df$CTA.t.FN,n=10000),col="#88CCEE")
lines(density(diff.perc.df$tuned_CTA.t.FN,n=10000),col="#88CCEE",lty=2)
lines(density(diff.perc.df$MARS.t.FN,n=10000),col="#44AA99")
lines(density(diff.perc.df$tuned_MARS.t.FN,n=10000),col="#44AA99",lty=2)
lines(density(diff.perc.df$RF.t.FN,n=10000),col="#CC6677")
lines(density(diff.perc.df$tuned_RF.t.FN,n=10000),col="#CC6677",lty=2)
lines(density(diff.perc.df$ANN.t.FN,n=10000),col="#882255")
lines(density(diff.perc.df$tuned_ANN.t.FN,n=10000),col="#882255",lty=2)
lines(density(diff.perc.df$GBM.t.FN,n=10000),col="#AA4499")
lines(density(diff.perc.df$tuned_GBM.t.FN,n=10000),col="#AA4499",lty=2)


### Overall inaccuracy percentage
plot(density(diff.perc.df$All.f.F.perc,n=10000),xlim = c(-5,200),ylim = c(0,0.033))
lines(density(diff.perc.df$All.t.F.perc,n=10000),lty=3)
lines(density(diff.perc.df$GLM.f.F.perc,n=10000),col="#332288")
lines(density(diff.perc.df$GLM.t.F.perc,n=10000),col="#332288",lty=3)
lines(density(diff.perc.df$CTA.f.F.perc,n=10000),col="#88CCEE")
lines(density(diff.perc.df$CTA.t.F.perc,n=10000),col="#88CCEE",lty=3)
lines(density(diff.perc.df$MARS.f.F.perc,n=10000),col="#44AA99")
lines(density(diff.perc.df$MARS.t.F.perc,n=10000),col="#44AA99",lty=3)
lines(density(diff.perc.df$MAXNET.f.F.perc,n=10000),col="#999933")
lines(density(diff.perc.df$MAXNET.t.F.perc,n=10000),col="#999933",lty=3)
lines(density(diff.perc.df$XGBOOST.f.F.perc,n=10000),col="#DDCC77")
lines(density(diff.perc.df$XGBOOST.t.F.perc,n=10000),col="#DDCC77",lty=3)
lines(density(diff.perc.df$RF.f.F.perc,n=10000),col="#CC6677")
lines(density(diff.perc.df$RF.t.F.perc,n=10000),col="#CC6677",lty=3)
lines(density(diff.perc.df$ANN.f.F.perc,n=10000),col="#882255")
lines(density(diff.perc.df$ANN.t.F.perc,n=10000),col="#882255",lty=3)
lines(density(diff.perc.df$GBM.f.F.perc,n=10000),col="#AA4499")
lines(density(diff.perc.df$GBM.t.F.perc,n=10000),col="#AA4499",lty=3)

plot(density(diff.perc.df$All.t.F.perc,n=10000),xlim = c(-5,150),ylim = c(0,0.028),lty=3)
lines(density(diff.perc.df$GLM.t.F.perc,n=10000),col="#332288",lty=3)
lines(density(diff.perc.df$CTA.t.F.perc,n=10000),col="#88CCEE",lty=3)
lines(density(diff.perc.df$MARS.t.F.perc,n=10000),col="#44AA99",lty=3)
lines(density(diff.perc.df$MAXNET.t.F.perc,n=10000),col="#999933",lty=3)
lines(density(diff.perc.df$XGBOOST.t.F.perc,n=10000),col="#DDCC77",lty=3)
lines(density(diff.perc.df$RF.t.F.perc,n=10000),col="#CC6677",lty=3)
lines(density(diff.perc.df$ANN.t.F.perc,n=10000),col="#882255",lty=3)
lines(density(diff.perc.df$GBM.t.F.perc,n=10000),col="#AA4499",lty=3)

plot(density(diff.perc.df$GLM.t.F.perc,n=10000),col="#332288",xlim = c(-0,150),ylim = c(0,0.028))
lines(density(diff.perc.df$tuned_GLM.t.F.perc,n=10000),col="#332288",lty=2)
lines(density(diff.perc.df$CTA.t.F.perc,n=10000),col="#88CCEE")
lines(density(diff.perc.df$tuned_CTA.t.F.perc,n=10000),col="#88CCEE",lty=2)
lines(density(diff.perc.df$MARS.t.F.perc,n=10000),col="#44AA99")
lines(density(diff.perc.df$tuned_MARS.t.F.perc,n=10000),col="#44AA99",lty=2)
lines(density(diff.perc.df$RF.t.F.perc,n=10000),col="#CC6677")
lines(density(diff.perc.df$tuned_RF.t.F.perc,n=10000),col="#CC6677",lty=2)
lines(density(diff.perc.df$ANN.t.F.perc,n=10000),col="#882255")
lines(density(diff.perc.df$tuned_ANN.t.F.perc,n=10000),col="#882255",lty=2)
lines(density(diff.perc.df$GBM.t.F.perc,n=10000),col="#AA4499")
lines(density(diff.perc.df$tuned_GBM.t.F.perc,n=10000),col="#AA4499",lty=2)

plot(density(diff.perc.df$GLM.f.F.perc,n=10000),col="#0077BB",xlim = c(-0,150),ylim = c(0,0.030))
lines(density(diff.perc.df$tuned_GLM.f.F.perc,n=10000),col="#88CCEE")
lines(density(diff.perc.df$CTA.f.F.perc,n=10000),col="#44AA99")
lines(density(diff.perc.df$tuned_CTA.f.F.perc,n=10000),col="#117733")
lines(density(diff.perc.df$MARS.f.F.perc,n=10000),col="#DDCC77")
lines(density(diff.perc.df$tuned_MARS.f.F.perc,n=10000),col="#999933")
lines(density(diff.perc.df$RF.f.F.perc,n=10000),col="#EE7733")
lines(density(diff.perc.df$tuned_RF.f.F.perc,n=10000),col="#CC6677")
lines(density(diff.perc.df$ANN.f.F.perc,n=10000),col="#AA4499")
lines(density(diff.perc.df$tuned_ANN.f.F.perc,n=10000),col="#882255")
lines(density(diff.perc.df$GBM.f.F.perc,n=10000),col="#BBBBBB")
lines(density(diff.perc.df$tuned_GBM.f.F.perc,n=10000),col="#332288")
lines(density(diff.perc.df$GLM.t.F.perc,n=10000),col="#0077BB",lty=2)
lines(density(diff.perc.df$tuned_GLM.t.F.perc,n=10000),col="#88CCEE",lty=2)
lines(density(diff.perc.df$CTA.t.F.perc,n=10000),col="#44AA99",lty=2)
lines(density(diff.perc.df$tuned_CTA.t.F.perc,n=10000),col="#117733",lty=2)
lines(density(diff.perc.df$MARS.t.F.perc,n=10000),col="#DDCC77",lty=2)
lines(density(diff.perc.df$tuned_MARS.t.F.perc,n=10000),col="#999933",lty=2)
lines(density(diff.perc.df$RF.t.F.perc,n=10000),col="#EE7733",lty=2)
lines(density(diff.perc.df$tuned_RF.t.F.perc,n=10000),col="#CC6677",lty=2)
lines(density(diff.perc.df$ANN.t.F.perc,n=10000),col="#AA4499",lty=2)
lines(density(diff.perc.df$tuned_ANN.t.F.perc,n=10000),col="#882255",lty=2)
lines(density(diff.perc.df$GBM.t.F.perc,n=10000),col="#BBBBBB",lty=2)
lines(density(diff.perc.df$tuned_GBM.t.F.perc,n=10000),col="#332288",lty=2)


plot(density(diff.perc.df$GLM.t.F.perc,n=10000),col="#0077BB",lty=2,xlim = c(-0,150),ylim = c(0,0.030))
lines(density(diff.perc.df$tuned_GLM.t.F.perc,n=10000),col="#88CCEE",lty=2)
lines(density(diff.perc.df$CTA.t.F.perc,n=10000),col="#44AA99",lty=2)
lines(density(diff.perc.df$tuned_CTA.t.F.perc,n=10000),col="#117733",lty=2)
lines(density(diff.perc.df$MARS.t.F.perc,n=10000),col="#DDCC77",lty=2)
lines(density(diff.perc.df$tuned_MARS.t.F.perc,n=10000),col="#999933",lty=2)
lines(density(diff.perc.df$RF.t.F.perc,n=10000),col="#EE7733",lty=2)
lines(density(diff.perc.df$tuned_RF.t.F.perc,n=10000),col="#CC6677",lty=2)
lines(density(diff.perc.df$ANN.t.F.perc,n=10000),col="#AA4499",lty=2)
lines(density(diff.perc.df$tuned_ANN.t.F.perc,n=10000),col="#882255",lty=2)
lines(density(diff.perc.df$GBM.t.F.perc,n=10000),col="#BBBBBB",lty=2)
lines(density(diff.perc.df$tuned_GBM.t.F.perc,n=10000),col="#332288",lty=2)


plot(density(diff.perc.df$tuned_GLM.f.F.perc,n=10000),col="#88CCEE",xlim = c(-0,150),ylim = c(0,0.030))
lines(density(diff.perc.df$tuned_CTA.f.F.perc,n=10000),col="#117733")
lines(density(diff.perc.df$tuned_MARS.f.F.perc,n=10000),col="#999933")
lines(density(diff.perc.df$tuned_RF.f.F.perc,n=10000),col="#CC6677")
lines(density(diff.perc.df$tuned_ANN.f.F.perc,n=10000),col="#882255")
lines(density(diff.perc.df$tuned_GBM.f.F.perc,n=10000),col="#332288")
lines(density(diff.perc.df$tuned_GLM.t.F.perc,n=10000),col="#88CCEE",lty=2)
lines(density(diff.perc.df$tuned_CTA.t.F.perc,n=10000),col="#117733",lty=2)
lines(density(diff.perc.df$tuned_MARS.t.F.perc,n=10000),col="#999933",lty=2)
lines(density(diff.perc.df$tuned_RF.t.F.perc,n=10000),col="#CC6677",lty=2)
lines(density(diff.perc.df$tuned_ANN.t.F.perc,n=10000),col="#882255",lty=2)
lines(density(diff.perc.df$tuned_GBM.t.F.perc,n=10000),col="#332288",lty=2)


### Overall inaccuracy number
plot(density(diff.perc.df$All.f.F,n=10000),xlim = c(-1500,22000),ylim = c(0,0.0002))
lines(density(diff.perc.df$All.t.F,n=10000),lty=3)
lines(density(diff.perc.df$GLM.f.F,n=10000),col="#332288")
lines(density(diff.perc.df$GLM.t.F,n=10000),col="#332288",lty=3)
lines(density(diff.perc.df$CTA.f.F,n=10000),col="#88CCEE")
lines(density(diff.perc.df$CTA.t.F,n=10000),col="#88CCEE",lty=3)
lines(density(diff.perc.df$MARS.f.F,n=10000),col="#44AA99")
lines(density(diff.perc.df$MARS.t.F,n=10000),col="#44AA99",lty=3)
lines(density(diff.perc.df$MAXNET.f.F,n=10000),col="#999933")
lines(density(diff.perc.df$MAXNET.t.F,n=10000),col="#999933",lty=3)
lines(density(diff.perc.df$XGBOOST.f.F,n=10000),col="#DDCC77")
lines(density(diff.perc.df$XGBOOST.t.F,n=10000),col="#DDCC77",lty=3)
lines(density(diff.perc.df$RF.f.F,n=10000),col="#CC6677")
lines(density(diff.perc.df$RF.t.F,n=10000),col="#CC6677",lty=3)
lines(density(diff.perc.df$ANN.f.F,n=10000),col="#882255")
lines(density(diff.perc.df$ANN.t.F,n=10000),col="#882255",lty=3)
lines(density(diff.perc.df$GBM.f.F,n=10000),col="#AA4499")
lines(density(diff.perc.df$GBM.t.F,n=10000),col="#AA4499",lty=3)

plot(density(diff.perc.df$All.t.F,n=10000),xlim = c(-1500,22000),ylim = c(0,0.00016),lty=3)
lines(density(diff.perc.df$GLM.t.F,n=10000),col="#332288",lty=3)
lines(density(diff.perc.df$CTA.t.F,n=10000),col="#88CCEE",lty=3)
lines(density(diff.perc.df$MARS.t.F,n=10000),col="#44AA99",lty=3)
lines(density(diff.perc.df$MAXNET.t.F,n=10000),col="#999933",lty=3)
lines(density(diff.perc.df$XGBOOST.t.F,n=10000),col="#DDCC77",lty=3)
lines(density(diff.perc.df$RF.t.F,n=10000),col="#CC6677",lty=3)
lines(density(diff.perc.df$ANN.t.F,n=10000),col="#882255",lty=3)
lines(density(diff.perc.df$GBM.t.F,n=10000),col="#AA4499",lty=3)

plot(density(diff.perc.df$GLM.t.F,n=10000),col="#332288",xlim = c(-1500,22000),ylim = c(0,0.00020))
lines(density(diff.perc.df$tuned_GLM.t.F,n=10000),col="#332288",lty=2)
lines(density(diff.perc.df$CTA.t.F,n=10000),col="#88CCEE")
lines(density(diff.perc.df$tuned_CTA.t.F,n=10000),col="#88CCEE",lty=2)
lines(density(diff.perc.df$MARS.t.F,n=10000),col="#44AA99")
lines(density(diff.perc.df$tuned_MARS.t.F,n=10000),col="#44AA99",lty=2)
lines(density(diff.perc.df$RF.t.F,n=10000),col="#CC6677")
lines(density(diff.perc.df$tuned_RF.t.F,n=10000),col="#CC6677",lty=2)
lines(density(diff.perc.df$ANN.t.F,n=10000),col="#882255")
lines(density(diff.perc.df$tuned_ANN.t.F,n=10000),col="#882255",lty=2)
lines(density(diff.perc.df$GBM.t.F,n=10000),col="#AA4499")
lines(density(diff.perc.df$tuned_GBM.t.F,n=10000),col="#AA4499",lty=2)
lines(density(diff.perc.df$XGBOOST.t.F,n=10000),col="#999933")
lines(density(diff.perc.df$MAXNET.t.F,n=10000),col="#DDCC77")



########## Stastistically test results ##########

ks.test(diff.perc.df$GBM.t.F,
        diff.perc.df$tuned_GBM.t.F,
        alternative = "two.sided",exact = NULL, 
        simulate.p.value=FALSE)

ks.test(diff.perc.df$GLM.f.F, #x
        diff.perc.df$tuned_GLM.f.F, #y
        alternative = "greater",exact = NULL, #Alternative: X>Y
        simulate.p.value=FALSE)

ks.test(diff.perc.df$MARS.f.F,
        diff.perc.df$tuned_MARS.f.F,
        alternative = "two.sided",exact = NULL, 
        simulate.p.value=FALSE)

ks.test(diff.perc.df$RF.t.F, #x
        diff.perc.df$MARS.t.F, #y
        alternative = "greater",exact = NULL, #Alternative: X>Y
        simulate.p.value=FALSE)

ks.test(diff.perc.df$MAXNET.t.F, #x
        diff.perc.df$tuned_GBM.t.F, #y
        alternative = "greater",exact = NULL, #Alternative: X>Y
        simulate.p.value=FALSE)








###### Model tuning results ######


plot(density(diff.perc.df$GLM.f.F,n=10000),col="#332288",xlim=c(-5000,40000),ylim=c(0,0.0003))
lines(density(diff.perc.df$tuned_GLM.f.F,n=10000),col="#88CCEE")
lines(density(diff.perc.df$GLM.t.F,n=10000),col="#332288",lty=2)
lines(density(diff.perc.df$tuned_GLM.t.F,n=10000),col="#88CCEE",lty=2)

plot(density(diff.perc.df$GBM.f.F,n=10000),col="#332288",xlim=c(-5000,40000),ylim=c(0,0.0003))
lines(density(diff.perc.df$tuned_GBM.f.F,n=10000),col="#88CCEE")
lines(density(diff.perc.df$GBM.t.F,n=10000),col="#332288",lty=2)
lines(density(diff.perc.df$tuned_GBM.t.F,n=10000),col="#88CCEE",lty=2)

plot(density(diff.perc.df$CTA.f.F,n=10000),col="#332288",xlim=c(-5000,40000),ylim=c(0,0.0003))
lines(density(diff.perc.df$tuned_CTA.f.F,n=10000),col="#88CCEE")
lines(density(diff.perc.df$CTA.t.F,n=10000),col="#332288",lty=2)
lines(density(diff.perc.df$tuned_CTA.t.F,n=10000),col="#88CCEE",lty=2)

plot(density(diff.perc.df$MARS.f.F,n=10000),col="#332288",xlim=c(-5000,40000),ylim=c(0,0.0003))
lines(density(diff.perc.df$tuned_MARS.f.F,n=10000),col="#88CCEE")
lines(density(diff.perc.df$MARS.t.F,n=10000),col="#332288",lty=2)
lines(density(diff.perc.df$tuned_MARS.t.F,n=10000),col="#88CCEE",lty=2)

plot(density(diff.perc.df$RF.f.F,n=10000),col="#332288",xlim=c(-5000,40000),ylim=c(0,0.0003))
lines(density(diff.perc.df$tuned_RF.f.F,n=10000),col="#88CCEE")
lines(density(diff.perc.df$RF.t.F,n=10000),col="#332288",lty=2)
lines(density(diff.perc.df$tuned_RF.t.F,n=10000),col="#88CCEE",lty=2)

plot(density(diff.perc.df$ANN.f.F,n=10000),col="#332288",xlim=c(-5000,40000),ylim=c(0,0.0003))
lines(density(diff.perc.df$tuned_ANN.f.F,n=10000),col="#88CCEE")
lines(density(diff.perc.df$ANN.t.F,n=10000),col="#332288",lty=2)
lines(density(diff.perc.df$tuned_ANN.t.F,n=10000),col="#88CCEE",lty=2)







plot(density(diff.perc.df$GLM.f.F.perc,n=10000),col="#332288",xlim=c(-10,170),ylim=c(0,0.036))
lines(density(diff.perc.df$tuned_GLM.f.F.perc,n=10000),col="#88CCEE")
lines(density(diff.perc.df$GLM.t.F.perc,n=10000),col="#332288",lty=2)
lines(density(diff.perc.df$tuned_GLM.t.F.perc,n=10000),col="#88CCEE",lty=2)

plot(density(diff.perc.df$GBM.f.F.perc,n=10000),col="#332288",xlim=c(-10,170),ylim=c(0,0.036))
lines(density(diff.perc.df$tuned_GBM.f.F.perc,n=10000),col="#88CCEE")
lines(density(diff.perc.df$GBM.t.F.perc,n=10000),col="#332288",lty=2)
lines(density(diff.perc.df$tuned_GBM.t.F.perc,n=10000),col="#88CCEE",lty=2)

plot(density(diff.perc.df$CTA.f.F.perc,n=10000),col="#332288",xlim=c(-10,170),ylim=c(0,0.036))
lines(density(diff.perc.df$tuned_CTA.f.F.perc,n=10000),col="#88CCEE")
lines(density(diff.perc.df$CTA.t.F.perc,n=10000),col="#332288",lty=2)
lines(density(diff.perc.df$tuned_CTA.t.F.perc,n=10000),col="#88CCEE",lty=2)

plot(density(diff.perc.df$MARS.f.F.perc,n=10000),col="#332288",xlim=c(-10,170),ylim=c(0,0.036))
lines(density(diff.perc.df$tuned_MARS.f.F.perc,n=10000),col="#88CCEE")
lines(density(diff.perc.df$MARS.t.F.perc,n=10000),col="#332288",lty=2)
lines(density(diff.perc.df$tuned_MARS.t.F.perc,n=10000),col="#88CCEE",lty=2)

plot(density(diff.perc.df$RF.f.F.perc,n=10000),col="#332288",xlim=c(-10,170),ylim=c(0,0.036))
lines(density(diff.perc.df$tuned_RF.f.F.perc,n=10000),col="#88CCEE")
lines(density(diff.perc.df$RF.t.F.perc,n=10000),col="#332288",lty=2)
lines(density(diff.perc.df$tuned_RF.t.F.perc,n=10000),col="#88CCEE",lty=2)

plot(density(diff.perc.df$ANN.f.F.perc,n=10000),col="#332288",xlim=c(-10,170),ylim=c(0,0.036))
lines(density(diff.perc.df$tuned_ANN.f.F.perc,n=10000),col="#88CCEE")
lines(density(diff.perc.df$ANN.t.F.perc,n=10000),col="#332288",lty=2)
lines(density(diff.perc.df$tuned_ANN.t.F.perc,n=10000),col="#88CCEE",lty=2)




