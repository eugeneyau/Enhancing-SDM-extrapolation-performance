
library(terra)
options(scipen = 999)

# Define the list of algorithms
algo <- sort(c("ANN","CTA","GLM","MARS","MAXNET","RF","XGBOOST","GBM","tuned_ANN","tuned_CTA","tuned_GLM","tuned_GBM","tuned_MARS","tuned_RF"))
#algo <- sort(c("ANN","CTA","GLM","MARS","MAXNET","RF","XGBOOST","GBM","tuned_ANN","tuned_CTA","tuned_GLM","tuned_GBM","tuned_0interact_MARS","tuned_RF"))

# Specify filenames
Virtualsp_dir <- "biomod424_BASE_REV" #Virtual distribution
All_algo_full_dir <- "bmod424_trunc_rev_SOUTH_trun0perc" #All algo FULL
All_algo_trun_dir <- "bmod424_trunc_rev_SOUTH_trun20perc" #All algo TRUN
Best_algo_full_dir <- "NA" #Best algo FULL
Best_algo_trun_dir <- "bmod424_trunc_rev_bestalgo2_SOUTH_trun20perc" #Best algo TRUN

# List filenames
filenames <- c(
  Virtualsp_dir,
  All_algo_full_dir,
  All_algo_trun_dir,
  sort(c(paste0("bmod424_trunc_rev_", algo, "_SOUTH_trun0perc"),
         paste0("bmod424_trunc_rev_", algo, "_SOUTH_trun20perc")))
)

if (!Best_algo_full_dir=="NA") {filenames <- c(filenames, Best_algo_full_dir)}
if (!Best_algo_trun_dir=="NA") {filenames <- c(filenames, Best_algo_trun_dir)}

# Settings
VirtualEM <- "EMmean"
EM <- "EMca"
AlgoEM <- "EMca"
TSSprojection <- "proj_TSSNA"

save.filename <- paste0("Algo_REV_",EM,"_",TSSprojection) # Name the df output file

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
  
  truerange_size <- subset(freq(truerange), value == 1)[, 3]
  
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
    Algo_full_dir <- paste0(base_dir,"bmod424_trunc_rev_",algo[k],"_SOUTH_trun0perc")
    Algo_trun_dir <- paste0(base_dir,"bmod424_trunc_rev_",algo[k],"_SOUTH_trun20perc")
    
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
