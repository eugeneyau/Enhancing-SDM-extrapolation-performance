
library(terra)
options(scipen = 999)


Analysis <- "EM" # "TSS" or "EM"

filename.base <- "biomod424_BASE_REV" #Base distribution
filename <- "bmod424_trunc_rev_bestalgo_SOUTH_trun20perc"

save.filename <- paste0(Analysis,"_",filename) # Name the df output file

VirtualEM <- "EMmean"
EM1 <- "EMca"
EM2 <- "EMmean"
EM3 <- "EMwmean"
EM4 <- "EMmedian"

TSSprojection1 <- "proj_TSSNA"
TSSprojection2 <- "proj_TSSNA"
TSSprojection3 <- "proj_TSSNA"
TSSprojection4 <- "proj_TSSNA"

projection1 <- paste0("/",TSSprojection1,"/",TSSprojection1,"_")
projection2 <- paste0("/",TSSprojection2,"/",TSSprojection2,"_")
projection3 <- paste0("/",TSSprojection3,"/",TSSprojection3,"_")
projection4 <- paste0("/",TSSprojection4,"/",TSSprojection4,"_")

BASEprojection <- "/proj_Current/proj_Current_"

species.list.out <- list.files(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename))
basepath <- paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename,"/")
#unlist(species.list.out1)
poss1 <- paste0(basepath,species.list.out,projection1,species.list.out,"_ensemble_TSSbin.tif")
tf1 <- file.exists(poss1)
splist1 <- species.list.out[tf1]

poss2 <- paste0(basepath,species.list.out,projection2,species.list.out,"_ensemble_TSSbin.tif")
tf2 <- file.exists(poss2)
splist2 <- species.list.out[tf2]

poss3 <- paste0(basepath,species.list.out,projection3,species.list.out,"_ensemble_TSSbin.tif")
tf3 <- file.exists(poss3)
splist3 <- species.list.out[tf3]

poss4 <- paste0(basepath,species.list.out,projection4,species.list.out,"_ensemble_TSSbin.tif")
tf4 <- file.exists(poss4)
splist4 <- species.list.out[tf4]


########## Common species list ##########
common.list1 <- intersect(splist1, splist2) 
common.list2 <- intersect(splist3, splist4) 
common.list <- intersect(common.list1, common.list2) 
#unlist(common.list)


### Clean memory
rm(poss1, poss2, poss3, poss4, tf1, tf2, tf3, tf4)
rm(splist1, splist2, splist3, splist4, common.list1, common.list2)
gc(full = TRUE)


########## FPOS FNEG analysis ##########

### Create df for documentation
if (Analysis == "EM") {
  columns <- c("Species",paste0(EM1,".FP"),paste0(EM1,".FP.perc"),paste0(EM1,".FN"),
               paste0(EM1,".FN.perc"),paste0(EM1,".F"),paste0(EM1,".F.perc"),
               paste0(EM2,".FP"),paste0(EM2,".FP.perc"),paste0(EM2,".FN"),
               paste0(EM2,".FN.perc"),paste0(EM2,".F"),paste0(EM2,".F.perc"),
               paste0(EM3,".FP"),paste0(EM3,".FP.perc"),paste0(EM3,".FN"),
               paste0(EM3,".FN.perc"),paste0(EM3,".F"),paste0(EM3,".F.perc"),
               paste0(EM4,".FP"),paste0(EM4,".FP.perc"),paste0(EM4,".FN"),
               paste0(EM4,".FN.perc"),paste0(EM4,".F"),paste0(EM4,".F.perc"))
}

diff.perc.df <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(diff.perc.df) <- columns

### Do the analysis

# Automated calculation of error for each species
Accuracy_Fun <- function(i) {
  
  # Load maps
  map1 <- terra::rast(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename,"/",common.list[[i]],projection1,common.list[[i]],"_ensemble_TSSbin.tif"))
  map1 <- map1[EM1]
  map2 <- terra::rast(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename,"/",common.list[[i]],projection2,common.list[[i]],"_ensemble_TSSbin.tif"))
  map2 <- map2[EM2]
  map3 <- terra::rast(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename,"/",common.list[[i]],projection3,common.list[[i]],"_ensemble_TSSbin.tif"))
  map3 <- map3[EM3]
  map4 <- terra::rast(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename,"/",common.list[[i]],projection4,common.list[[i]],"_ensemble_TSSbin.tif"))
  map4 <- map4[EM4]
  
  truerange <- terra::rast(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDM_BASE/results_",filename.base,"/",common.list[[i]],BASEprojection,common.list[[i]],"_ensemble_TSSbin.tif"))
  truerange <- truerange[VirtualEM]
  truerange <- terra::crop(truerange,terra::ext(map1))
  
  
  # Get range
  truerange_size <- subset(freq(truerange), value==1)[,3]
  
  # False Positives
  FPOS1 <- sum(values(truerange<map1), na.rm=TRUE)
  FPOS.perc1 <- FPOS1/truerange_size*100
  
  FPOS2 <- sum(values(truerange<map2), na.rm=TRUE)
  FPOS.perc2 <- FPOS2/truerange_size*100
  
  FPOS3 <- sum(values(truerange<map3), na.rm=TRUE)
  FPOS.perc3 <- FPOS3/truerange_size*100
  
  FPOS4 <- sum(values(truerange<map4), na.rm=TRUE)
  FPOS.perc4 <- FPOS4/truerange_size*100
  
  # False Negatives
  FNEG1 <- sum(values(truerange>map1), na.rm=TRUE)
  FNEG.perc1 <- FNEG1/truerange_size*100
  
  FNEG2 <- sum(values(truerange>map2), na.rm=TRUE)
  FNEG.perc2 <- FNEG2/truerange_size*100
  
  FNEG3 <- sum(values(truerange>map3), na.rm=TRUE)
  FNEG.perc3 <- FNEG3/truerange_size*100
  
  FNEG4 <- sum(values(truerange>map4), na.rm=TRUE)
  FNEG.perc4 <- FNEG4/truerange_size*100
  
  # Total wrong predictions
  FSUM1 <- FPOS1+FNEG1
  FSUM.perc1 <- FPOS.perc1+FNEG.perc1
  
  FSUM2 <- FPOS2+FNEG2
  FSUM.perc2 <- FPOS.perc2+FNEG.perc2
  
  FSUM3 <- FPOS3+FNEG3
  FSUM.perc3 <- FPOS.perc3+FNEG.perc3
  
  FSUM4 <- FPOS4+FNEG4
  FSUM.perc4 <- FPOS.perc4+FNEG.perc4
  
  # Add to df
  sp.data <- list(common.list[[i]],FPOS1,FPOS.perc1,FNEG1,FNEG.perc1,FSUM1,FSUM.perc1,
                  FPOS2,FPOS.perc2,FNEG2,FNEG.perc2,FSUM2,FSUM.perc2,
                  FPOS3,FPOS.perc3,FNEG3,FNEG.perc3,FSUM3,FSUM.perc3,
                  FPOS4,FPOS.perc4,FNEG4,FNEG.perc4,FSUM4,FSUM.perc4)
  
  diff.perc.df[nrow(diff.perc.df)+ 1,] <<- sp.data
}

i=1:length(common.list)

lapply(i,Accuracy_Fun)

### Save output df
saveRDS(diff.perc.df,paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",save.filename,".RData"))

### Load output df
diff.perc.df <- readRDS(paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",save.filename,".RData"))






########## Stastistically test results ##########
wd <- "/lustre1/g/sbs_bonebrake/Eugene/Tempdata/"
setwd(wd)

stat.filename.allalgo <- "EM_bmod424_trunc_rev_SOUTH_trun20perc" # All untuned algo
stat.filename.selectedalgo <- "EM_bmod424_trunc_rev_bestalgo_SOUTH_trun20perc" # Selected untuned algo

allalgo.df <- readRDS(paste0(stat.filename.allalgo,".RData"))
selected.df <- readRDS(paste0(stat.filename.selectedalgo,".RData"))

ks.test(allalgo.df$EMca.F,
        selected.df$EMca.F,
        alternative = "two.sided",exact = NULL, 
        simulate.p.value=FALSE)

ks.test(allalgo.df$EMmean.F,
        selected.df$EMmean.F,
        alternative = "two.sided",exact = NULL, 
        simulate.p.value=FALSE)

ks.test(allalgo.df$EMmedian.F,
        selected.df$EMmedian.F,
        alternative = "two.sided",exact = NULL, 
        simulate.p.value=FALSE)

ks.test(selected.df$EMca.F,
        allalgo.df$EMca.F,
        alternative = "greater",exact = NULL,
        simulate.p.value=FALSE)

ks.test(selected.df$EMmean.F,
        allalgo.df$EMmean.F,
        alternative = "greater",exact = NULL,
        simulate.p.value=FALSE)

ks.test(selected.df$EMmedian.F,
        allalgo.df$EMmedian.F,
        alternative = "greater",exact = NULL,
        simulate.p.value=FALSE)
