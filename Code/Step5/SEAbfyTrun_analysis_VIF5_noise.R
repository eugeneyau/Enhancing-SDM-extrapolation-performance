
library(terra)
options(scipen = 999)


filename.base <- "biomod424_BASE_VIF5" #Virtual distribution
runID <- "VIF5_SOUTH_trun20" #Standard
#filename.base <- "biomod424_BASE_PC_V2" #Virtual distribution
#runID <- "dataV1.1_PCvsp_SOUTH_trun20" #PC

filename1 <- paste0("bmod424_trunc_",runID,"perc") #All algo trun
filename2 <- paste0("bmod424_trunc_bio_phchndvi_",runID,"perc") #Bioclim + pH + ch trun
filename3 <- paste0("bmod424_trunc_varnoise_",runID,"perc") #All + Noise trun

VirtualEM <- "EMmean"
EM <- "EMca"
TSSprojection <- "proj_TSSNA"

save.filename <- paste0("Biophchndvi_",EM,"_",TSSprojection,"_",runID) # Name the df output file

projection <- paste0("/",TSSprojection,"/",TSSprojection,"_")
BASEprojection <- "/proj_Current/proj_Current_"

basepath1 <- paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename1)
species.list.out <- list.files(basepath1)
poss1 <- paste0(basepath1,"/",species.list.out,projection,species.list.out,"_ensemble_TSSbin.tif")
tf1 <- file.exists(poss1)
splist1 <- species.list.out[tf1]

basepath2 <- paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename2)
species.list.out <- list.files(basepath2)
poss2 <- paste0(basepath2,"/",species.list.out,projection,species.list.out,"_ensemble_TSSbin.tif")
tf2 <- file.exists(poss2)
splist2 <- species.list.out[tf2]

basepath3 <- paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename3)
species.list.out <- list.files(basepath3)
poss3 <- paste0(basepath3,"/",species.list.out,projection,species.list.out,"_ensemble_TSSbin.tif")
tf3 <- file.exists(poss3)
splist3 <- species.list.out[tf3]


########## Common species list ##########
common.list1 <- intersect(splist1, splist2) 
common.list <- intersect(common.list1, splist3) 
#unlist(common.list)

### Clean memory
rm(poss1, poss2, poss3, tf1, tf2, tf3)
rm(splist1, splist2, splist3, common.list1)
gc(full = TRUE)


########## FPOS FNEG analysis ##########

### Create df for documentation
columns <- c("Species","All.t.FP","All.t.FP.perc","All.t.FN","All.t.FN.perc","All.t.F","All.t.F.perc",
             "BIO.t.FP","BIO.t.FP.perc","BIO.t.FN","BIO.t.FN.perc","BIO.t.F","BIO.t.F.perc",
             "NOI.t.FP","NOI.t.FP.perc","NOI.t.FN","NOI.t.FN.perc","NOI.t.F","NOI.t.F.perc")
diff.perc.df <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(diff.perc.df) <- columns

### Do the analysis

# Automated calculation of error for each species
Accuracy_Fun <- function(i) {
  
  # Load maps
  map1 <- terra::rast(paste0(basepath1,"/",common.list[[i]],projection,common.list[[i]],"_ensemble_TSSbin.tif"))
  map1 <- map1[EM]
  map2 <- terra::rast(paste0(basepath2,"/",common.list[[i]],projection,common.list[[i]],"_ensemble_TSSbin.tif"))
  map2 <- map2[EM]
  map3 <- terra::rast(paste0(basepath3,"/",common.list[[i]],projection,common.list[[i]],"_ensemble_TSSbin.tif"))
  map3 <- map3[EM]
  
  truerange <- terra::rast(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDM_BASE/results_",filename.base,"/",common.list[[i]],BASEprojection,common.list[[i]],"_ensemble_TSSbin.tif"))
  truerange <- truerange[VirtualEM]
  truerange <- terra::crop(truerange,terra::ext(map1))
  
  # Get range
  truerange_size <- subset(terra::freq(truerange), value==1)[,3]
  
  # False Positives
  FPOS1 <- sum(values(truerange<map1), na.rm=TRUE)
  FPOS.perc1 <- FPOS1/truerange_size*100
  
  FPOS2 <- sum(values(truerange<map2), na.rm=TRUE)
  FPOS.perc2 <- FPOS2/truerange_size*100
  
  FPOS3 <- sum(values(truerange<map3), na.rm=TRUE)
  FPOS.perc3 <- FPOS3/truerange_size*100
  
  
  # False Negatives
  FNEG1 <- sum(values(truerange>map1), na.rm=TRUE)
  FNEG.perc1 <- FNEG1/truerange_size*100
  
  FNEG2 <- sum(values(truerange>map2), na.rm=TRUE)
  FNEG.perc2 <- FNEG2/truerange_size*100
  
  FNEG3 <- sum(values(truerange>map3), na.rm=TRUE)
  FNEG.perc3 <- FNEG3/truerange_size*100
  
  
  # Total wrong predictions
  FSUM1 <- FPOS1+FNEG1
  FSUM.perc1 <- FPOS.perc1+FNEG.perc1
  
  FSUM2 <- FPOS2+FNEG2
  FSUM.perc2 <- FPOS.perc2+FNEG.perc2
  
  FSUM3 <- FPOS3+FNEG3
  FSUM.perc3 <- FPOS.perc3+FNEG.perc3
  
  # Add to df
  sp.data <- list(common.list[[i]],FPOS1,FPOS.perc1,FNEG1,FNEG.perc1,FSUM1,FSUM.perc1,
                  FPOS2,FPOS.perc2,FNEG2,FNEG.perc2,FSUM2,FSUM.perc2,
                  FPOS3,FPOS.perc3,FNEG3,FNEG.perc3,FSUM3,FSUM.perc3)
  
  diff.perc.df[nrow(diff.perc.df)+ 1,] <<- sp.data
}

i=1:length(common.list)

lapply(i,Accuracy_Fun)

### Save output df
saveRDS(diff.perc.df,paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",save.filename,".RData"))


########## Visualise results ##########

### Load output df
diff.perc.df <- readRDS(paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",save.filename,".RData"))

### False positive percentage
plot(density(diff.perc.df$All.t.FP.perc,n=10000),xlim = c(-20,170),ylim = c(0,0.03))
lines(density(diff.perc.df$BIO.t.FP.perc,n=10000),col="red")
lines(density(diff.perc.df$NOI.t.FP.perc,n=10000),col="blue")

### False positive number
plot(density(diff.perc.df$All.t.FP,n=10000),xlim = c(-2000,20000),ylim = c(0,0.00038))
lines(density(diff.perc.df$BIO.t.FP,n=10000),col="red")
lines(density(diff.perc.df$NOI.t.FP,n=10000),col="blue")


### False negative percentage
plot(density(diff.perc.df$All.t.FN.perc,n=10000),xlim = c(-5,60),ylim = c(0,0.036))
lines(density(diff.perc.df$BIO.t.FN.perc,n=10000),col="red")
lines(density(diff.perc.df$NOI.t.FN.perc,n=10000),col="blue")

### False negative number
plot(density(diff.perc.df$All.t.FN,n=10000),xlim = c(-1500,9000),ylim = c(0,0.00025))
lines(density(diff.perc.df$BIO.t.FN,n=10000),col="red")
lines(density(diff.perc.df$NOI.t.FN,n=10000),col="blue")


### Overall inaccuracy percentage
plot(density(diff.perc.df$All.t.F.perc,n=10000),xlim = c(0,100),ylim = c(0,0.035))
lines(density(diff.perc.df$BIO.t.F.perc,n=10000),col="red")
lines(density(diff.perc.df$NOI.t.F.perc,n=10000),col="blue")

### Overall inaccuracy number
plot(density(diff.perc.df$All.t.F,n=10000),xlim = c(-1000,16000),ylim = c(0,0.00015))
lines(density(diff.perc.df$BIO.t.F,n=10000),col="red")
lines(density(diff.perc.df$NOI.t.F,n=10000),col="blue")


########## Stastistically test results ##########

ks.test(diff.perc.df$All.t.F, #x
        diff.perc.df$BIO.t.F, #y
        alternative = "greater",exact = NULL, #Alternative == X>Y
        simulate.p.value=FALSE)

ks.test(diff.perc.df$NOI.t.F, 
        diff.perc.df$All.t.F, 
        alternative = "two.sided",exact = NULL, 
        simulate.p.value=FALSE)



