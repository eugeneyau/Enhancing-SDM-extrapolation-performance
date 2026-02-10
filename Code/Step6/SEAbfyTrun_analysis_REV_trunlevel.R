

library(terra)
options(scipen = 999)

ncores <- 8

filename.base <- "biomod424_BASE_REV" #Base distribution
filename1 <- "bmod424_RAWvspEMmean_trunc_rev_lat_0-36" #Full
filename2 <- "bmod424_RAWvspEMmean_trunc_rev_lat_5-36" #5-36
filename3 <- "bmod424_RAWvspEMmean_trunc_rev_lat_10-36" #10-36
filename4 <- "bmod424_RAWvspEMmean_trunc_rev_lat_15-36" #15-36

VirtualEM <- "EMmean"
EM <- "EMca"

save.filename <- paste0("Trunlevel_REV_",EM) # Name the df output file
TSSprojection <- "proj_TSSNA"
projection <- paste0("/",TSSprojection,"/",TSSprojection,"_")

BASEprojection <- "/proj_Current/proj_Current_"

species.list.out1 <- list.files(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename1))
basepath1 <- paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename1,"/")
#unlist(species.list.out1)
poss1 <- paste0(basepath1,species.list.out1,projection,species.list.out1,"_ensemble_TSSbin.tif")
tf1 <- file.exists(poss1)
splist1 <- species.list.out1[tf1]

species.list.out2 <- list.files(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename2))
basepath2 <- paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename2,"/")
#unlist(species.list.out2)
poss2 <- paste0(basepath2,species.list.out2,projection,species.list.out2,"_ensemble_TSSbin.tif")
tf2 <- file.exists(poss2)
splist2 <- species.list.out2[tf2]

species.list.out3 <- list.files(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename3))
basepath3 <- paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename3,"/")
#unlist(species.list.out3)
poss3 <- paste0(basepath3,species.list.out3,projection,species.list.out3,"_ensemble_TSSbin.tif")
tf3 <- file.exists(poss3)
splist3 <- species.list.out3[tf3]

species.list.out4 <- list.files(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename4))
basepath4 <- paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename4,"/")
#unlist(species.list.out4)
poss4 <- paste0(basepath4,species.list.out4,projection,species.list.out4,"_ensemble_TSSbin.tif")
tf4 <- file.exists(poss4)
splist4 <- species.list.out4[tf4]


########## Common species list ##########
common.list1 <- intersect(splist1, splist2) 
common.list2 <- intersect(splist3, splist4) 
common.list <- intersect(common.list1, common.list2) 
#unlist(common.list)

### Clean memory
rm(species.list.out1, species.list.out2, species.list.out3, species.list.out4)
rm(poss1, poss2, poss3, poss4, tf1, tf2, tf3, tf4)
rm(splist1, splist2, splist3, splist4, common.list1, common.list2)
gc(full = TRUE)


########## Projected alpha diversity maps ##########

### Diversity pattern
species.list <- list.files(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDM_BASE/results_",filename.base))
my_stack_base_allvsp <- terra::rast(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDM_BASE/results_",filename.base,"/",species.list,BASEprojection,species.list,"_ensemble_TSSbin.tif"))
my_stack_base_allvsp <- my_stack_base_allvsp[VirtualEM]
alpha_base_allvsp <- terra::app(my_stack_base_allvsp, fun=sum, cores=ncores)
plot(alpha_base_allvsp)

my_stack_base <- terra::rast(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDM_BASE/results_",filename.base,"/",common.list,BASEprojection,common.list,"_ensemble_TSSbin.tif"))
my_stack_base <- my_stack_base[VirtualEM]
alpha_base <- terra::app(my_stack_base, fun=sum, cores=ncores)
plot(alpha_base)

my_stack1 <- terra::rast(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename1,"/",common.list,projection,common.list,"_ensemble_TSSbin.tif"))
my_stack1 <- my_stack1[EM]
alpha1 <- terra::app(my_stack1, fun=sum, cores=ncores)
plot(alpha1)

my_stack2 <- terra::rast(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename2,"/",common.list,projection,common.list,"_ensemble_TSSbin.tif"))
my_stack2 <- my_stack2[EM]
alpha2 <- terra::app(my_stack2, fun=sum, cores=ncores)
plot(alpha2)

my_stack3 <- terra::rast(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename3,"/",common.list,projection,common.list,"_ensemble_TSSbin.tif"))
my_stack3 <- my_stack3[EM]
alpha3 <- terra::app(my_stack3, fun=sum, cores=ncores)
plot(alpha3)

my_stack4 <- terra::rast(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename4,"/",common.list,projection,common.list,"_ensemble_TSSbin.tif"))
my_stack4 <- my_stack4[EM]
alpha4 <- terra::app(my_stack4, fun=sum, cores=ncores)
plot(alpha4)

writeRaster(alpha_base_allvsp, paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",filename.base,"_diversity_allvsp.tif"), overwrite=FALSE)
writeRaster(alpha_base, paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",filename.base,"_diversity.tif"), overwrite=FALSE)
writeRaster(alpha1, paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",filename1,"_diversity_",EM,".tif"), overwrite=FALSE)
writeRaster(alpha2, paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",filename2,"_diversity_",EM,".tif"), overwrite=FALSE)
writeRaster(alpha3, paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",filename3,"_diversity_",EM,".tif"), overwrite=FALSE)
writeRaster(alpha4, paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",filename4,"_diversity_",EM,".tif"), overwrite=FALSE)


### Diversity difference
alpha_base <- terra::rast(paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",filename.base,"_diversity.tif"))
alpha1 <- terra::rast(paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",filename1,"_diversity_",EM,".tif"))
alpha2 <- terra::rast(paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",filename2,"_diversity_",EM,".tif"))
alpha3 <- terra::rast(paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",filename3,"_diversity_",EM,".tif"))
alpha4 <- terra::rast(paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",filename4,"_diversity_",EM,".tif"))

alpha_base <- terra::crop(alpha_base,terra::ext(alpha1))

diff.full <- alpha1-alpha_base
#plot(diff.full)
diff.trun1 <- alpha2-alpha_base
#plot(diff.trun1)
diff.trun2 <- alpha3-alpha_base
#plot(diff.trun2)
diff.trun3 <- alpha4-alpha_base
#plot(diff.trun3)

writeRaster(diff.full, paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",filename1,"_diff_",EM,".tif"), overwrite=FALSE)
writeRaster(diff.trun1, paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",filename2,"_diff_",EM,".tif"), overwrite=FALSE)
writeRaster(diff.trun2, paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",filename3,"_diff_",EM,".tif"), overwrite=FALSE)
writeRaster(diff.trun3, paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",filename4,"_diff_",EM,".tif"), overwrite=FALSE)


########## Projected alpha diversity maps ##########

diff.perc.full <- diff.full/alpha_base
diff.perc.trun1 <- diff.trun1/alpha_base
diff.perc.trun2 <- diff.trun2/alpha_base
diff.perc.trun3 <- diff.trun3/alpha_base

writeRaster(diff.perc.full, paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",filename1,"_diffperc_",EM,".tif"), overwrite=FALSE)
writeRaster(diff.perc.trun1, paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",filename2,"_diffperc_",EM,".tif"), overwrite=FALSE)
writeRaster(diff.perc.trun2, paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",filename3,"_diffperc_",EM,".tif"), overwrite=FALSE)
writeRaster(diff.perc.trun3, paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",filename4,"_diffperc_",EM,".tif"), overwrite=FALSE)


########## Stack clamping mask ##########

### Uncertainty pattern
mask_stack1 <- terra::rast(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename1,"/",common.list,projection,"ClampingMask.grd"))
mask_alpha1 <- terra::app(mask_stack1, fun=sum, cores=ncores)
plot(mask_alpha1)

mask_stack2 <- terra::rast(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename2,"/",common.list,projection,"ClampingMask.grd"))
mask_alpha2 <- terra::app(mask_stack2, fun=sum, cores=ncores)
plot(mask_alpha2)

mask_stack3 <- terra::rast(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename3,"/",common.list,projection,"ClampingMask.grd"))
mask_alpha3 <- terra::app(mask_stack3, fun=sum, cores=ncores)
plot(mask_alpha3)

mask_stack4 <- terra::rast(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename4,"/",common.list,projection,"ClampingMask.grd"))
mask_alpha4 <- terra::app(mask_stack4, fun=sum, cores=ncores)
plot(mask_alpha4)


writeRaster(mask_alpha1, paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",filename1,"_uncertaintymask_",EM,".tif"), overwrite=FALSE)
writeRaster(mask_alpha2, paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",filename2,"_uncertaintymask_",EM,".tif"), overwrite=FALSE)
writeRaster(mask_alpha3, paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",filename3,"_uncertaintymask_",EM,".tif"), overwrite=FALSE)
writeRaster(mask_alpha4, paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",filename4,"_uncertaintymask_",EM,".tif"), overwrite=FALSE)


########## FPOS FNEG analysis ##########

# Create df for documentation
columns <- c("Species","Full.FP","Full.FP.perc","Full.FN","Full.FN.perc","Full.F","Full.F.perc",
             "Trun5.FP","Trun5.FP.perc","Trun5.FN","Trun5.FN.perc","Trun5.F","Trun5.F.perc",
             "Trun10.FP","Trun10.FP.perc","Trun10.FN","Trun10.FN.perc","Trun10.F","Trun10.F.perc",
             "Trun15.FP","Trun15.FP.perc","Trun15.FN","Trun15.FN.perc","Trun15.F","Trun15.F.perc")
diff.perc.df <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(diff.perc.df) <- columns

### Do the analysis

# Automated calculation of error for each species
Accuracy_Fun <- function(i) {
  
  # Load maps
  map1 <- terra::rast(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename1,"/",common.list[[i]],projection,common.list[[i]],"_ensemble_TSSbin.tif"))
  map1 <- map1[EM]
  map2 <- terra::rast(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename2,"/",common.list[[i]],projection,common.list[[i]],"_ensemble_TSSbin.tif"))
  map2 <- map2[EM]
  map3 <- terra::rast(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename3,"/",common.list[[i]],projection,common.list[[i]],"_ensemble_TSSbin.tif"))
  map3 <- map3[EM]
  map4 <- terra::rast(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename4,"/",common.list[[i]],projection,common.list[[i]],"_ensemble_TSSbin.tif"))
  map4 <- map4[EM]
  
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


########## Visualise results ##########

### Load output df
diff.perc.df <- readRDS(paste0("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/",save.filename,".RData"))
diff.perc.df <- readRDS("/lustre1/g/sbs_bonebrake/Eugene/Tempdata/Trunlevel_REV_EMca_2025-09-05.RData")

### False positive percentage
plot(density(diff.perc.df$Full.FP.perc,n=10000),xlim = c(-20,250))
lines(density(diff.perc.df$Trun5.FP.perc,n=10000),col="#364B9A")
lines(density(diff.perc.df$Trun10.FP.perc,n=10000),col="#4A7BB7")
lines(density(diff.perc.df$Trun15.FP.perc,n=10000),col="#6EA6CD")

### False positive number
plot(density(diff.perc.df$Full.FP,n=10000),xlim = c(-1000,10000))
lines(density(diff.perc.df$Trun5.FP,n=10000),col="#364B9A")
lines(density(diff.perc.df$Trun10.FP,n=10000),col="#4A7BB7")
lines(density(diff.perc.df$Trun15.FP,n=10000),col="#6EA6CD")


### False negative percentage
plot(density(diff.perc.df$Full.FN.perc,n=10000),xlim = c(-5,70))
lines(density(diff.perc.df$Trun5.FN.perc,n=10000),col="#364B9A")
lines(density(diff.perc.df$Trun10.FN.perc,n=10000),col="#4A7BB7")
lines(density(diff.perc.df$Trun15.FN.perc,n=10000),col="#6EA6CD")

### False negative number
plot(density(diff.perc.df$Full.FN,n=10000),xlim = c(-1000,8000))
lines(density(diff.perc.df$Trun5.FN,n=10000),col="#364B9A")
lines(density(diff.perc.df$Trun10.FN,n=10000),col="#4A7BB7")
lines(density(diff.perc.df$Trun15.FN,n=10000),col="#6EA6CD")


### Overall inaccuracy percentage
plot(density(diff.perc.df$Full.F.perc,n=10000),xlim = c(5,120))
lines(density(diff.perc.df$Trun5.F.perc,n=10000),col="#364B9A")
lines(density(diff.perc.df$Trun10.F.perc,n=10000),col="#4A7BB7")
lines(density(diff.perc.df$Trun15.F.perc,n=10000),col="#6EA6CD")

### Overall inaccuracy number
plot(density(diff.perc.df$Full.F,n=10000),xlim = c(-1000,15000))
lines(density(diff.perc.df$Trun5.F,n=10000),col="#364B9A")
lines(density(diff.perc.df$Trun10.F,n=10000),col="#4A7BB7")
lines(density(diff.perc.df$Trun15.F,n=10000),col="#6EA6CD")


########## Stastistically test results ##########

ks.test(diff.perc.df$Full.FN,
        diff.perc.df$Trun15.FN,
        alternative = "two.sided",exact = NULL, 
        simulate.p.value=FALSE)

ks.test(diff.perc.df$Trun5.FN,
        diff.perc.df$Trun15.FN,
        alternative = "two.sided",exact = NULL, 
        simulate.p.value=FALSE)
