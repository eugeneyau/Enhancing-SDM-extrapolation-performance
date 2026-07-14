#library(dplyr)

########## Produce Variable Importance Table ###########


### Virtual species ###

filepath <- "biomod424_BASE_VIF5"
basepath <- paste0("/lustre1/g/sbs_bonebrake/Eugene/SDM_BASE/results_",filepath)
BASEprojection <- "/proj_Current/proj_Current_"

species.list.out <- list.files(basepath)
poss <- paste0(basepath,"/",species.list.out,BASEprojection,species.list.out,"_ensemble_TSSbin.tif")
tf <- file.exists(poss)
sp.list <- species.list.out[tf]
sp_list <- gsub("\\."," ",sp.list)

# Producing var imp table by adding rows
varimp <- read.csv(paste0(basepath, "/",sp.list[[1]],"/",sp_list[[1]],"_formal_models_variables_importance.csv"))
for(i in 2:length(sp.list)) {
  new <- read.csv(paste0(basepath, "/",sp.list[[i]],"/",sp_list[[i]],"_formal_models_variables_importance.csv"))
  varimp <- rbind(varimp, new)
}
species.name <- gsub("_.*","",varimp$full.name)
varimp$full.name <- species.name

# Get mean var importance for each species
listed.sp <- unique(species.name)
var.list <- unique(varimp$expl.var)

columns <- c("Species","Variable","SPmean") 
varimp.weighted <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(varimp.weighted) <- columns
for(x in 1:length(listed.sp)) { sp.data <- subset(varimp, varimp$full.name==listed.sp[[x]])
for(y in 1:length(var.list)) {
  sp.vardata <- subset(sp.data, sp.data$expl.var==var.list[[y]])
  SPmean <- mean(sp.vardata$var.imp, na.rm=TRUE)
  varimp.weighted <- rbind(varimp.weighted, data.frame("Species"=listed.sp[[x]], "Variable"=var.list[[y]], "SPmean"=SPmean))
}}

# Summarise (mean of mean, median of mean, sd of mean) varaible imp for all species 
columns <- c("Variable","Importance.mean", "Importance.median", "Importance.sd") 
varimp.breakdown <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(varimp.breakdown) <- columns

for(i in 1:length(var.list)) {
  data <- subset(varimp.weighted, varimp.weighted$Variable==var.list[[i]])
  Imp.mean <- mean(data$SPmean, na.rm=TRUE)
  Imp.sd <- sd(data$SPmean, na.rm=TRUE)
  Importance.median <- median(data$SPmean, na.rm=TRUE)
  varimp.breakdown <- rbind(varimp.breakdown, data.frame("Variable"=var.list[[i]], "Importance.mean"=Imp.mean, "Importance.median"=Importance.median, "Importance.sd"=Imp.sd))
}

varimp.breakdown

                       Variable Importance.mean Importance.median Importance.sd
1               WorldClim_bio04      0.45051264        0.44371887    0.16892352
2               WorldClim_bio05      0.08210124        0.06257169    0.06788813
3               WorldClim_bio12      0.14471552        0.11203727    0.10312002
4               WorldClim_bio14      0.09427080        0.07576918    0.06329786
5               WorldClim_bio15      0.07358022        0.05870793    0.05528694
6                      NDVImean      0.05075679        0.03798134    0.04681380
7                 Canopy_Height      0.10045851        0.08055759    0.08016126
8 Soilgrids_soil_organic_carbon      0.06538010        0.05584953    0.03823088
9      Soilgrids_total_nitrogen      0.05082278        0.03991583    0.03962405


### SDM output ###

filepath <- "bmod424_trunc_varnoise_VIF5_SOUTH_trun20perc"
#filepath <- "bmod424_trunc_VIF5_SOUTH_trun20perc"
basepath <- paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filepath)
BASEprojection <- "/proj_TSSNA/proj_TSSNA_"

species.list.out <- list.files(basepath)
poss <- paste0(basepath,"/",species.list.out,BASEprojection,species.list.out,"_ensemble_TSSbin.tif")
tf <- file.exists(poss)
sp.list <- species.list.out[tf]
sp_list <- gsub("\\."," ",sp.list)

# Producing var imp table by adding rows
varimp <- read.csv(paste0(basepath, "/",sp.list[[1]],"/",sp.list[[1]],"_formal_models_variables_importance.csv"))
for(i in 2:length(sp.list)) {
  new <- read.csv(paste0(basepath, "/",sp.list[[i]],"/",sp.list[[i]],"_formal_models_variables_importance.csv"))
  varimp <- rbind(varimp, new)
}
species.name <- gsub("_.*","",varimp$full.name)
varimp$full.name <- species.name

# Get mean var importance for each species
listed.sp <- unique(species.name)
var.list <- unique(varimp$expl.var)

columns <- c("Species","Variable","SPmean") 
varimp.weighted <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(varimp.weighted) <- columns
for(x in 1:length(listed.sp)) { sp.data <- subset(varimp, varimp$full.name==listed.sp[[x]])
for(y in 1:length(var.list)) {
  sp.vardata <- subset(sp.data, sp.data$expl.var==var.list[[y]])
  SPmean <- mean(sp.vardata$var.imp, na.rm=TRUE)
  varimp.weighted <- rbind(varimp.weighted, data.frame("Species"=listed.sp[[x]], "Variable"=var.list[[y]], "SPmean"=SPmean))
}}

# Summarise varaible imp for all species
columns <- c("Variable","Importance.mean", "Importance.median", "Importance.sd") 
varimp.breakdown <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(varimp.breakdown) <- columns

for(i in 1:length(var.list)) {
  data <- subset(varimp.weighted, varimp.weighted$Variable==var.list[[i]])
  Imp.mean <- mean(data$SPmean, na.rm=TRUE)
  Imp.sd <- sd(data$SPmean, na.rm=TRUE)
  Importance.median <- median(data$SPmean, na.rm=TRUE)
  varimp.breakdown <- rbind(varimp.breakdown, data.frame("Variable"=var.list[[i]], "Importance.mean"=Imp.mean, "Importance.median"=Importance.median, "Importance.sd"=Imp.sd))
}

varimp.breakdown

"bmod424_trunc_varnoise_VIF5_SOUTH_trun20perc"
                             Variable Importance.mean Importance.median Importance.sd
1                     WorldClim_bio04     0.413866780       0.421699977   0.248554924
2                     WorldClim_bio05     0.061478481       0.029927174   0.077256443
3                     WorldClim_bio12     0.145676681       0.087222004   0.152468596
4                     WorldClim_bio14     0.077293587       0.029470973   0.114275367
5                     WorldClim_bio15     0.050253228       0.026397002   0.064189050
6                            NDVImean     0.041074988       0.019524656   0.065556273
7                       Canopy_Height     0.090592963       0.051011822   0.106437792
8       Soilgrids_soil_organic_carbon     0.032913368       0.021345161   0.040769548
9            Soilgrids_total_nitrogen     0.030766589       0.019590483   0.037621436
10 Soilgrids_cation_exchange_capacity     0.013890589       0.009665617   0.013034065
11         Soilgrids_coarse_fragments     0.014563010       0.011035600   0.012321471
12                     Soilgrids_sand     0.009784963       0.007654700   0.007676871
13                Soilgrids_soil_silt     0.010121702       0.008414279   0.007336779

"bmod424_trunc_VIF5_SOUTH_trun20perc"
                       Variable Importance.mean Importance.median Importance.sd
1               WorldClim_bio04      0.42491740        0.43244367    0.24876600
2               WorldClim_bio05      0.07061840        0.03821048    0.08122285
3               WorldClim_bio12      0.15618517        0.09839435    0.15843288
4               WorldClim_bio14      0.08465620        0.03391715    0.12019168
5               WorldClim_bio15      0.05861966        0.03208270    0.07183116
6                      NDVImean      0.04496516        0.02266711    0.06677404
7                 Canopy_Height      0.09976162        0.06400840    0.10901628
8 Soilgrids_soil_organic_carbon      0.03821319        0.02512742    0.04580416
9      Soilgrids_total_nitrogen      0.03438498        0.02284287    0.04198774



#write.csv(varimp, file.path(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/A_", filepath, "_variables_importance.csv")))


