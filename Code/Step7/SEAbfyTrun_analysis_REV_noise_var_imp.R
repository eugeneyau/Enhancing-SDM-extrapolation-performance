#library(dplyr)

########## Produce Variable Importance Table ###########


### Virtual species ###

filepath <- "biomod424_BASE_REV"
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

                        Variable Importance.mean Importance.median Importance.sd
1                WorldClim_bio01      0.09837335        0.06702474    0.09279478
2                WorldClim_bio04      0.30505821        0.26276491    0.16645432
3                WorldClim_bio05      0.05410877        0.04264950    0.04351210
4                WorldClim_bio06      0.15474720        0.12891267    0.11165036
5                WorldClim_bio12      0.10033447        0.08689469    0.05993580
6                WorldClim_bio13      0.06388941        0.04382835    0.06173088
7                WorldClim_bio14      0.07345707        0.06011292    0.05014602
8                WorldClim_bio15      0.05956143        0.04823076    0.04509388
9                       NDVImean      0.03839830        0.02715324    0.03911292
10                 Canopy_Height      0.07043672        0.05373891    0.05840215
11             Soilgrids_soil_pH      0.07592832        0.05326312    0.06449548
12 Soilgrids_soil_organic_carbon      0.05250929        0.04600407    0.02952119
13      Soilgrids_total_nitrogen      0.04113706        0.03158127    0.03298534




### SDM output ###

filepath <- "bmod424_trunc_varnoise_rev_SOUTH_trun20perc"
#filepath <- "bmod424_trunc_rev_SOUTH_trun20perc"
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

"bmod424_trunc_varnoise_rev_SOUTH_trun20perc"
                             Variable Importance.mean Importance.median Importance.sd
1                     WorldClim_bio01     0.057424625       0.041588245   0.057463542
2                     WorldClim_bio04     0.310995476       0.250417716   0.236326089
3                     WorldClim_bio05     0.039524290       0.024691437   0.043419955
4                     WorldClim_bio06     0.120011327       0.072937610   0.126350733
5                     WorldClim_bio12     0.097065834       0.059046216   0.106457591
6                     WorldClim_bio13     0.037191019       0.024605105   0.047946323
7                     WorldClim_bio14     0.063417778       0.025983845   0.096936627
8                     WorldClim_bio15     0.043370783       0.022799445   0.058833605
9                            NDVImean     0.030858031       0.013811320   0.053519654
10                      Canopy_Height     0.059505295       0.029339748   0.075040149
11                  Soilgrids_soil_pH     0.091394180       0.046239788   0.105484321
12      Soilgrids_soil_organic_carbon     0.026701348       0.018970208   0.028288406
13           Soilgrids_total_nitrogen     0.025324139       0.017253492   0.027304034
14 Soilgrids_cation_exchange_capacity     0.011739610       0.008856117   0.010236602
15         Soilgrids_coarse_fragments     0.012962368       0.009830524   0.011433024
16                     Soilgrids_sand     0.008700338       0.006645694   0.007015520
17                Soilgrids_soil_silt     0.008720327       0.006884131   0.007035937


"bmod424_trunc_rev_SOUTH_trun20perc"
                        Variable Importance.mean Importance.median Importance.sd
1                WorldClim_bio01      0.06345497        0.04705455    0.06278627
2                WorldClim_bio04      0.31521917        0.25859488    0.23242967
3                WorldClim_bio05      0.04369770        0.02834054    0.04613513
4                WorldClim_bio06      0.12774238        0.08072867    0.13219402
5                WorldClim_bio12      0.09928616        0.06313438    0.10534783
6                WorldClim_bio13      0.03814801        0.02407628    0.04639277
7                WorldClim_bio14      0.06545899        0.02825402    0.09778422
8                WorldClim_bio15      0.04748737        0.02523324    0.06110740
9                       NDVImean      0.03208419        0.01419339    0.05484398
10                 Canopy_Height      0.06258506        0.03109697    0.07901109
11             Soilgrids_soil_pH      0.09629734        0.04767183    0.10937315
12 Soilgrids_soil_organic_carbon      0.02883777        0.02093270    0.02956849
13      Soilgrids_total_nitrogen      0.02824239        0.01846818    0.03228194




#write.csv(varimp, file.path(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/A_", filepath, "_variables_importance.csv")))



















