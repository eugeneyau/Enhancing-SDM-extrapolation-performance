library(dplyr)

########## Produce Variable Importance Table ###########
filepath <- "biomod424_BASE_69-161.6_-10-36"
basepath <- paste0("/lustre1/g/sbs_bonebrake/Eugene/SDM_BASE/results_",filepath)
species.list.out <- list.files(basepath)

BASEprojection <- "/proj_Current/proj_Current_"
poss <- paste0(basepath,"/",species.list.out,BASEprojection,species.list.out,"_ensemble_TSSbin.tif")
tf <- file.exists(poss)
sp.list <- species.list.out[tf]
sp_list <- gsub("\\."," ",sp.list)


# Producing var imp table by adding rows
varimp <- read.csv(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDM_BASE/results_", filepath, "/",sp.list[[1]],"/",sp_list[[1]],"_formal_models_variables_importance.csv"))
for(i in 2:length(sp.list)) {
  new <- read.csv(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDM_BASE/results_", filepath, "/",sp.list[[i]],"/",sp_list[[i]],"_formal_models_variables_importance.csv"))
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
#write.csv(varimp, file.path(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/A_", filepath, "_variables_importance.csv")))

#Filter by algorithm
algorithm <- "MAXNET" #From "GLM","MAXNET","SRE","MARS"
var.list <- unique(varimp$expl.var)

columns <- c("Variable","Importance.mean", "Importance.sd") 
varimp.breakdown <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(varimp.breakdown) <- columns

for(i in 1:length(var.list)) {
  data <- subset(varimp,expl.var==var.list[[i]] & algo==algorithm)
  Imp.mean <- mean(data$var.imp, na.rm=TRUE)
  Imp.sd <- sd(data$var.imp, na.rm=TRUE)
  varimp.breakdown = rbind(varimp.breakdown, data.frame("Variable"=var.list[[i]], "Importance.mean"=Imp.mean, "Importance.sd"=Imp.sd))
}

varimp.breakdown
