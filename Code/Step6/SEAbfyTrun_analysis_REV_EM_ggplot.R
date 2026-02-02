extrafont::loadfonts(device="win")
extrafont::loadfonts(device="postscript")
library(ggplot2)
library(ggpubr)

#wd <- "/lustre1/g/sbs_bonebrake/Eugene/Tempdata"  # HPC
wdpc <- "C:/Users/Eugene/Desktop/Paper1_Vis"  # PC
wdG15 <- "C:/Users/skywa/Desktop/Paper1_Vis"  # Laptop

setwd(wdpc)

options(scipen = 999)

save.filename <- "EM_TSSNA_REV_Trunc20" # Name the df output file

diff.perc.df <- readRDS(paste0(save.filename,".RData"))


########## Format data for ggplot ##########

# Create new DF 
columns <- c("Species","setup",
             "FP","FP.perc","FN","FN.perc","F","F.perc")
ggplot.df <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(ggplot.df) <- columns

# Write new DF
Accuracy_Fun <- function(i) {
  
  sp.data <- diff.perc.df[i,]
  species <- sp.data$Species
  
  ###
  EMca.setup <- "EMca"
  EMca.FP <- sp.data$EMca.FP
  EMca.FP.perc <- sp.data$EMca.FP.perc
  EMca.FN <- sp.data$EMca.FN
  EMca.FN.perc <- sp.data$EMca.FN.perc
  EMca.F <- sp.data$EMca.F
  EMca.F.perc <- sp.data$EMca.F.perc
  
  EMca.data <- list(species,EMca.setup,EMca.FP,EMca.FP.perc,EMca.FN,EMca.FN.perc,EMca.F,EMca.F.perc)
  
  ggplot.df[nrow(ggplot.df)+1,] <<- EMca.data
  
  ###
  EMmean.setup <- "EMmean"
  EMmean.FP <- sp.data$EMmean.FP
  EMmean.FP.perc <- sp.data$EMmean.FP.perc
  EMmean.FN <- sp.data$EMmean.FN
  EMmean.FN.perc <- sp.data$EMmean.FN.perc
  EMmean.F <- sp.data$EMmean.F
  EMmean.F.perc <- sp.data$EMmean.F.perc
  
  EMmean.data <- list(species,EMmean.setup,EMmean.FP,EMmean.FP.perc,EMmean.FN,EMmean.FN.perc,EMmean.F,EMmean.F.perc)
  
  ggplot.df[nrow(ggplot.df)+1,] <<- EMmean.data
  
  ###
  EMwmean.setup <- "EMwmean"
  EMwmean.FP <- sp.data$EMwmean.FP
  EMwmean.FP.perc <- sp.data$EMwmean.FP.perc
  EMwmean.FN <- sp.data$EMwmean.FN
  EMwmean.FN.perc <- sp.data$EMwmean.FN.perc
  EMwmean.F <- sp.data$EMwmean.F
  EMwmean.F.perc <- sp.data$EMwmean.F.perc
  
  EMwmean.data <- list(species,EMwmean.setup,EMwmean.FP,EMwmean.FP.perc,EMwmean.FN,EMwmean.FN.perc,EMwmean.F,EMwmean.F.perc)
  
  ggplot.df[nrow(ggplot.df)+1,] <<- EMwmean.data
  
  ###
  EMmedian.setup <- "EMmedian"
  EMmedian.FP <- sp.data$EMmedian.FP
  EMmedian.FP.perc <- sp.data$EMmedian.FP.perc
  EMmedian.FN <- sp.data$EMmedian.FN
  EMmedian.FN.perc <- sp.data$EMmedian.FN.perc
  EMmedian.F <- sp.data$EMmedian.F
  EMmedian.F.perc <- sp.data$EMmedian.F.perc
  
  EMmedian.data <- list(species,EMmedian.setup,EMmedian.FP,EMmedian.FP.perc,EMmedian.FN,EMmedian.FN.perc,EMmedian.F,EMmedian.F.perc)
  
  ggplot.df[nrow(ggplot.df)+1,] <<- EMmedian.data
  
}

i=1:nrow(diff.perc.df)

lapply(i,Accuracy_Fun)






###### Count results ######

### False 

Fplot <- ggplot(ggplot.df) + 
  stat_density(aes(x=F, colour=setup), geom="line", position="identity", linewidth=0.55) +
  #stat_density(aes(x=t.F, colour=setup), geom="line", linetype=5, position="identity", size=0.55) +
  #stat_density(aes(x=f.F, colour=setup), geom="line", position="identity", size=0.55) +
  geom_hline(yintercept = 0, color = "dark grey", linewidth=0.5) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth=0.7),
        ###
        plot.margin = margin(0.5, 0.5, 0.6, 0.5, unit = "cm"),
        plot.title    = element_text(size=22, family = "Calibri", face = "bold", hjust = 0.5, vjust = 1),
        plot.subtitle = element_text(size=20, family = "Calibri"),
        ###
        axis.title.x  = element_text(size=17, family = "Calibri", vjust=-0.3),
        axis.title.y  = element_text(size=17, family = "Calibri", vjust=1.5),
        axis.text.x   = element_text(size=12.5, family = "Calibri"),
        axis.text.y   = element_blank(),
        axis.ticks.y  = element_blank(),
        ###
        axis.ticks = element_line(linewidth = 0.6) , 
        axis.ticks.length = unit(0.2, "cm"),
        ###
        legend.key = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank(),
        legend.position = "none") +
  labs(x = "Total false predictions", y = "Density", fill="xyz") +
  scale_color_manual(values=c("#999933","#0077BB","#AA4499","#33BBEE")) +
  scale_x_continuous(limits = c(-2000,19000),expand=c(0,0), breaks = seq(0,15000,5000)) #+
  #scale_y_continuous(expand=c(0.025,0)) # c(a,b): a=relative expansion(%) b=absolute expansion


### False positive 

FPplot <- ggplot(ggplot.df) + 
  stat_density(aes(x=FP, colour=setup), geom="line", position="identity", linewidth=0.55) +
  #stat_density(aes(x=t.F, colour=setup), geom="line", linetype=5, position="identity", size=0.55) +
  #stat_density(aes(x=f.F, colour=setup), geom="line", position="identity", size=0.55) +
  geom_hline(yintercept = 0, color = "dark grey", linewidth=0.5) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth=0.7),
        ###
        plot.margin = margin(0.5, 0.5, 0.6, 0.5, unit = "cm"),
        plot.title    = element_text(size=22, family = "Calibri", face = "bold", hjust = 0.5, vjust = 1),
        plot.subtitle = element_text(size=20, family = "Calibri"),
        axis.title.x  = element_text(size=17, family = "Calibri", vjust=-0.3),
        axis.title.y  = element_text(size=17, family = "Calibri", vjust=1.5),
        axis.text.x   = element_text(size=12.5, family = "Calibri"),
        axis.text.y   = element_blank(),
        axis.ticks.y  = element_blank(),
        ###
        axis.ticks = element_line(linewidth = 0.6) , 
        axis.ticks.length = unit(0.2, "cm"),
        ###
        legend.key = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank(),
        legend.position = "none") +
  labs(x = "False positive predictions", y = "Density", fill="xyz") +
  scale_color_manual(values=c("#999933","#0077BB","#AA4499","#33BBEE")) +
  scale_x_continuous(limits = c(-1700,11000),expand=c(0,0), breaks = seq(0,10000,2500)) #+
  #scale_y_continuous(expand=c(0.025,0)) # c(a,b): a=relative expansion(%) b=absolute expansion


### False negative 

FNplot <- ggplot(ggplot.df) + 
  stat_density(aes(x=FN, colour=setup), geom="line", position="identity", linewidth=0.55) +
  #stat_density(aes(x=t.F, colour=setup), geom="line", linetype=5, position="identity", size=0.55) +
  #stat_density(aes(x=f.F, colour=setup), geom="line", position="identity", size=0.55) +
  geom_hline(yintercept = 0, color = "dark grey", linewidth=0.5) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth=0.7),
        ###
        plot.margin = margin(0.5, 0.5, 0.6, 0.5, unit = "cm"),
        plot.title    = element_text(size=22, family = "Calibri", face = "bold", hjust = 0.5, vjust = 1),
        plot.subtitle = element_text(size=20, family = "Calibri"),
        axis.title.x  = element_text(size=17, family = "Calibri", vjust=-0.3),
        axis.title.y  = element_text(size=17, family = "Calibri", vjust=1.5),
        axis.text.x   = element_text(size=12.5, family = "Calibri"),
        axis.text.y   = element_blank(),
        axis.ticks.y  = element_blank(),
        ###
        axis.ticks = element_line(linewidth = 0.6) , 
        axis.ticks.length = unit(0.2, "cm"),
        ###
        legend.key = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank(),
        legend.position = "none") +
  labs(x = "False negative predictions", y = "Density", fill="xyz") +
  scale_color_manual(values=c("#999933","#0077BB","#AA4499","#33BBEE")) +
  scale_x_continuous(limits = c(-1200,8000),expand=c(0,0), breaks = seq(0,8500,2500)) #+
  #scale_y_continuous(expand=c(0.025,0)) # c(a,b): a=relative expansion(%) b=absolute expansion


ggarrange(Fplot, FPplot, FNplot,
          nrow = 3, ncol = 1, align = "v")

ggsave("EM.png", width = 24.5, height = 35, units = "cm", dpi = 1000)




###### Legend ######
ggplot(ggplot.df) + 
  stat_density(aes(x=FN.perc, colour=setup), geom="line", position="identity", size=6) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        ###
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm"),
        axis.title.y  = element_blank(),
        axis.text.y   = element_blank(),
        axis.ticks.y  = element_blank(),
        legend.key = element_blank(),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.text = element_text(size = 14, family = "Calibri"),
        legend.title = element_text(size = 14, family = "Calibri"
                                    # , margin = margin(b=20)
        ),
        legend.key.width = unit(1,"cm"),
        legend.key.spacing.x = unit(0.5,'cm')) +
  labs(y = "Density", x = "False predictions", fill="xyz") +
  guides(color = guide_legend(title = "Ensemble algorithm      ")) +
  scale_color_manual(values=c("#999933","#0077BB","#AA4499","#33BBEE"), 
                     labels=c("Committee averaging","Mean","Median","Weighted mean"))

ggsave("EM_legend.png", width = 26, height = 12, units = "cm", dpi = 1200)








###### Bar plot ######
Fthresh <- 5000
FPthresh <- 3000
FNthresh <- 1500

### Create new DF 
columns <- c("setup","compare","t.F.sp","t.F.sp.perc","t.FP.sp","t.FP.sp.perc","t.FN.sp","t.FN.sp.perc")
ggbar.df <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(ggbar.df) <- columns

### Specify names
setup1 <- "EMca"
name1 <- "EMca"
setup2 <- "EMmean"
name2 <- "EMmean"
setup3 <- "EMwmean"
name3 <- "EMwmean"
setup4 <- "EMmedian"
name4 <- "EMmedian"

### Get data
if(!is.na(setup1)) {
  setup1.F.more <- length(diff.perc.df[[paste0(setup1, ".F")]][diff.perc.df[[paste0(setup1, ".F")]] >= Fthresh])
  setup1.F.more.perc <- setup1.F.more / length(diff.perc.df[[paste0(setup1, ".F")]])*100
  setup1.F.less <- length(diff.perc.df[[paste0(setup1, ".F")]][diff.perc.df[[paste0(setup1, ".F")]] < Fthresh])
  setup1.F.less.perc <- setup1.F.less / length(diff.perc.df[[paste0(setup1, ".F")]])*100
  
  setup1.FP.more <- length(diff.perc.df[[paste0(setup1, ".FP")]][diff.perc.df[[paste0(setup1, ".FP")]] >= FPthresh])
  setup1.FP.more.perc <- setup1.FP.more / length(diff.perc.df[[paste0(setup1, ".FP")]])*100
  setup1.FP.less <- length(diff.perc.df[[paste0(setup1, ".FP")]][diff.perc.df[[paste0(setup1, ".FP")]] < FPthresh])
  setup1.FP.less.perc <- setup1.FP.less / length(diff.perc.df[[paste0(setup1, ".FP")]])*100
  
  setup1.FN.more <- length(diff.perc.df[[paste0(setup1, ".FN")]][diff.perc.df[[paste0(setup1, ".FN")]] >= FNthresh])
  setup1.FN.more.perc <- setup1.FN.more / length(diff.perc.df[[paste0(setup1, ".FN")]])*100
  setup1.FN.less <- length(diff.perc.df[[paste0(setup1, ".FN")]][diff.perc.df[[paste0(setup1, ".FN")]] < FNthresh])
  setup1.FN.less.perc <- setup1.FN.less / length(diff.perc.df[[paste0(setup1, ".FN")]])*100
} 

if(!is.na(setup2)) {
  setup2.F.more <- length(diff.perc.df[[paste0(setup2, ".F")]][diff.perc.df[[paste0(setup2, ".F")]] >= Fthresh])
  setup2.F.more.perc <- setup2.F.more / length(diff.perc.df[[paste0(setup2, ".F")]])*100
  setup2.F.less <- length(diff.perc.df[[paste0(setup2, ".F")]][diff.perc.df[[paste0(setup2, ".F")]] < Fthresh])
  setup2.F.less.perc <- setup2.F.less / length(diff.perc.df[[paste0(setup2, ".F")]])*100
  
  setup2.FP.more <- length(diff.perc.df[[paste0(setup2, ".FP")]][diff.perc.df[[paste0(setup2, ".FP")]] >= FPthresh])
  setup2.FP.more.perc <- setup2.FP.more / length(diff.perc.df[[paste0(setup2, ".FP")]])*100
  setup2.FP.less <- length(diff.perc.df[[paste0(setup2, ".FP")]][diff.perc.df[[paste0(setup2, ".FP")]] < FPthresh])
  setup2.FP.less.perc <- setup2.FP.less / length(diff.perc.df[[paste0(setup2, ".FP")]])*100
  
  setup2.FN.more <- length(diff.perc.df[[paste0(setup2, ".FN")]][diff.perc.df[[paste0(setup2, ".FN")]] >= FNthresh])
  setup2.FN.more.perc <- setup2.FN.more / length(diff.perc.df[[paste0(setup2, ".FN")]])*100
  setup2.FN.less <- length(diff.perc.df[[paste0(setup2, ".FN")]][diff.perc.df[[paste0(setup2, ".FN")]] < FNthresh])
  setup2.FN.less.perc <- setup2.FN.less / length(diff.perc.df[[paste0(setup2, ".FN")]])*100
}

if(!is.na(setup3)) {
  setup3.F.more <- length(diff.perc.df[[paste0(setup3, ".F")]][diff.perc.df[[paste0(setup3, ".F")]] >= Fthresh])
  setup3.F.more.perc <- setup3.F.more / length(diff.perc.df[[paste0(setup3, ".F")]])*100
  setup3.F.less <- length(diff.perc.df[[paste0(setup3, ".F")]][diff.perc.df[[paste0(setup3, ".F")]] < Fthresh])
  setup3.F.less.perc <- setup3.F.less / length(diff.perc.df[[paste0(setup3, ".F")]])*100
  
  setup3.FP.more <- length(diff.perc.df[[paste0(setup3, ".FP")]][diff.perc.df[[paste0(setup3, ".FP")]] >= FPthresh])
  setup3.FP.more.perc <- setup3.FP.more / length(diff.perc.df[[paste0(setup3, ".FP")]])*100
  setup3.FP.less <- length(diff.perc.df[[paste0(setup3, ".FP")]][diff.perc.df[[paste0(setup3, ".FP")]] < FPthresh])
  setup3.FP.less.perc <- setup3.FP.less / length(diff.perc.df[[paste0(setup3, ".FP")]])*100
  
  setup3.FN.more <- length(diff.perc.df[[paste0(setup3, ".FN")]][diff.perc.df[[paste0(setup3, ".FN")]] >= FNthresh])
  setup3.FN.more.perc <- setup3.FN.more / length(diff.perc.df[[paste0(setup3, ".FN")]])*100
  setup3.FN.less <- length(diff.perc.df[[paste0(setup3, ".FN")]][diff.perc.df[[paste0(setup3, ".FN")]] < FNthresh])
  setup3.FN.less.perc <- setup3.FN.less / length(diff.perc.df[[paste0(setup3, ".FN")]])*100
}

if(!is.na(setup4)) {
  setup4.F.more <- length(diff.perc.df[[paste0(setup4, ".F")]][diff.perc.df[[paste0(setup4, ".F")]] >= Fthresh])
  setup4.F.more.perc <- setup4.F.more / length(diff.perc.df[[paste0(setup4, ".F")]])*100
  setup4.F.less <- length(diff.perc.df[[paste0(setup4, ".F")]][diff.perc.df[[paste0(setup4, ".F")]] < Fthresh])
  setup4.F.less.perc <- setup4.F.less / length(diff.perc.df[[paste0(setup4, ".F")]])*100
  
  setup4.FP.more <- length(diff.perc.df[[paste0(setup4, ".FP")]][diff.perc.df[[paste0(setup4, ".FP")]] >= FPthresh])
  setup4.FP.more.perc <- setup4.FP.more / length(diff.perc.df[[paste0(setup4, ".FP")]])*100
  setup4.FP.less <- length(diff.perc.df[[paste0(setup4, ".FP")]][diff.perc.df[[paste0(setup4, ".FP")]] < FPthresh])
  setup4.FP.less.perc <- setup4.FP.less / length(diff.perc.df[[paste0(setup4, ".FP")]])*100
  
  setup4.FN.more <- length(diff.perc.df[[paste0(setup4, ".FN")]][diff.perc.df[[paste0(setup4, ".FN")]] >= FNthresh])
  setup4.FN.more.perc <- setup4.FN.more / length(diff.perc.df[[paste0(setup4, ".FN")]])*100
  setup4.FN.less <- length(diff.perc.df[[paste0(setup4, ".FN")]][diff.perc.df[[paste0(setup4, ".FN")]] < FNthresh])
  setup4.FN.less.perc <- setup4.FN.less / length(diff.perc.df[[paste0(setup4, ".FN")]])*100
}


### Write data
if(!is.na(setup1)) {
  setup1.more.data <- list(name1,"more",setup1.F.more,setup1.F.more.perc,setup1.FP.more,setup1.FP.more.perc,setup1.FN.more,setup1.FN.more.perc)
  setup1.less.data <- list(name1,"less",setup1.F.less,setup1.F.less.perc,setup1.FP.less,setup1.FP.less.perc,setup1.FN.less,setup1.FN.less.perc)
  ggbar.df[1,] <- setup1.more.data
  ggbar.df[2,] <- setup1.less.data
}

if(!is.na(setup2)) {
  setup2.more.data <- list(name2,"more",setup2.F.more,setup2.F.more.perc,setup2.FP.more,setup2.FP.more.perc,setup2.FN.more,setup2.FN.more.perc)
  setup2.less.data <- list(name2,"less",setup2.F.less,setup2.F.less.perc,setup2.FP.less,setup2.FP.less.perc,setup2.FN.less,setup2.FN.less.perc)
  ggbar.df[3,] <- setup2.more.data
  ggbar.df[4,] <- setup2.less.data
}

if(!is.na(setup3)) {
  setup3.more.data <- list(name3,"more",setup3.F.more,setup3.F.more.perc,setup3.FP.more,setup3.FP.more.perc,setup3.FN.more,setup3.FN.more.perc)
  setup3.less.data <- list(name3,"less",setup3.F.less,setup3.F.less.perc,setup3.FP.less,setup3.FP.less.perc,setup3.FN.less,setup3.FN.less.perc)
  ggbar.df[5,] <- setup3.more.data
  ggbar.df[6,] <- setup3.less.data
}

if(!is.na(setup4)) {
  setup4.more.data <- list(name4,"more",setup4.F.more,setup4.F.more.perc,setup4.FP.more,setup4.FP.more.perc,setup4.FN.more,setup4.FN.more.perc)
  setup4.less.data <- list(name4,"less",setup4.F.less,setup4.F.less.perc,setup4.FP.less,setup4.FP.less.perc,setup4.FN.less,setup4.FN.less.perc)
  ggbar.df[7,] <- setup4.more.data
  ggbar.df[8,] <- setup4.less.data
}


ggbar.df[,4] <- round(ggbar.df[,4], digits = 0)
ggbar.df[,6] <- round(ggbar.df[,6], digits = 0)
ggbar.df[,8] <- round(ggbar.df[,8], digits = 0)
if(ncol(ggbar.df)>=10) {
  ggbar.df[,10] <- round(ggbar.df[,10], digits = 0)}

#outline <- rep(100, times=nrow(ggbar.df))
ggbar.df$outline <- rep(100, times=nrow(ggbar.df))

library(plyr)
ggbar.df <- ddply(ggbar.df, .(setup), transform, pos.t.F.sp.perc = 100 - (cumsum(t.F.sp.perc) - (0.5 * t.F.sp.perc)))
ggbar.df <- ddply(ggbar.df, .(setup), transform, pos.t.FP.sp.perc = 100 - (cumsum(t.FP.sp.perc) - (0.5 * t.FP.sp.perc)))
ggbar.df <- ddply(ggbar.df, .(setup), transform, pos.t.FN.sp.perc = 100 - (cumsum(t.FN.sp.perc) - (0.5 * t.FN.sp.perc)))
if(ncol(ggbar.df)>=14) {
  ggbar.df <- ddply(ggbar.df, .(setup), transform, pos.f.F.sp.perc = 100 - (cumsum(f.F.sp.perc) - (0.5 * f.F.sp.perc)))}





### F
ggplot(ggbar.df, aes(x=setup, y=t.F.sp.perc)) +
  geom_bar(data = subset(ggbar.df, compare == "less"), aes(fill = setup), 
           stat = "identity", position = "stack", width = 0.75) +
  #geom_bar(data = subset(ggbar.df, compare == "more"), aes(y = outline, color = setup), 
  #         stat = "identity", fill = NA, size = 0.3, color = "black", width = 0.75) +
  geom_bar(data = subset(ggbar.df, compare == "more"), aes(y = outline, color = setup), 
           stat = "identity", fill = NA, size = 0.3, color = c("#999933","#0077BB","#AA4499","#33BBEE"), width = 0.75) + # Coloured version
  geom_text(data = subset(ggbar.df, compare == "less"), aes(x = setup, y = pos.t.F.sp.perc, label = paste0(t.F.sp.perc,"%")),
            colour="white", family="Calibri", size=4.5, angle = 270) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        ###
        plot.margin = margin(0.5, 0.5, 0.6, 0.5, unit = "cm"),
        plot.title    = element_text(size=22, family = "Calibri", face = "bold", hjust = 0.5, vjust = 1),
        plot.subtitle = element_text(size=20, family = "Calibri"),
        axis.title.x  = element_text(size=17, family = "Calibri", vjust=-0.3),
        axis.title.y  = element_text(size=17, family = "Calibri", vjust=2.5),
        axis.text.x   = element_text(size=12.5, family = "Calibri"),
        axis.text.y   = element_blank(),
        ###
        axis.ticks = element_blank(), 
        legend.key = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank(),
        legend.position = "none") +
  labs(x="TSS threshold", y=paste0("Percentage of species with total
false prediction less than ",Fthresh," cells"), fill="xyz") +
  scale_fill_manual(values = c("#999933","#0077BB","#AA4499","#33BBEE")) +
  scale_color_manual(values = c("#999933","#0077BB","#AA4499","#33BBEE")) +
  scale_y_continuous(position = "right", expand=c(0.005,0))

ggsave("EMbar_F.png", width = 5.65, height = 11.7, units = "cm", dpi = 800)


### FP
ggplot(ggbar.df, aes(x=setup, y=t.FP.sp.perc)) +
  geom_bar(data = subset(ggbar.df, compare == "less"), aes(fill = setup), 
           stat = "identity", position = "stack", width = 0.75) +
  #geom_bar(data = subset(ggbar.df, compare == "more"), aes(y = outline, color = setup), 
  #         stat = "identity", fill = NA, size = 0.3, color = "black", width = 0.75) +
  geom_bar(data = subset(ggbar.df, compare == "more"), aes(y = outline, color = setup), 
           stat = "identity", fill = NA, size = 0.3, color = c("#999933","#0077BB","#AA4499","#33BBEE"), width = 0.75) + # Coloured version
  geom_text(data = subset(ggbar.df, compare == "less"), aes(x = setup, y = pos.t.FP.sp.perc, label = paste0(t.FP.sp.perc,"%")),
            colour="white", family="Calibri", size=4.5, angle = 270) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        ###
        plot.margin = margin(0.5, 0.5, 0.6, 0.5, unit = "cm"),
        plot.title    = element_text(size=22, family = "Calibri", face = "bold", hjust = 0.5, vjust = 1),
        plot.subtitle = element_text(size=20, family = "Calibri"),
        axis.title.x  = element_text(size=17, family = "Calibri", vjust=-0.3),
        axis.title.y  = element_text(size=17, family = "Calibri", vjust=2.5),
        axis.text.x   = element_text(size=12.5, family = "Calibri"),
        axis.text.y   = element_blank(),
        ###
        axis.ticks = element_blank(), 
        legend.key = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank(),
        legend.position = "none") +
  labs(x="TSS threshold", y=paste0("Percentage of species with false
positive prediction less than ",FPthresh," cells"), fill="xyz") +
  scale_fill_manual(values = c("#999933","#0077BB","#AA4499","#33BBEE")) +
  scale_color_manual(values = c("#999933","#0077BB","#AA4499","#33BBEE")) +
  scale_y_continuous(position = "right", expand=c(0.005,0))

ggsave("EMbar_FP.png", width = 5.65, height = 11.7, units = "cm", dpi = 800)


### FN
ggplot(ggbar.df, aes(x=setup, y=t.FN.sp.perc)) +
  geom_bar(data = subset(ggbar.df, compare == "less"), aes(fill = setup), 
           stat = "identity", position = "stack", width = 0.75) +
  #geom_bar(data = subset(ggbar.df, compare == "more"), aes(y = outline, color = setup), 
  #         stat = "identity", fill = NA, size = 0.3, color = "black", width = 0.75) +
  geom_bar(data = subset(ggbar.df, compare == "more"), aes(y = outline, color = setup), 
           stat = "identity", fill = NA, size = 0.3, color = c("#999933","#0077BB","#AA4499","#33BBEE"), width = 0.75) + # Coloured version
  geom_text(data = subset(ggbar.df, compare == "less"), aes(x = setup, y = pos.t.FN.sp.perc, label = paste0(t.FN.sp.perc,"%")),
            colour="white", family="Calibri", size=4.5, angle = 270) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        ###
        plot.margin = margin(0.5, 0.5, 0.6, 0.5, unit = "cm"),
        plot.title    = element_text(size=22, family = "Calibri", face = "bold", hjust = 0.5, vjust = 1),
        plot.subtitle = element_text(size=20, family = "Calibri"),
        axis.title.x  = element_text(size=17, family = "Calibri", vjust=-0.3),
        axis.title.y  = element_text(size=17, family = "Calibri", vjust=2.5),
        axis.text.x   = element_text(size=12.5, family = "Calibri"),
        axis.text.y   = element_blank(),
        ###
        axis.ticks = element_blank(), 
        legend.key = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank(),
        legend.position = "none") +
  labs(x="TSS threshold", y=paste0("Percentage of species with false
negative prediction less than ",FNthresh," cells"), fill="xyz") +
  scale_fill_manual(values = c("#999933","#0077BB","#AA4499","#33BBEE")) +
  scale_color_manual(values = c("#999933","#0077BB","#AA4499","#33BBEE")) +
  scale_y_continuous(position = "right", expand=c(0.005,0))

ggsave("EMbar_FN.png", width = 5.65, height = 11.7, units = "cm", dpi = 800)

