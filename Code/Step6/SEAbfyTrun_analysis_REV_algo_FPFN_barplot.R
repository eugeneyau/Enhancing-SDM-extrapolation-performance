extrafont::loadfonts(device="win")
extrafont::loadfonts(device="postscript")
library(ggplot2)
library(ggpubr)
library(dplyr)

wd <- "/lustre1/g/sbs_bonebrake/Eugene/Paper1_Vis"  # HPC
wdpc <- "C:/Users/Eugene/Desktop/Paper1_Vis"  # PC
wdG15 <- "C:/Users/skywa/Desktop/Paper1_Vis"  # Laptop

setwd(wdG15)

options(scipen = 999)

save.filename <- "Algo_REV_EMca_proj_TSSNA" # Name the df output file
#save.filename <- "Algo_REVPC_EMca_proj_TSSNA" # Name the df output file

diff.perc.df <- readRDS(paste0(save.filename,".RData"))


########## Format data for ggplot ##########

# Define the list of algorithms
algorithms <- c("All", "ANN", "CTA", "GBM", "GLM", "MARS", "MAXNET", "RF", "XGBOOST",
                "tuned_ANN", "tuned_CTA", "tuned_GBM", "tuned_GLM",
                "tuned_MARS", "tuned_RF")


# Create an empty data frame to store the results
ggplot_columns <- c("Species","setup",
                    "f.FP","f.FP.perc","f.FN","f.FN.perc","f.F","f.F.perc",
                    "t.FP","t.FP.perc","t.FN","t.FN.perc","t.F","t.F.perc")
ggplot.df <- data.frame(matrix(nrow = 0, ncol = length(ggplot_columns)))
colnames(ggplot.df) <- ggplot_columns

# Analysis function
Accuracy_Fun <- function(i, diff.perc.df) { # Pass diff.perc.df as an argument
  sp.data <- diff.perc.df[i, ]
  species <- sp.data$Species
  
  for (algorithm in algorithms) {
    # Construct the column name prefixes
    prefix <- paste0(algorithm, ".") # Handle "All" case
    
    # Extract data based on the current algorithm
    f.FP <- sp.data[[paste0(prefix, "f.FP")]]
    f.FP.perc <- sp.data[[paste0(prefix, "f.FP.perc")]]
    f.FN <- sp.data[[paste0(prefix, "f.FN")]]
    f.FN.perc <- sp.data[[paste0(prefix, "f.FN.perc")]]
    f.F <- sp.data[[paste0(prefix, "f.F")]]
    f.F.perc <- sp.data[[paste0(prefix, "f.F.perc")]]
    t.FP <- sp.data[[paste0(prefix, "t.FP")]]
    t.FP.perc <- sp.data[[paste0(prefix, "t.FP.perc")]]
    t.FN <- sp.data[[paste0(prefix, "t.FN")]]
    t.FN.perc <- sp.data[[paste0(prefix, "t.FN.perc")]]
    t.F <- sp.data[[paste0(prefix, "t.F")]]
    t.F.perc <- sp.data[[paste0(prefix, "t.F.perc")]]
    
    # Create the data list
    algorithm.data <- list(species, algorithm, f.FP, f.FP.perc, f.FN, f.FN.perc, f.F, f.F.perc,
                           t.FP, t.FP.perc, t.FN, t.FN.perc, t.F, t.F.perc)
    
    # Add the data to the ggplot data frame
    ggplot.df[nrow(ggplot.df) + 1,] <<- algorithm.data
  }
}

# Apply the function to each virtual species
lapply(1:nrow(diff.perc.df), Accuracy_Fun, diff.perc.df = diff.perc.df) # Pass diff.perc.df

# Write bar plot DF
columns <- c("Algo","Metric","Mean","SD","Median")

barplot.full.df <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(barplot.full.df) <- columns
barplot.full.perc.df <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(barplot.full.perc.df) <- columns

barplot.df <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(barplot.df) <- columns
barplot.perc.df <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(barplot.perc.df) <- columns

# Loop through the algorithms and create subsets and calculate statistics
for (algorithm in algorithms) {
  # Create the subset name
  subset_name <- paste0("data.", algorithm)
  
  # Create the subset (handling "All" case)
  if (algorithm == "All") {
    assign(subset_name, subset(ggplot.df, setup == "All"))
  } else {
    assign(subset_name, subset(ggplot.df, setup == algorithm))
  }
  
  # Correct algo name
  if (algorithm == "ANN"||algorithm == "CTA"||algorithm == "GBM"||algorithm == "GLM"||algorithm == "MARS"||algorithm == "RF") {
    algo_name <- algorithm
    } else if (algorithm == "All") {algo_name <- "ALL algo
untuned"} else if (algorithm == "MAXNET") {algo_name <- "MaxEnt"
} else if (algorithm == "XGBOOST") {algo_name <- "XGBoost"
} else if (algorithm == "tuned_ANN") {algo_name <- "ANN
tuned"} else if (algorithm == "tuned_CTA") {algo_name <- "CTA
tuned"} else if (algorithm == "tuned_GBM") {algo_name <- "GBM
tuned"} else if (algorithm == "tuned_GLM") {algo_name <- "GLM
tuned"} else if (algorithm == "tuned_MARS") {algo_name <- "MARS
tuned"} else if (algorithm == "tuned_RF") {algo_name <- "RF
tuned"}
  
  # Calculate statistics and add to barplot data frames
  barplot.full.df[nrow(barplot.full.df) + 1,] <- list(algo_name, "fFP", mean(get(subset_name)$f.FP), sd(get(subset_name)$f.FP), median(get(subset_name)$f.FP))
  barplot.full.df[nrow(barplot.full.df) + 1,] <- list(algo_name, "fFN", mean(get(subset_name)$f.FN), sd(get(subset_name)$f.FN), median(get(subset_name)$f.FN))
  barplot.full.perc.df[nrow(barplot.full.perc.df) + 1,] <- list(algo_name, "fFP.perc", mean(get(subset_name)$f.FP.perc), sd(get(subset_name)$f.FP.perc), median(get(subset_name)$f.FP.perc))
  barplot.full.perc.df[nrow(barplot.full.perc.df) + 1,] <- list(algo_name, "fFN.perc", mean(get(subset_name)$f.FN.perc), sd(get(subset_name)$f.FN.perc), median(get(subset_name)$f.FN.perc))
  
  barplot.df[nrow(barplot.df) + 1,] <- list(algo_name, "tFP", mean(get(subset_name)$t.FP), sd(get(subset_name)$t.FP), median(get(subset_name)$t.FP))
  barplot.df[nrow(barplot.df) + 1,] <- list(algo_name, "tFN", mean(get(subset_name)$t.FN), sd(get(subset_name)$t.FN), median(get(subset_name)$t.FN))
  barplot.perc.df[nrow(barplot.perc.df) + 1,] <- list(algo_name, "tFP.perc", mean(get(subset_name)$t.FP.perc), sd(get(subset_name)$t.FP.perc), median(get(subset_name)$t.FP.perc))
  barplot.perc.df[nrow(barplot.perc.df) + 1,] <- list(algo_name, "tFN.perc", mean(get(subset_name)$t.FN.perc), sd(get(subset_name)$t.FN.perc), median(get(subset_name)$t.FN.perc))
}

barplot.full.df$Metric <- factor(barplot.full.df$Metric, levels = c("fFP", "fFN"))
barplot.full.perc.df$Metric <- factor(barplot.full.perc.df$Metric, levels = c("fFP.perc", "fFN.perc"))

barplot.df$Metric <- factor(barplot.df$Metric, levels = c("tFP", "tFN"))
barplot.perc.df$Metric <- factor(barplot.perc.df$Metric, levels = c("tFP.perc", "tFN.perc"))


### Active: Mean dodged (truncated data) ###
ggplot(data = barplot.df) + 
  geom_bar(aes(x = Algo, y = Mean, fill = Metric),stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0, color = "dark grey", linewidth=0.5) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth=0.7),
        ###
        plot.margin = margin(0.5, 0.5, 0, 0.5, unit = "cm"),
        plot.title    = element_text(size=22, family = "Calibri", face = "bold", hjust = 0.5, vjust = 1),
        plot.subtitle = element_text(size=20, family = "Calibri"),
        axis.title.x  = element_text(size=19, family = "Calibri", vjust=-1.5),
        axis.title.y  = element_text(size=19, family = "Calibri", vjust=2),
        axis.text.x   = element_text(size=15, family = "Calibri"),
        axis.text.y   = element_text(size=15, family = "Calibri"),
        ###
        axis.ticks = element_line(linewidth = 0.6) , 
        axis.ticks.length = unit(0.2, "cm"),
        ###
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.background = element_blank(),
        legend.text = element_text(size = 15, family = "Calibri"),
        legend.title = element_text(size = 15, family = "Calibri", margin = margin(b=25))) +
  labs(x = "SDM algorithm", y = "Mean false predictions", fill=" ") +
  scale_fill_manual(values=c("#33BBEE","#882255"), 
                    labels=c("False positive      ","False negative      ")) #+
#scale_y_continuous(breaks = seq(0, 250, 50))

ggsave("algo_FPFNbar_count_mean_dodged_EMca.png", width = 30, height = 18, units = "cm", dpi = 800)


### Active: Mean dodged (full data) ###
ggplot(data = barplot.full.df) + 
  geom_bar(aes(x = Algo, y = Mean, fill = Metric),stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0, color = "dark grey", linewidth=0.5) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth=0.7),
        ###
        plot.margin = margin(0.5, 0.5, 0, 0.5, unit = "cm"),
        plot.title    = element_text(size=22, family = "Calibri", face = "bold", hjust = 0.5, vjust = 1),
        plot.subtitle = element_text(size=20, family = "Calibri"),
        axis.title.x  = element_text(size=19, family = "Calibri", vjust=-1.5),
        axis.title.y  = element_text(size=19, family = "Calibri", vjust=2),
        axis.text.x   = element_text(size=15, family = "Calibri"),
        axis.text.y   = element_text(size=15, family = "Calibri"),
        ###
        axis.ticks = element_line(linewidth = 0.6) , 
        axis.ticks.length = unit(0.2, "cm"),
        ###
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.background = element_blank(),
        legend.text = element_text(size = 15, family = "Calibri"),
        legend.title = element_text(size = 15, family = "Calibri", margin = margin(b=25))) +
  labs(x = "SDM algorithm", y = "Mean false predictions", fill=" ") +
  scale_fill_manual(values=c("#33BBEE","#882255"), 
                    labels=c("False positive      ","False negative      ")) #+
#scale_y_continuous(breaks = seq(0, 250, 50))

ggsave("algo_FPFNbar_count_fullmean_dodged_EMca.png", width = 30, height = 18, units = "cm", dpi = 800)


### Active: Mean stacked ###
ggplot(data = barplot.df) + 
  geom_bar(aes(x = Algo, y = Mean, fill = Metric),stat = "identity", position = "stack") +
  geom_hline(yintercept = 0, color = "dark grey", linewidth=0.5) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth=0.7),
        ###
        plot.margin = margin(0.5, 0.5, 0, 0.5, unit = "cm"),
        plot.title    = element_text(size=22, family = "Calibri", face = "bold", hjust = 0.5, vjust = 1),
        plot.subtitle = element_text(size=20, family = "Calibri"),
        axis.title.x  = element_text(size=19, family = "Calibri", vjust=-1.5),
        axis.title.y  = element_text(size=19, family = "Calibri", vjust=2),
        axis.text.x   = element_text(size=15, family = "Calibri"),
        axis.text.y   = element_text(size=15, family = "Calibri"),
        ###
        axis.ticks = element_line(linewidth = 0.6) , 
        axis.ticks.length = unit(0.2, "cm"),
        ###
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.background = element_blank(),
        legend.text = element_text(size = 15, family = "Calibri"),
        legend.title = element_text(size = 15, family = "Calibri", margin = margin(b=25))) +
  labs(x = "SDM algorithm", y = "False predictions (mean)", fill=" ") +
  scale_fill_manual(values=c("#33BBEE","#882255"), 
                    labels=c("False positive      ","False negative      ")) #+
#scale_y_continuous(breaks = seq(0, 250, 50))

ggsave("algo_FPFNbar_count_mean_stacked_EMca.png", width = 28, height = 18, units = "cm", dpi = 800)



