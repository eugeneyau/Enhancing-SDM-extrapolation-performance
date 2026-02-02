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

#save.filename <- "EM_bmod424_PCvspEMmean_trunc_rev_SOUTH_trun20perc" # All untuned algo
#save.filename.selectedalgo <- "EM_bmod424_PCvspEMmean_trunc_rev_best4algo2_SOUTH_trun20perc" # Selected untuned algo

save.filename.selectedalgo <- "NA" # Selected untuned algo

save.filename <- "EM_bmod424_trunc_rev_SOUTH_trun20perc" # All untuned algo
save.filename.selectedalgo <- "EM_bmod424_trunc_rev_bestalgo_SOUTH_trun20perc" # Selected untuned algo


diff.perc.df <- readRDS(paste0(save.filename,".RData"))


########## Format data for ggplot ##########

# Define the list of items
items <- c("EMca", "EMmean", "EMwmean", "EMmedian")

# Create an empty data frame to store the results
ggplot_columns <- c("Species","setup",
                    "FP","FP.perc","FN","FN.perc","F","F.perc")
ggplot.df <- data.frame(matrix(nrow = 0, ncol = length(ggplot_columns)))
colnames(ggplot.df) <- ggplot_columns

# Analysis function
Accuracy_Fun <- function(i, diff.perc.df) { # Pass diff.perc.df as an argument
  sp.data <- diff.perc.df[i, ]
  species <- sp.data$Species
  
  for (item in items) {
    # Construct the column name prefixes
    prefix <- paste0(item, ".") # Handle "All" case
    
    # Extract data based on the current item
    t.FP <- sp.data[[paste0(prefix, "FP")]]
    t.FP.perc <- sp.data[[paste0(prefix, "FP.perc")]]
    t.FN <- sp.data[[paste0(prefix, "FN")]]
    t.FN.perc <- sp.data[[paste0(prefix, "FN.perc")]]
    t.F <- sp.data[[paste0(prefix, "F")]]
    t.F.perc <- sp.data[[paste0(prefix, "F.perc")]]
    
    # Create the data list
    item.data <- list(species, item, t.FP, t.FP.perc, t.FN, t.FN.perc, t.F, t.F.perc)
    
    # Add the data to the ggplot data frame
    ggplot.df[nrow(ggplot.df) + 1,] <<- item.data
  }
}

# Apply the function to each virtual species
lapply(1:nrow(diff.perc.df), Accuracy_Fun, diff.perc.df = diff.perc.df) 


if (!save.filename.selectedalgo=="NA") {
  full.ggplot.df <- ggplot.df
  
  diff.perc.df <- readRDS(paste0(save.filename.selectedalgo,".RData"))
  
  ggplot.df <- data.frame(matrix(nrow = 0, ncol = length(ggplot_columns)))
  colnames(ggplot.df) <- ggplot_columns
  
  lapply(1:nrow(diff.perc.df), Accuracy_Fun, diff.perc.df = diff.perc.df) 
}


# Write bar plot DF
columns <- c("Item","Metric","Mean","SD","Median")
barplot.df <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(barplot.df) <- columns
barplot.perc.df <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(barplot.perc.df) <- columns

# Loop through the items and create subsets and calculate statistics
for (item in items) {
  # Create the subset name
  subset_name <- paste0("data.", item)
  
  # Create the subset 
  assign(subset_name, subset(ggplot.df, setup == item))
  
  # Correct item name
  item_name <- item

  # Calculate statistics and add to barplot data frames
  barplot.df[nrow(barplot.df) + 1,] <- list(item_name, "tFP", mean(get(subset_name)$FP), sd(get(subset_name)$FP), median(get(subset_name)$FP))
  barplot.df[nrow(barplot.df) + 1,] <- list(item_name, "tFN", mean(get(subset_name)$FN), sd(get(subset_name)$FN), median(get(subset_name)$FN))
  barplot.perc.df[nrow(barplot.perc.df) + 1,] <- list(item_name, "tFP.perc", mean(get(subset_name)$FP.perc), sd(get(subset_name)$FP.perc), median(get(subset_name)$FP.perc))
  barplot.perc.df[nrow(barplot.perc.df) + 1,] <- list(item_name, "tFN.perc", mean(get(subset_name)$FN.perc), sd(get(subset_name)$FN.perc), median(get(subset_name)$FN.perc))
}

if (!save.filename.selectedalgo=="NA") {
  # Loop through the items and create subsets and calculate statistics
  for (item in items) {
    # Create the subset name
    subset_name <- paste0("data.", item)
    
    # Create the subset 
    assign(subset_name, subset(full.ggplot.df, setup == item))
    
    # Correct item name
    item_name <- paste0(item,"_all_untuned")
    
    # Calculate statistics and add to barplot data frames
    barplot.df[nrow(barplot.df) + 1,] <- list(item_name, "tFP", mean(get(subset_name)$FP), sd(get(subset_name)$FP), median(get(subset_name)$FP))
    barplot.df[nrow(barplot.df) + 1,] <- list(item_name, "tFN", mean(get(subset_name)$FN), sd(get(subset_name)$FN), median(get(subset_name)$FN))
    barplot.perc.df[nrow(barplot.perc.df) + 1,] <- list(item_name, "tFP.perc", mean(get(subset_name)$FP.perc), sd(get(subset_name)$FP.perc), median(get(subset_name)$FP.perc))
    barplot.perc.df[nrow(barplot.perc.df) + 1,] <- list(item_name, "tFN.perc", mean(get(subset_name)$FN.perc), sd(get(subset_name)$FN.perc), median(get(subset_name)$FN.perc))
  }
}

barplot.df$Metric <- factor(barplot.df$Metric, levels = c("tFP", "tFN"))
barplot.perc.df$Metric <- factor(barplot.perc.df$Metric, levels = c("tFP.perc", "tFN.perc"))

replacement_vector <- c(
  "EMca" = "Committee
averaging",
  "EMmean" = "Mean",
  "EMwmean" = "Weighted
mean",
  "EMmedian" = "Median")

if (!save.filename.selectedalgo=="NA") {
  replacement_vector <- c(
  "EMca" = "Committee
averaging
(selected algo)",
  "EMmean" = "Mean
(selected algo)",
  "EMwmean" = "Weighted
mean
(selected algo)",
  "EMmedian" = "Median
(selected algo)",
  "EMca_all_untuned" = "Committee
averaging",
  "EMmean_all_untuned" = "Mean",
  "EMwmean_all_untuned" = "Weighted
mean",
  "EMmedian_all_untuned" = "Median")
}

barplot.df$Item <- replacement_vector[barplot.df$Item]


# Compare best algo (stacked)
ggplot(data = barplot.df) + 
  geom_bar(aes(x = Item, y = Mean, fill = Metric),stat = "identity") +
  geom_hline(yintercept = 0, color = "dark grey", linewidth=0.5) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth=0.7),
        ###
        plot.margin = margin(0.5, 0.5, 0, 0.5, unit = "cm"),
        plot.title    = element_text(size=22, family = "Calibri", face = "bold", hjust = 0.5, vjust = 1),
        plot.subtitle = element_text(size=22, family = "Calibri"),
        axis.title.x  = element_text(size=22, family = "Calibri", vjust=-1.5),
        axis.title.y  = element_text(size=22, family = "Calibri", vjust=2),
        axis.text.x   = element_text(size=18, family = "Calibri"),
        axis.text.y   = element_text(size=18, family = "Calibri"),
        ###
        axis.ticks = element_line(linewidth = 0.6) , 
        axis.ticks.length = unit(0.2, "cm"),
        ###
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.background = element_blank(),
        legend.text = element_text(size = 18, family = "Calibri"),
        legend.title = element_text(size = 18, family = "Calibri", margin = margin(b=25))) +
  labs(x = "SDM ensemble algorithm", y = "Mean false predictions", fill=" ") +
  scale_fill_manual(values=c("#33BBEE","#882255"), 
                    labels=c("False positive      ","False negative      "))

ggsave("EM_best_FPFNbar_count_mean_staked.png", width = 28, height = 22, units = "cm", dpi = 800)

