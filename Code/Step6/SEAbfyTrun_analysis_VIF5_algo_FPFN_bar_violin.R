extrafont::loadfonts(device="win")
extrafont::loadfonts(device="postscript")
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(see)

wd <- "/lustre1/g/sbs_bonebrake/Eugene/Paper1_Vis_VIF5"  # HPC
wdpc <- "C:/Users/Eugene/Desktop/Paper1_Vis_VIF5"  # PC
wdG15 <- "C:/Users/skywa/Desktop/Paper1_Vis_VIF5"  # Laptop

setwd(wdG15)

options(scipen = 999)

save.filename <- "Algo_VIF5_EMca_proj_TSSNA_trunSOUTH" # Name the df output file

diff.perc.df <- readRDS(paste0(save.filename,".RData"))


########## Format data for ggplot ##########
totalcells <- 125339

# Define the list of algorithms
algorithms <- c("All", "ANN", "CTA", "GBM", "GLM", "MARS", "MAXNET", "RF", "XGBOOST",
                "tuned_ANN", "tuned_CTA", "tuned_GBM", "tuned_GLM",
                "tuned_MARS", "tuned_RF")

# Create an empty data frame to store the results
ggplot_columns <- c("Species","setup",
                    "f.FP","f.FP.siteperc","f.FP.rangeperc","f.FN","f.FN.siteperc","f.FN.rangeperc","f.F","f.F.siteperc","f.F.rangeperc",
                    "t.FP","t.FP.siteperc","t.FP.rangeperc","t.FN","t.FN.siteperc","t.FN.rangeperc","t.F","t.F.siteperc","t.F.rangeperc")
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
    f.FP.siteperc <- f.FP/totalcells*100
    f.FP.rangeperc <- sp.data[[paste0(prefix, "f.FP.perc")]]
    f.FN <- sp.data[[paste0(prefix, "f.FN")]]
    f.FN.siteperc <- f.FN/totalcells*100
    f.FN.rangeperc <- sp.data[[paste0(prefix, "f.FN.perc")]]
    f.F <- sp.data[[paste0(prefix, "f.F")]]
    f.F.siteperc <- f.F/totalcells*100
    f.F.rangeperc <- sp.data[[paste0(prefix, "f.F.perc")]]
    
    t.FP <- sp.data[[paste0(prefix, "t.FP")]]
    t.FP.siteperc <- t.FP/totalcells*100
    t.FP.rangeperc <- sp.data[[paste0(prefix, "t.FP.perc")]]
    t.FN <- sp.data[[paste0(prefix, "t.FN")]]
    t.FN.siteperc <- t.FN/totalcells*100
    t.FN.rangeperc <- sp.data[[paste0(prefix, "t.FN.perc")]]
    t.F <- sp.data[[paste0(prefix, "t.F")]]
    t.F.siteperc <- t.F/totalcells*100
    t.F.rangeperc <- sp.data[[paste0(prefix, "t.F.perc")]]
    
    # Create the data list
    algorithm.data <- list(species, algorithm, 
                           f.FP, f.FP.siteperc, f.FP.rangeperc, f.FN, f.FN.siteperc, f.FN.rangeperc, f.F, f.F.siteperc, f.F.rangeperc,
                           t.FP, t.FP.siteperc, t.FP.rangeperc, t.FN, t.FN.siteperc, t.FN.rangeperc, t.F, t.F.siteperc, t.F.rangeperc)
    
    # Add the data to the ggplot data frame
    ggplot.df[nrow(ggplot.df) + 1,] <<- algorithm.data
  }
}

# Apply the function to each virtual species
lapply(1:nrow(diff.perc.df), Accuracy_Fun, diff.perc.df = diff.perc.df) # Pass diff.perc.df

algo_levels <- c(
  "ALL algo\nuntuned", "ANN", "ANN\ntuned", "CTA", "CTA\ntuned", "GBM", "GBM\ntuned",
  "GLM", "GLM\ntuned", "MARS", "MARS\ntuned", "MaxEnt", "RF", "RF\ntuned", "XGBoost")

algo_clean_map <- c(
  "All" = "ALL algo\nuntuned", "ANN" = "ANN", "CTA" = "CTA", "GBM" = "GBM", "GLM" = "GLM", 
  "MARS" = "MARS", "MAXNET" = "MaxEnt", "RF" = "RF", "XGBOOST" = "XGBoost",
  "tuned_ANN" = "ANN\ntuned", "tuned_CTA" = "CTA\ntuned", "tuned_GBM" = "GBM\ntuned",
  "tuned_GLM" = "GLM\ntuned", "tuned_MARS" = "MARS\ntuned", "tuned_RF" = "RF\ntuned")

violin_raw.df <- ggplot.df %>%
  mutate(Algo = factor(algo_clean_map[setup], levels = algo_levels)) %>%
  select(Species, Algo, t.FP.siteperc, t.FN.siteperc) %>%
  pivot_longer(cols = c(t.FP.siteperc, t.FN.siteperc), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = c("t.FP.siteperc", "t.FN.siteperc")))

columns <- c("Algo","Metric","Mean","SD","Median")
barplot.full.df <- data.frame(matrix(nrow = 0, ncol = length(columns))); colnames(barplot.full.df) <- columns
barplot.full.siteperc.df <- data.frame(matrix(nrow = 0, ncol = length(columns))); colnames(barplot.full.siteperc.df) <- columns
barplot.full.rangeperc.df <- data.frame(matrix(nrow = 0, ncol = length(columns))); colnames(barplot.full.rangeperc.df) <- columns
barplot.df <- data.frame(matrix(nrow = 0, ncol = length(columns))); colnames(barplot.df) <- columns
barplot.siteperc.df <- data.frame(matrix(nrow = 0, ncol = length(columns))); colnames(barplot.siteperc.df) <- columns
barplot.rangeperc.df <- data.frame(matrix(nrow = 0, ncol = length(columns))); colnames(barplot.rangeperc.df) <- columns

for (algorithm in algorithms) {
  subset_name <- paste0("data.", algorithm)
  if (algorithm == "All") { assign(subset_name, subset(ggplot.df, setup == "All")) } else { assign(subset_name, subset(ggplot.df, setup == algorithm)) }
  if (algorithm == "ANN"||algorithm == "CTA"||algorithm == "GBM"||algorithm == "GLM"||algorithm == "MARS"||algorithm == "RF") { algo_name <- algorithm
  } else if (algorithm == "All") {algo_name <- "ALL algo\nuntuned"} else if (algorithm == "MAXNET") {algo_name <- "MaxEnt"
  } else if (algorithm == "XGBOOST") {algo_name <- "XGBoost"} else if (algorithm == "tuned_ANN") {algo_name <- "ANN\ntuned"
  } else if (algorithm == "tuned_CTA") {algo_name <- "CTA\ntuned"} else if (algorithm == "tuned_GBM") {algo_name <- "GBM\ntuned"
  } else if (algorithm == "tuned_GLM") {algo_name <- "GLM\ntuned"} else if (algorithm == "tuned_MARS") {algo_name <- "MARS\ntuned"
  } else if (algorithm == "tuned_RF") {algo_name <- "RF\ntuned"}
  
  barplot.full.df[nrow(barplot.full.df) + 1,] <- list(algo_name, "fFP", mean(get(subset_name)$f.FP), sd(get(subset_name)$f.FP), median(get(subset_name)$f.FP))
  barplot.full.df[nrow(barplot.full.df) + 1,] <- list(algo_name, "fFN", mean(get(subset_name)$f.FN), sd(get(subset_name)$f.FN), median(get(subset_name)$f.FN))
  barplot.full.siteperc.df[nrow(barplot.full.siteperc.df) + 1,] <- list(algo_name, "fFP.siteperc", mean(get(subset_name)$f.FP.siteperc), sd(get(subset_name)$f.FP.siteperc), median(get(subset_name)$f.FP.siteperc))
  barplot.full.siteperc.df[nrow(barplot.full.siteperc.df) + 1,] <- list(algo_name, "fFN.siteperc", mean(get(subset_name)$f.FN.siteperc), sd(get(subset_name)$f.FN.siteperc), median(get(subset_name)$f.FN.siteperc))
  barplot.full.rangeperc.df[nrow(barplot.full.rangeperc.df) + 1,] <- list(algo_name, "fFP.rangeperc", mean(get(subset_name)$f.FP.rangeperc), sd(get(subset_name)$f.FP.rangeperc), median(get(subset_name)$f.FP.rangeperc))
  barplot.full.rangeperc.df[nrow(barplot.full.rangeperc.df) + 1,] <- list(algo_name, "fFN.rangeperc", mean(get(subset_name)$f.FN.rangeperc), sd(get(subset_name)$f.FN.rangeperc), median(get(subset_name)$f.FN.rangeperc))
  barplot.df[nrow(barplot.df) + 1,] <- list(algo_name, "tFP", mean(get(subset_name)$t.FP), sd(get(subset_name)$t.FP), median(get(subset_name)$t.FP))
  barplot.df[nrow(barplot.df) + 1,] <- list(algo_name, "tFN", mean(get(subset_name)$t.FN), sd(get(subset_name)$t.FN), median(get(subset_name)$t.FN))
  barplot.siteperc.df[nrow(barplot.siteperc.df) + 1,] <- list(algo_name, "tFP.siteperc", mean(get(subset_name)$t.FP.siteperc), sd(get(subset_name)$t.FP.siteperc), median(get(subset_name)$t.FP.siteperc))
  barplot.siteperc.df[nrow(barplot.siteperc.df) + 1,] <- list(algo_name, "tFN.siteperc", mean(get(subset_name)$t.FN.siteperc), sd(get(subset_name)$t.FN.siteperc), median(get(subset_name)$t.FN.siteperc))
  barplot.rangeperc.df[nrow(barplot.rangeperc.df) + 1,] <- list(algo_name, "tFP.rangeperc", mean(get(subset_name)$t.FP.rangeperc), sd(get(subset_name)$t.FP.rangeperc), median(get(subset_name)$t.FP.rangeperc))
  barplot.rangeperc.df[nrow(barplot.rangeperc.df) + 1,] <- list(algo_name, "tFN.rangeperc", mean(get(subset_name)$t.FN.rangeperc), sd(get(subset_name)$t.FN.rangeperc), median(get(subset_name)$t.FN.rangeperc))
}

# Set violin fill color
violin_raw.df$FP_ViolinColor <- ifelse(violin_raw.df$Metric == "t.FP.siteperc", "#A9E5F7", "transparent")
violin_raw.df$FN_ViolinColor <- ifelse(violin_raw.df$Metric == "t.FN.siteperc", "#D996B7", "transparent")

# Set violin line color
violin_raw.df$FP_LineColor   <- ifelse(violin_raw.df$Metric == "t.FP.siteperc", "black", "transparent")
violin_raw.df$FN_LineColor   <- ifelse(violin_raw.df$Metric == "t.FN.siteperc", "black", "transparent")

# Set bar color
violin_raw.df$BarColor       <- ifelse(violin_raw.df$Metric == "t.FP.siteperc", "#33BBEE", "#882255")

################################################################################
################################################################################

### Format legend item
draw_key_custom_box <- function(data, params, size) {
  fill_val  <- if (!is.null(data$fill)) data$fill else "white"
  alpha_val <- if (!is.null(data$alpha)) data$alpha else 1
  col_val   <- if (!is.null(data$colour)) data$colour else "black"
  lty_val   <- if (!is.null(data$linetype)) data$linetype else 1
  lwd_val   <- if (!is.null(data$linewidth)) data$linewidth else 0.5
  fatten_val <- if (!is.null(params$fatten)) params$fatten else 2.5
  
  # Create box
  box_grob <- rectGrob(
    x = unit(0.5, "npc"), y = unit(0.5, "npc"),
    width = unit(0.5, "npc"), height = unit(0.6, "npc"),
    gp = gpar(
      fill = alpha(fill_val, alpha_val),
      col = col_val,
      lty = lty_val,
      lwd = lwd_val * .pt))
  
  # Create median line
  median_grob <- linesGrob(
    x = unit(c(0.25, 0.75), "npc"),
    y = unit(c(0.5, 0.5), "npc"),
    gp = gpar(
      col = col_val,
      lwd = lwd_val * fatten_val * .pt,
      lineend = "butt"))
  
  grobTree(box_grob, median_grob)
}

### Actual ggplot
ggplot(data = violin_raw.df, aes(x = Algo, y = Value)) + 
  
  # Violin FP (left)
  geom_violinhalf(aes(fill = FP_ViolinColor, color = FP_LineColor), flip = TRUE, 
                  linewidth = 0.4, position = position_nudge(x = 0), width = 1.5) +
  
  # Violin FN (right)
  geom_violinhalf(aes(fill = FN_ViolinColor, color = FN_LineColor), linewidth = 0.4, 
                  position = position_nudge(x = 0), width = 1.5) +
  
  # Boxplot
  stat_summary(aes(fill = BarColor, group = interaction(Algo, BarColor)),
               fun.data = function(x) {
                 stats <- c(
                   ymin   = as.numeric(quantile(x, 0.25)), # Hide whisker
                   lower  = as.numeric(quantile(x, 0.25)),
                   middle = as.numeric(median(x)),
                   upper  = as.numeric(quantile(x, 0.75)),
                   ymax   = as.numeric(quantile(x, 0.75))) # Hide whisker
                 return(stats)},
               geom = "boxplot", position = position_dodge(0.3), width = 0.3,
               color = "black", linewidth = 0.45, fatten = 2.5, # Bold line for median
               key_glyph = draw_key_custom_box) + # Specify custom legend format
  
  # Log transform Y-axis
  scale_y_continuous(trans = "log1p", breaks = c(0, 1, 5, 25, 100)) +
  
  # Set color
  scale_fill_identity(name = " ", breaks = c("#33BBEE", "#882255"), guide = "legend",
                      labels = c("False positive      ", "False negative      ")) +
  scale_color_identity(guide = "none") +
  
  # Line at Y=0
  geom_hline(yintercept = 0, color = "dark grey", linewidth = 0.5) +
  
  # Set theme
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
        ###
        plot.margin = margin(0.5, 0.5, 0, 0.5, unit = "cm"),
        plot.title    = element_text(size = 22, family = "Calibri", face = "bold", hjust = 0.5, vjust = 1),
        plot.subtitle = element_text(size = 20, family = "Calibri"),
        axis.title.x  = element_text(size = 19, family = "Calibri", vjust = -1.5),
        axis.title.y  = element_text(size = 19, family = "Calibri", vjust = 2),
        axis.text.x   = element_text(size = 15, family = "Calibri"), 
        axis.text.y   = element_text(size = 15, family = "Calibri"),
        ###
        axis.ticks = element_line(linewidth = 0.6) , 
        axis.ticks.length = unit(0.2, "cm"),
        ###
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.background = element_blank(),
        legend.box.spacing = unit(0.5, "cm"),
        legend.text = element_text(size = 19, family = "Calibri"),
        legend.title = element_text(size = 19, family = "Calibri", margin = margin(b = 25)),
        legend.key.width  = unit(1.1, "cm"), 
        legend.key.height = unit(1.5, "cm")) +
  labs(x = "SDM algorithm", y = "False prediction percentage (log1p)")

ggsave("algo_violin_box_EMca.png", width = 33, height = 18.5, units = "cm", dpi = 800)

















ggplot(data = violin_raw.df, aes(x = Algo, y = Value)) + 
  
  # Violin FP (left)
  geom_violinhalf(aes(fill = FP_ViolinColor, color = FP_LineColor), flip = TRUE, 
                  linewidth = 0.4, position = position_nudge(x = 0), width = 1.5) +
  
  # Violin FN (right)
  geom_violinhalf(aes(fill = FN_ViolinColor, color = FN_LineColor), linewidth = 0.4, 
                  position = position_nudge(x = 0), width = 1.5) +
  
  # Boxplot
  stat_summary(aes(fill = BarColor, group = interaction(Algo, BarColor)),
               fun.data = function(x) {
                 stats <- c(
                   ymin   = as.numeric(quantile(x, 0.25)), # Hide whisker
                   lower  = as.numeric(quantile(x, 0.25)),
                   middle = as.numeric(median(x)),
                   upper  = as.numeric(quantile(x, 0.75)),
                   ymax   = as.numeric(quantile(x, 0.75))) # Hide whisker
                 return(stats)},
               geom = "boxplot", position = position_dodge(0.3), width = 0.3,
               color = "black", linewidth = 0.45, fatten = 2.5, # Bold line for median
               key_glyph = draw_key_rect) + 
  
  # Log transform Y-axis
  scale_y_continuous(trans = "log1p", breaks = c(0, 1, 5, 25, 100)) +
  
  # Set color
  scale_fill_identity(name = " ", breaks = c("#33BBEE", "#882255"), guide = "legend",
                      labels = c("False positive      ", "False negative      ")) +
  scale_color_identity(guide = "none") +
  
  # Line at Y=0
  geom_hline(yintercept = 0, color = "dark grey", linewidth = 0.5) +
  
  # Set theme
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
        ###
        plot.margin = margin(0.5, 0.5, 0, 0.5, unit = "cm"),
        plot.title    = element_text(size = 22, family = "Calibri", face = "bold", hjust = 0.5, vjust = 1),
        plot.subtitle = element_text(size = 20, family = "Calibri"),
        axis.title.x  = element_text(size = 19, family = "Calibri", vjust = -1.5),
        axis.title.y  = element_text(size = 19, family = "Calibri", vjust = 2),
        axis.text.x   = element_text(size = 15, family = "Calibri"), 
        axis.text.y   = element_text(size = 15, family = "Calibri"),
        ###
        axis.ticks = element_line(linewidth = 0.6) , 
        axis.ticks.length = unit(0.2, "cm"),
        ###
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.background = element_blank(),
        legend.box.spacing = unit(0.5, "cm"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15, margin = margin(b = 25))) +
  labs(x = "SDM algorithm", y = "False prediction percentage (log1p)")

ggsave("algo_violin_box_EMca.png", width = 33, height = 18.5, units = "cm", dpi = 800)

