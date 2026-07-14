library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(PMCMRplus) 

wd <- "/lustre1/g/sbs_bonebrake/Eugene/Paper1_Vis"
wdpc <- "C:/Users/Eugene/Desktop/Paper1_Vis_VIF5"
wdG15 <- "C:/Users/skywa/Desktop/Paper1_Vis_VIF5"

setwd(wdpc)

options(scipen = 999)

save.filename.EMca <- "BINMETH_VIF5_TruncS20_EMca"
save.filename.EMmean <- "BINMETH_VIF5_TruncS20_EMmean"

diff.perc.df.EMca <- readRDS(paste0(save.filename.EMca,".RData"))
diff.perc.df.EMmean <- readRDS(paste0(save.filename.EMmean,".RData"))

colnames(diff.perc.df.EMca)[colnames(diff.perc.df.EMca) != "Species"] <- paste0("EMca.", colnames(diff.perc.df.EMca)[colnames(diff.perc.df.EMca) != "Species"])
colnames(diff.perc.df.EMmean)[colnames(diff.perc.df.EMmean) != "Species"] <- paste0("EMmean.", colnames(diff.perc.df.EMmean)[colnames(diff.perc.df.EMmean) != "Species"])

diff.perc.df <- inner_join(diff.perc.df.EMca, diff.perc.df.EMmean, by = "Species")

bin_methods <- c("EMca.TSS", "EMca.ROC", "EMca.POD", "EMca.POFD",
               "EMmean.TSS", "EMmean.ROC", "EMmean.POD", "EMmean.POFD")

extracted_list <- list()

for (bin in bin_methods) {
  
  FP_col <- paste0(bin, ".FP")
  FN_col <- paste0(bin, ".FN")
  
  if (FP_col %in% colnames(diff.perc.df) && FN_col %in% colnames(diff.perc.df)) {
    
    bin_display <- case_when(
      bin == "EMca.TSS" ~ "TSS\n(Committee\naveraging)",
      bin == "EMca.ROC" ~ "ROC\n(Committee\naveraging)",
      bin == "EMca.POD" ~ "POD\n(Committee\naveraging)",
      bin == "EMca.POFD" ~ "POFD\n(Committee\naveraging)",
      bin == "EMmean.TSS" ~ "TSS\n(Mean)",
      bin == "EMmean.ROC" ~ "ROC\n(Mean)",
      bin == "EMmean.POD" ~ "POD\n(Mean)",
      bin == "EMmean.POFD" ~ "POFD\n(Mean)",
      TRUE ~ bin
    )
    
    fp_vector <- as.vector(as.numeric(unlist(diff.perc.df[[FP_col]])))
    fn_vector <- as.vector(as.numeric(unlist(diff.perc.df[[FN_col]])))
    
    tmp <- data.frame(
      Species = as.character(diff.perc.df$Species),
      BIN = bin_display,
      FP = fp_vector,
      FN = fn_vector,
      stringsAsFactors = FALSE
    )
    
    tmp$F <- tmp$FP + tmp$FN
    
    extracted_list[[bin]] <- tmp
  }
}

analysis_df <- do.call(rbind, extracted_list)
analysis_df$BIN <- as.factor(analysis_df$BIN)
analysis_df$Species <- as.factor(analysis_df$Species)

target_metrics <- c("F", "FP", "FN")
titles <- c("(a) Overall performance (all false predictions)",
            "(b) Overestimation (false positive)            ",
            "(c) Underestimation (false negative)           ")
filenames <- c("F", "FP", "FN")

# Create empty list to store the figures for ggarrange
fig_list <- list()

for (m in 1:3) {
  metric_name <- target_metrics[m]
  form <- as.formula(paste(metric_name, "~ BIN | Species"))
  
  friedman_res <- friedman.test(form, data = analysis_df)
  cat("\n=====================================\n", titles[m], "\nGlobal Friedman p-value:", friedman_res$p.value, "\n")
  
  nemenyi_res <- frdAllPairsNemenyiTest(form, data = analysis_df)
  p_matrix <- nemenyi_res$p.value
  
  matrix_df <- as.data.frame(p_matrix) %>%
    mutate(BIN1 = rownames(p_matrix)) %>%
    pivot_longer(-BIN1, names_to = "BIN2", values_to = "p_val") %>%
    filter(!is.na(p_val)) 
  
  matrix_df <- matrix_df %>%
    mutate(
      ColorKey = case_when(
        p_val < 0.001 ~ "p001",
        p_val < 0.01  ~ "p01",
        p_val < 0.05  ~ "p05",
        TRUE          ~ "ns"
      ),
      ColorKey = factor(ColorKey, levels = c("ns", "p05", "p01", "p001")),
      Label = case_when(
        p_val < 0.001 ~ "< 0.001",
        p_val < 0.01  ~ "< 0.01",
        p_val < 0.05  ~ "< 0.05",
        TRUE          ~ "N.S."
      ),
      LabelValue = case_when(
        p_val < 0.05 ~ formatC(p_val, digits = 4, format = "f"),
        TRUE         ~ "N.S."
      ),
      LabelAllValue = formatC(p_val, digits = 4, format = "f")
    )
  
  fig <- 
    ggplot(matrix_df, aes(x = BIN2, y = BIN1, fill = ColorKey)) +
    geom_tile(color = "white", linewidth = 0.5, show.legend = TRUE) +
    geom_text(aes(label = LabelAllValue), color = "black", fontface = "plain", size = 4.5, family = "Calibri") +
    scale_fill_manual(
      values = c("grey85", "#6EA6CD", "#4A7BB7", "#364B9A"),
      labels = c("Non-significant    ", "p < 0.05    ", "p < 0.01    ", "p < 0.001    "),
      drop = FALSE
    ) +
    scale_x_discrete(position = "bottom") +
    scale_y_discrete(limits = rev) +
    labs(subtitle = titles[m], x = "Binary method", y = "Binary method", fill = " ") + 
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.margin = margin(0.3, 0.3, 0.3, 0.3, unit = "cm"),
      #plot.title.position = "plot", #Aligns titles to plot left edge
      plot.subtitle = element_text(size = 18, family = "Calibri", hjust = 0.025), 
      axis.title.x  = element_blank(),
      axis.title.y  = element_blank(),
      axis.text.x   = element_text(size = 12.5, family = "Calibri"),
      axis.text.y   = element_text(size = 12.5, family = "Calibri"),
      axis.ticks = element_blank(), 
      legend.key = element_blank(),
      legend.text = element_blank(),
      legend.title = element_blank(),
      legend.position = "none"
    ) +
    guides(fill = guide_legend(
      keywidth = unit(1, "cm"),
      keyheight = unit(0.6, "cm"))
    )
  
  #print(fig)
  #ggsave(paste0(save.filename,"_",filenames[m],".png"), plot = fig, width = 30, height = 28.5, units = "cm", dpi = 800)
  fig_list[[m]] <- fig
  
}

combined_fig <- ggarrange(fig_list[[1]], fig_list[[2]], fig_list[[3]],
                          nrow = 3, ncol = 1, align = "v"#, 
                          #common.legend = TRUE, legend = "bottom"
                          )

print(combined_fig)
ggsave("BINmeth_VIF5_sig_combined.png", plot = combined_fig, width = 20, height = 35, units = "cm", dpi = 800)




### Plot legend
legend <- ggplot(matrix_df, aes(x = Algo2, y = Algo1, fill = ColorKey)) +
  geom_tile(color = "white", linewidth = 0.5, show.legend = TRUE) +
  geom_text(aes(label = LabelAllValue), color = "black", fontface = "plain", size = 3.5, family = "Calibri") +
  #geom_text(aes(label = LabelValue), color = "black", fontface = "plain", size = 3.5, family = "Calibri") +
  #geom_text(aes(label = Label), color = "black", fontface = "plain", size = 3.5, family = "Calibri") +
  scale_fill_manual(
    values = c("grey85", "#6EA6CD", "#4A7BB7", "#364B9A"),
    labels = c("Non-significant    ", "p < 0.05    ", "p < 0.01    ", "p < 0.001    "),
    #labels = c("Non-significant (N.S.)    ", "p < 0.05    ", "p < 0.01    ", "p < 0.001    "),
    drop = FALSE
  ) +
  scale_x_discrete(position = "bottom") +
  scale_y_discrete(limits = rev) +
  labs(title = titles[m], x = "SDM algorithm", y = "SDM algorithm", fill = " ") +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.margin = margin(0.5, 0.5, 0, 0.5, unit = "cm"),
    plot.title    = element_text(size = 22, family = "Calibri", face = "bold", hjust = 0.5, vjust = 1),
    plot.subtitle = element_text(size = 20, family = "Calibri"),
    axis.title.x  = element_blank(),
    axis.title.y  = element_blank(),
    axis.text.x   = element_text(size = 13, family = "Calibri", angle = 45, hjust = 1),
    axis.text.y   = element_text(size = 13, family = "Calibri"),
    axis.ticks = element_blank(), 
    legend.direction = "horizontal",
    legend.position = "bottom",
    legend.background = element_blank(),
    legend.box.spacing = unit(0.5, "cm"),
    legend.text = element_text(size = 15, family = "Calibri", vjust = 0.6),
    legend.title = element_text(size = 15, family = "Calibri", margin = margin(b = 25))
  ) +
  guides(fill = guide_legend(
    keywidth = unit(1, "cm"),
    keyheight = unit(0.6, "cm"))
  )

ggsave("Sig_matrix_legend.png", plot = legend, width = 25, height = 18, units = "cm", dpi = 800)

