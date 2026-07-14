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

save.filename <- "EM_TSSNA_VIF5_TruncS20"
#save.filename <- "EM_TSSNA_VIF5_TruncN20"

diff.perc.df <- readRDS(paste0(save.filename,".RData"))

ensembles <- c("EMca", "EMmean", "EMwmean", "EMmedian")

extracted_list <- list()

for (em in ensembles) {
  
  t_FP_col <- paste0(em, ".FP")
  t_FN_col <- paste0(em, ".FN")
  
  if (t_FP_col %in% colnames(diff.perc.df) && t_FN_col %in% colnames(diff.perc.df)) {
    
    em_display <- case_when(
      em == "EMca" ~ "Committee\naveraging",
      em == "EMmean" ~ "Mean",
      em == "EMwmean" ~ "Weighted\nmean",
      em == "EMmedian" ~ "Median",
      TRUE ~ em
    )
    
    fp_vector <- as.vector(as.numeric(unlist(diff.perc.df[[t_FP_col]])))
    fn_vector <- as.vector(as.numeric(unlist(diff.perc.df[[t_FN_col]])))
    
    tmp <- data.frame(
      Species = as.character(diff.perc.df$Species),
      EM = em_display,
      FP = fp_vector,
      FN = fn_vector,
      stringsAsFactors = FALSE
    )
    
    tmp$F <- tmp$FP + tmp$FN
    
    extracted_list[[em]] <- tmp
  }
}

analysis_df <- do.call(rbind, extracted_list)
analysis_df$EM <- as.factor(analysis_df$EM)
analysis_df$Species <- as.factor(analysis_df$Species)

target_metrics <- c("F", "FP", "FN")
titles <- c("(a) Overall performance",
            "(b) False positive     ",
            "(c) False negative     ")
filenames <- c("F", "FP", "FN")

# Create empty list to store the figures for ggarrange
fig_list <- list()

for (m in 1:3) {
  metric_name <- target_metrics[m]
  form <- as.formula(paste(metric_name, "~ EM | Species"))
  
  friedman_res <- friedman.test(form, data = analysis_df)
  cat("\n=====================================\n", titles[m], "\nGlobal Friedman p-value:", friedman_res$p.value, "\n")
  
  nemenyi_res <- frdAllPairsNemenyiTest(form, data = analysis_df)
  p_matrix <- nemenyi_res$p.value
  
  matrix_df <- as.data.frame(p_matrix) %>%
    mutate(EM1 = rownames(p_matrix)) %>%
    pivot_longer(-EM1, names_to = "EM2", values_to = "p_val") %>%
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
    ggplot(matrix_df, aes(x = EM2, y = EM1, fill = ColorKey)) +
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
    labs(subtitle = titles[m], x = "Ensemble method", y = "Ensemble method", fill = " ") +
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.margin = margin(0.3, 0.3, 0.3, 0.3, unit = "cm"),
      #plot.title    = element_text(size = 22, family = "Calibri", face = "bold", hjust = 0.5, vjust = 1),
      plot.title    = element_text(size = 16, family = "Calibri", hjust = 0, vjust = 1),
      plot.title.position = "plot",
      plot.subtitle = element_text(size = 16, family = "Calibri", hjust = 0.45),
      axis.title.x  = element_blank(),
      axis.title.y  = element_blank(),
      axis.text.x   = element_text(size = 13, family = "Calibri"),
      axis.text.y   = element_text(size = 13, family = "Calibri"),
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

combined_fig <- ggarrange(fig_list[[1]],fig_list[[2]],fig_list[[3]],
                          nrow = 2, ncol = 2, align = "v"#, 
                          #common.legend = TRUE, legend = "bottom"
)

print(combined_fig)
ggsave(paste0(save.filename,"_Combined.png"), plot = combined_fig, width = 20, height = 12, units = "cm", dpi = 800)



### Plot legend
legend <- ggplot(matrix_df, aes(x = EM2, y = EM1, fill = ColorKey)) +
  geom_tile(color = "white", linewidth = 0.5, show.legend = TRUE) +
  geom_text(aes(label = LabelAllValue), color = "black", fontface = "plain", size = 3.5, family = "Calibri") +
  scale_fill_manual(
    values = c("grey85", "#6EA6CD", "#4A7BB7", "#364B9A"),
    labels = c("Non-significant    ", "p < 0.05    ", "p < 0.01    ", "p < 0.001    "),
    drop = FALSE
  ) +
  scale_x_discrete(position = "bottom") +
  scale_y_discrete(limits = rev) +
  labs(subtitle = titles[m], x = "Ensemble method", y = "Ensemble method", fill = " ") +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.margin = margin(0.3, 0.3, 0.3, 0.3, unit = "cm"),
    plot.title    = element_text(size = 16, family = "Calibri", hjust = 0, vjust = 1),
    plot.title.position = "plot",
    plot.subtitle = element_text(size = 16, family = "Calibri", hjust = 0),
    axis.title.x  = element_blank(),
    axis.title.y  = element_blank(),
    axis.text.x   = element_text(size = 13, family = "Calibri"),
    axis.text.y   = element_text(size = 13, family = "Calibri"),
    axis.ticks = element_blank(), 
    legend.direction = "vertical",
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

ggsave("Sig_matrix_legend_vert.png", plot = legend, width = 25, height = 25, units = "cm", dpi = 800)

