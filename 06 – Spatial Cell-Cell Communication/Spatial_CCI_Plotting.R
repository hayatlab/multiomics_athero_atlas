library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)

misty <- readRDS("./Atherosclerosis/Spatial_Misty.rds") # For Ligand as Predictor
misty$LR_pair <- paste(misty$Predictor, misty$Target, sep="->")
misty <- misty %>% group_by(LR_pair, view) %>% add_count(LR_pair, name = "LR_count")

########.Cololaization ranking by cell types involved

plot_spatial_colocalization <- function(spatial_interactions) {
  # Aggregating and summarizing data
  spatial_summary <- spatial_interactions %>%
    group_by(Predictor, Target, Disease) %>%
    summarize(mean_importance = mean(mean_importance), sd_importance = sd(mean_importance)) %>%
    ungroup() %>%
    arrange(desc(mean_importance)) %>%
    filter(mean_importance > 0) %>%
    mutate(LR_pair = paste(Predictor, Target, sep = "->"))
  
  # Aggregating and reordering LR pairs
  LR_pair_agg <- spatial_summary %>%
    group_by(LR_pair) %>%
    summarize(order_mean = median(mean_importance)) %>%
    arrange(order_mean)
  
  # Reordering levels of 'LR_pair' based on aggregated 'order_mean'
  spatial_summary$LR_pair <- factor(spatial_summary$LR_pair, levels = unique(LR_pair_agg$LR_pair))
  
  # Create plot of the ranked CCI colocalizations
  ggplot(spatial_summary, aes(x = mean_importance, y = LR_pair, fill = Disease)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_classic() +
    theme(
      legend.position = "top",
      axis.text = element_text(size = 6, color = "black"),
      strip.text.y = element_text(size = 12),
      axis.title = element_text(size = 12)
    ) +
    labs(x = "LR co-expression [Mean importance]", y = "Spatial LR interaction of VSMCs and Fibroblasts") +
    scale_fill_manual(values = c('#CD113B', '#52006A', '#FF7600', "gray")) +
    scale_x_continuous(expand = c(0, 0))
}



## Fibroblast and VSMC
spatial_interactions <- misty %>%
  filter(Target_Main_Cell_Types %in% c("VSMC", "Fibroblast") | Predictor_Main_Cell_Types %in% c("VSMC", "Fibroblast"))

filtered_summary <- spatial_interactions %>%
  filter(Disease != "Intermediate\nlesion", mean_importance > 0.5) # Cut-off only for visulization purpose

plot_spatial_colocalization(filtered_summary)

## Macrophages

spatial_interactions <- misty %>%
  filter(Target_Main_Cell_Types=="Macrophage" | Predictor_Main_Cell_Types=="Macrophage")
filtered_summary <- spatial_interactions %>%
  filter(Disease != "Intermediate\nlesion", mean_importance > 0) # Cut-off only for visulization purpose
plot_spatial_colocalization(filtered_summary)



########. Co-localization and importance plot

# Filter the resulting CCI colocalizations.
# The value is abitery, to reduce the number of CCI for the plot
filtered_data <- misty %>%
  filter(Predictor_Main_Cell_Types_Specificity > 1.3 & Target_Main_Cell_Types_Specificity > 1.3)

# Aggregate mean values for 'Predictor' column
agg_predictor <- filtered_data %>%
  group_by(Predictor) %>%
  summarise(order_mean = mean(Predictor_Main_Cell_Types_Specificity)) %>%
  arrange(order_mean)

# Reorder levels of 'Predictor' based on aggregated 'order_mean'
filtered_data$Predictor <- factor(filtered_data$Predictor, levels = agg_predictor$Predictor)

# Aggregate mean values for 'Target' column
agg_target <- filtered_data %>%
  group_by(Target) %>%
  summarise(order_mean = mean(Target_Main_Cell_Types_Specificity)) %>%
  arrange(order_mean)

# Reorder levels of 'Target' based on aggregated 'order_mean'
filtered_data$Target <- factor(filtered_data$Target, levels = agg_target$Target)

# Define colors for 'Predictor_Main_Cell_Types'
color_df <- data.frame(
  "Predictor_Main_Cell_Types" = c("B cells", "Endothelial_1", "Endothelial_2", "Fibroblast", 
                                  "Lymphatic Endothelial", "Macrophage", "Mast cells", "Neuronal cells", 
                                  "Neutrophils", "Osteoblastic cells", "Pericytes", "Plasma cells",
                                  "T cells", "VSMC", "cDC1", "pDC"),
  "colors" = c('#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b',
               '#e377c2', '#b5bd61', '#17becf', '#aec7e8', '#ffbb78', '#98df8a',
               '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2')
)

# Join color information with the dataset
filtered_data <- inner_join(filtered_data, color_df, by = "Predictor_Main_Cell_Types")

# Create DotPlot (Seurat) and Scatter Plot (Co-localization).
# Genes are ordered by Cell Type Specificty Score
dot_plot <- DotPlot(Atlas, features = levels(filtered_data$Target), group.by = "Main_Cell_Types") +
  coord_flip() + RotatedAxis() +
  theme(legend.position = "left") +
  theme(axis.text = element_text(size = 3), legend.title = element_text(size = 10))

colocalization_plot <- ggplot(filtered_data, aes(x = Predictor, y = Target, color = Predictor_Main_Cell_Types, size = mean_importance)) +
  geom_point() +
  theme(panel.background = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 3), legend.title = element_text(size = 10))

# Display combined plot
dot_plot + colocalization_plot