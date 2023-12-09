library(CrossTalkeR)
library(CrossTalkeR)
library(igraph)
library(ggraph)
library(ggplot2)
library(tidyverse)


CrossTalkeR <- readRDS("./CellCellCommunication/CCI_Atherosclerosis_CrosstalkeR.rds")

colors = c("B cells"= '#1f77b4', "Endothelial1" = '#ff7f0e', "Endothelial2" = '#279e68', "Fibroblast" = '#d62728', "Lymphatic Endothelial"='#aa40fc',
           "Macrophage" = '#8c564b', "Mast cells" = '#e377c2', "Neuronal cells"='#b5bd61', "Neutrophils" = '#17becf', "Osteoblastic cells" = '#aec7e8',
           "Pericytes" = '#ffbb78', "Plasma cells" = '#98df8a', "T cells" = '#ff9896', "VSMC" = '#c5b0d5', "cDC1" = '#c49c94', "pDC" = '#f7b6d2')

CrossTalkeR@colors <- colors
plot_cci(graph = CrossTalkeR@graphs$Athero_x_Control,
         colors = CrossTalkeR@colors,
         plt_name = 'CCI Atherosclerosis',
         coords = CrossTalkeR@coords[V(CrossTalkeR@graphs$Athero_x_Control)$name,],
         emax = NULL,
         leg = TRUE,
         low = 0,
         high = 0,
         ignore_alpha = FALSE,
         log = FALSE,
         efactor = 8,
         vfactor = 12,pg=CrossTalkeR@rankings$Athero_x_Control$Influencer)

# Extracting the CrosstalkeR results
CCI <- CrossTalkeR@tables$Athero_x_Control %>%
  mutate(
    LR = paste(Ligand, Receptor, sep = "_"),
    Interaction = paste(Ligand.Cluster, Receptor.Cluster, sep = "_")
  ) %>%
  group_by(LR) %>%
  add_count(LR, name = "LR_count")

# The LR count defines the number of interactions. Hence a lower number of interactions represents
# less borad interactions in the data set and will be applied for filtering.

# Define excluded cell types
exclude_cell_types <- c("Neuronal cells", "VSMC", "Fibroblast", "Lymphatic Endothelial")

# Filter and plot for Ligand expressing cells
Endo_Ligand <- CCI %>%
  filter(str_detect(Ligand.Cluster, str_c(c("Endothelial", "Pericytes"), collapse = "|")),
         !str_detect(Interaction, str_c(exclude_cell_types, collapse = "|")),
         !str_detect(Receptor.Cluster, str_c(c("Endothelial", "Pericytes"), collapse = "|"))) %>%
  group_by(Interaction) %>%
  filter(LR_count < 40, LRScore > 0) %>%
  ggplot(aes(y = LRScore, x = Ligand.Cluster, fill = Ligand.Cluster)) +
  geom_boxplot() +
  facet_grid(. ~ Receptor.Cluster, scales = "free") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "", legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Ligand expressing cells") +
  scale_fill_manual(values = c('#ff7f0e', '#279e68', '#ffbb78'))

# Filter and plot for Receptor expressing cells
Endo_Receptor <- CCI %>%
  filter(str_detect(Receptor.Cluster, str_c(c("Endothelial", "Pericytes"), collapse = "|")),
         !str_detect(Interaction, str_c(exclude_cell_types, collapse = "|")),
         !str_detect(Ligand.Cluster, str_c(c("Endothelial", "Pericytes"), collapse = "|"))) %>%
  group_by(Interaction) %>%
  filter(LR_count < 40, LRScore > 0) %>%
  ggplot(aes(y = LRScore, x = Receptor.Cluster, fill = Receptor.Cluster)) +
  geom_boxplot() +
  facet_grid(. ~ Ligand.Cluster, scales = "free") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "", legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Receptor expressing cells") +
  scale_fill_manual(values = c('#ff7f0e', '#279e68', '#ffbb78'))

pdf("Lymphocyte_EC_Pericyite_CCI.pdf")
print(Endo_Ligand)
print(Endo_Receptor)
dev.off()