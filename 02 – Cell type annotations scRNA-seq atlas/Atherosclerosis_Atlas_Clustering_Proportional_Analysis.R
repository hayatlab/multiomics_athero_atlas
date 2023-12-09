library(Seurat)
library(ggrepel)
library(ggplot2)
library(tidyverse)

# Load the atherosclerosis atlas
Atlas <- readRDS("Human_Atherosclerosis_Atlas.rds")
data <- Atlas@meta.data

# Diseaste State compositions of the main cluster
color <- data.frame("Main_Cell_Types" = c("B cells", "Endothelial_1", "Endothelial_2", "Fibroblast", 
                                          "Lymphatic Endothelial", "Macrophage","Mast cells", "Neuronal cells", 
                                          "Neutrophils", "Osteoblastic cells", "Pericytes", "Plasma cells",
                                          "T cells", "VSMC", "cDC1", "pDC"),
                    "colors" = c('#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b',
                                 '#e377c2', '#b5bd61', '#17becf', '#aec7e8', '#ffbb78', '#98df8a',
                                 '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2'))
data <- inner_join(data, color, by ="Main_Cell_Types")

order <- data %>% group_by(Main_Cell_Types, colors) %>% summarise(Cells = n()) %>% arrange(Cells) %>% mutate(Main_Cell_Types=factor(Main_Cell_Types, levels=Main_Cell_Types))
data$Disease_State <- factor(data$Disease_State, levels = c("Control", "Proximal", "Atherosclerosis"))

pdf("Atheroscleriosis_Main_Cell_Cluster_Disease_State.pdf")
ggplot(data, aes(y=factor(Main_Cell_Types, level=c(order$Main_Cell_Types)), fill=Disease_State)) + geom_bar(position="fill")+
  theme_classic()+theme(legend.position = "top", line = element_blank(), axis.text = element_text(colour = "black")) +
  labs(x="Propositional Composition", y="", fill="Donor Background")+
  scale_x_continuous(labels = scales::percent)+scale_fill_brewer(palette = "Blues")+
  theme(axis.text.y = element_text(colour = order$colors, size = 12))
dev.off()


# Group by the different categories used for proportions
grouped_counts <- Atlas@meta.data %>%
  group_by(Artery, Disease_State, Subcluster_Annotation, Study) %>%
  summarise(Counts = n()) %>%
  group_by(Study) %>%
  mutate(Study_Total = sum(Counts)) %>%
  mutate(Normalized_Counts = Counts / Study_Total)

# Calculate averages for Artery
averages_artery <- grouped_counts %>%
  group_by(Artery, Subcluster_Annotation) %>%
  summarise(Average_Normalized_Counts = mean(Normalized_Counts)) %>%
  pivot_wider(names_from = Artery, values_from = Average_Normalized_Counts, values_fill = 0) %>%
  mutate(Ratio_Artery = log2((`Coronary artery` / `Carotid artery`) + 0.01))

# Calculate averages for Disease_State
averages_disease <- grouped_counts %>%
  group_by(Disease_State, Subcluster_Annotation) %>%
  summarise(Average_Normalized_Counts = mean(Normalized_Counts)) %>%
  pivot_wider(names_from = Disease_State, values_from = Average_Normalized_Counts, values_fill = 0.0001) %>%
  mutate(Ratio_Disease = log2((`Atherosclerosis` / `Control`) + 0.01))

# Merge the datasets
cluster_cell_counts <- Atlas@meta.data %>%
  group_by(Subcluster_Annotation) %>%
  summarise(Counts = n())

merged_averages <- merge(averages_disease, averages_artery, by = "Subcluster_Annotation")
subcluster_proportion <- merge(merged_averages, cluster_cell_counts, by = "Subcluster_Annotation")
# subcluster_proportion <- subset(subcluster_proportion, Counts > 200)  # If filtering by Counts is required


# Plotting the Subclustering proportions
pdf("Atherosclerosis_Atlas_Subclustering_Proportions.pdf")
ggplot(subcluster_proportion, aes(y=Ratio_Artery, x=Ratio_Disease, color=Subcluster_Annotation, label=Subcluster_Annotation, size=Counts))+geom_vline(xintercept = 0, col="gray", linetype = "dashed")+geom_hline(yintercept = 0, col="gray", linetype = "dashed")+geom_point()+
  theme(panel.background = element_blank(), legend.position = "none")+xlab(expression(Atherosclerosis %<->% Control))+ylab(expression(Carotid %<->% Coronary))+theme(axis.title = element_text(size = 16), axis.text = element_text(size = 12))+
  scale_color_manual(values = c('#ffff00', '#1ce6ff', '#ff34ff', '#ff4a46', '#008941', '#006fa6',
                                '#a30059', '#ffdbe5', '#7a4900', '#0000a6', '#63ffac', '#b79762',
                                '#004d43', '#8fb0ff', '#997d87', '#5a0007', '#809693', '#6a3a4c',
                                '#1b4400', '#4fc601', '#3b5dff', '#4a3b53', '#ff2f80', '#61615a',
                                '#ba0900', '#6b7900', '#00c2a0', '#ffaa92', '#ff90c9', '#b903aa',
                                '#d16100', '#ddefff'))+
  geom_text_repel(nudge_x = .15, box.padding = 0.5, nudge_y = 1, segment.curvature = -0.1, segment.ncp = 3, segment.angle = 20)


subcluster_proportion  <- subcluster_proportion %>% arrange(Ratio_Disease) %>% mutate(Subcluster_Annotation=factor(Subcluster_Annotation, levels=Subcluster_Annotation))
ggplot(subcluster_proportion, aes(x=Ratio_Disease, y=Subcluster_Annotation)) +
  geom_segment( aes(y=Subcluster_Annotation, yend=Subcluster_Annotation, x=0, xend=Ratio_Disease), color="grey") +
  geom_point(aes(color=Ratio_Artery, size=Counts))+scale_colour_gradient(low = "orange", high = "darkblue", breaks=c(-6,0,9),labels=c("Carotid", "0", "Coronary"), limits=c(-7,10))+
  labs(y="", x="Atherosclerosis Score", color="Artery \nScore")+theme_classic()+theme(axis.text = element_text(colour = "black", size =9))
dev.off()