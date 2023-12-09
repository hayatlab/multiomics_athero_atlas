library(Seurat)
library(SCP)
library(ggplot2)
library(tidyverse)

setwd("/data/uv525372/Atherosclerosis_Atlas")
Atlas <- readRDS("Human_Atherosclerosis_Atlas.rds")

#subset for the VSMC cluster
VSMC <- subset(Atlas, subset = Main_Cell_Types == "VSMC")
VSMC <- NormalizeData(VSMC)
VSMC <- RunSlingshot(srt = VSMC, group.by = "Subcluster_Annotation", end="Fibromyocytes", reduction = "scVI")


VSMC <- RunDynamicFeatures(srt = VSMC, lineages = "Lineage1", n_candidates = 4000)

ht <- DynamicHeatmap(
  srt = VSMC, lineages = "Lineage1",
  use_fitted = TRUE, n_split = 5,
  species = "Homo_sapiens", db = "GO_BP", anno_terms = TRUE, anno_keys = TRUE, anno_features = TRUE,
  heatmap_palette = "viridis",
  separate_annotation = list("Subcluster_Annotation", c("ACTC1", "TNFRSF11B")), separate_annotation_palette = c("Paired", "Set1"),
  pseudotime_label = 25, pseudotime_label_color = "red",
  height = 5, width = 2)

#Prepare the data for the Pseudotime analysis by Disease State
df <- data.frame(Disease=VSMC$Disease_State, 
                 Subclustering = VSMC$Subcluster_Annotation,
                 PseudoTime = VSMC$Lineage1)
df$Subclustering <- factor(df$Subclustering, levels=c("VSMC", "Myofibroblast", "Fibromyocytes"))


pdf("VSMC_Trajectory_Atherosclerosis.pdf", width = 15)
ggplot(df %>% filter(Disease!="Proximal"), aes(x=PseudoTime, y=Disease, color=Subclustering))+geom_jitter(size=0.3)+
  theme(panel.background = element_blank())+theme(axis.title = element_text(size = 21), axis.text = element_text(size = 16), legend.position = "top")+
  scale_color_brewer(palette = "Paired")+
  labs(y="", color="")+guides(color = guide_legend(override.aes = list(size = 8)))
CellDimPlot(VSMC, group.by = "Subcluster_Annotation", reduction = "umap", lineages = paste0("Lineage", 1), lineages_span = 0.1)
print(ht$plot)
dev.off()