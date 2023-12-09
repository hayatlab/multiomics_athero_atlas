library(Seurat)
library(tidyverse)
library(fgsea)
library(msigdbr)
library(stringr)
library(patchwork)
library(ggplot2)

#Define Reactome Pathways
pathwaysDF <- msigdbr("human", category="C2", subcategory = "CP:REACTOME")
pathway_Reactome <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)

#Function for pathway enrichment of DEG from scVI integrated data
pathway_enrichment_scVI <- function(cell_types, pathways) {
  pathway_scVI <- tibble()
  DE <- tibble()
  
  for (type in cell_types) {
    print(type)
    tmp <- read.table(paste0("DE_scVI_", type, ".csv"), sep = ",", header = TRUE, row.names = 1)
    tmp <- tmp[order(tmp$lfc_mean), ]
    ranks <- tmp$lfc_mean
    names(ranks) <- rownames(tmp)
    fgseaRes <- fgsea(pathways = pathways, stats = ranks, minSize = 15, maxSize = 500)
    fgseaRes$Main_Cell_Types <- type
    tmp$Main_Cell_Types <- type
    tmp$Gene <- rownames(tmp)
    rownames(tmp) <- NULL
    DE <- rbind(tmp, DE)
    pathway_scVI <- rbind(pathway_scVI, fgseaRes)
  }
  
  pathway_scVI$pathway <- str_remove(pathway_scVI$pathway, "^[^_]+_")
  return(list(DE = DE, pathway_scVI = pathway_scVI))
}


#Reactome Pathway of the main cluster
setwd("./DE_scVI_Main_Cell_Types/")
cell_types <- c("Fibroblast", "Macrophage", "Pericytes", "T cells", "VSMC", "Endothelial_1", "Endothelial_2", "Mast cells", "Neuronal cells")
pathway_scVI <- pathway_enrichment_scVI(cell_types = cell_types, pathways = pathway_Reactome)

#Convert the results to csv
pathway_scVI$Pathway_Genes <- sapply(pathway_scVI$leadingEdge, paste, collapse = ",")
pathway_scVI$leadingEdge <- NULL
write.table(as.data.frame(pathway_scVI), "Reactome_pathways_FGSEA_DE_scVI.csv", quote = F, row.names = F)

signif <- subset(pathway_scVI$pathway_scVI, pval<0.05)

subset_df <- pathway_scVI$pathway_scVI[pathway_scVI$pathway_scVI$pathway %in% signif$pathway, ]
subset_df$pathway[subset_df$pathway == "REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_TRANSPORT_AND_UPTAKE_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS"] <- "TRANSPORT_AND_UPTAKE_BY_IGFBPS"
subset_df$NES <- subset_df$NES*(-1)
path_plot <- ggplot(subset_df, aes(y=fct_inorder(pathway), x=Main_Cell_Types,color=NES, size=pval))+geom_point()+coord_flip()+
  scale_color_distiller(palette = "RdBu") +scale_size(range = c(0.0005, 5), breaks = c( 0.0005, 0.005, 0.05, 0.5), trans = 'reverse')+
  theme_classic()+theme(plot.title = element_text(colour = "black", size =12), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, colour = "black"), axis.text.y = element_text(size=9, colour = "black"))+
  labs(x="", y="", title = "Reactome Pathways - Control vs. Atherosclerosis")

DE_genes <- pathway_scVI$DE %>% group_by(Main_Cell_Types) %>%
  summarize(up_Atheroscloris = sum(lfc_mean > 0),
            up_ctrl = sum(lfc_mean < 0)) %>% 
  mutate(up_ctrl = up_ctrl*(-1)) %>%
  pivot_longer(cols = c(up_Atheroscloris, up_ctrl),
               names_to = "Count_Type",
               values_to = "Count"
  )
DE_plot <- ggplot(DE_genes, aes(x=Count, y=Main_Cell_Types, fill=Count_Type))+geom_bar(stat = "identity")+NoLegend()+labs(y="", x="", title="DE Genes")+theme_classic()+theme(plot.title = element_text(colour = "black", size =12), axis.text = element_text(colour = "black", size =9), line = element_blank(), legend.position = "")+scale_fill_manual(values = c("#D6604D", "#4393C3"))
pdf("Atherosclerosis_DEG_Reactome_Pathways.pdf")
wrap_plots(DE_plot, path_plot, widths = c(1,5))
dev.off()



## Subclustering in which clusters are compred to each other
pathway_enrichment_by_subcluster <- function(data, pathways) {
  pathway_scVI <- tibble()
  DE <- tibble()
  
  for (i in unique(data$group1)) {
    tmp <- data %>% filter(group1 == i)
    rownames(tmp) <- tmp$X
    tmp$X <- NULL
    tmp <- tmp[order(tmp$lfc_mean), ]
    ranks <- tmp$lfc_mean * tmp$bayes_factor
    names(ranks) <- rownames(tmp)
    fgseaRes <- fgsea(pathways = pathways, stats = ranks, minSize = 15, maxSize = 500)
    fgseaRes$Cluster <- i
    tmp$Gene <- rownames(tmp)
    rownames(tmp) <- NULL
    DE <- rbind(tmp, DE)
    pathway_scVI <- rbind(pathway_scVI, fgseaRes)
  }
  
  pathway_scVI$pathway <- str_remove(pathway_scVI$pathway, "^[^_]+_")
  return(list(DE = DE, pathway_scVI = pathway_scVI))
}

Subcluster_type <- "Macrophage"
Subcluster <- read.table(paste0("./DE_scVI_Atherosclerosis_",Subcluster_type,"_Subclustering.csv", sep=""), sep=",", h=T)
pathway_scVI <- pathway_enrichment_by_subcluster(Macrophages, pathway_Reactome)

p <- subset(pathway_scVI$pathway_scVI, pval<0.05)
subset_df <- pathway_scVI$pathway_scVI[pathway_scVI$pathway_scVI$pathway %in% p$pathway, ]

pdf(paste0("Atherosclerosis_Subclustering_", Subcluster_type, "Reactome_Pathways.pdf", sep=""))
ggplot(subset_df, aes(x=fct_inorder(pathway), y=Cluster,color=NES, size=pval))+geom_point()+coord_flip() + scale_color_distiller(palette = "RdBu") +scale_size(range = c(0.0005, 5), breaks = c( 0.0005, 0.005, 0.05, 0.5), trans = 'reverse')+theme_classic()+theme(plot.title = element_text(colour = "black", size =12), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, colour = "black"), axis.text.y = element_text(size=9, colour = "black"))+labs(x="", y="", title = "Reactome Pathways")
dev.off()
