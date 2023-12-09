library(DESeq2)
library(tibble)
library(dplyr)
library(tidyr)
library(ggrepel)
library(EnhancedVolcano)

# Load the atherosclerosis atlas
Atlas <- readRDS("Human_Atherosclerosis_Atlas.rds")
data <- Atlas@meta.data

Atlas <- subset(Atlas, subset = Disease_State != "Proximal")
subset<- subset(Atlas, subset = (Main_Cell_Types == "Endothelial_1" | Main_Cell_Types == "Endothelial_2"))

pseudobulk_DESeq2 <- function(subset, Entity, reference) {
  cts <- AggregateExpression(subset, group.by = Entity, slot = "counts", return.seurat = FALSE)
  cts <- cts$RNA
  colData <- data.frame(samples = colnames(cts))
  
  # Extract the "condition" column
  colData <- colData %>% mutate(condition = str_split(samples, "_", simplify = TRUE)[, 1])
  colData$condition <- sub("_[^_]+$", "", colData$samples)
  rownames(colData) <- colData$samples
  colData$samples <- NULL
  
  dds <- DESeqDataSetFromMatrix(countData = cts, colData = colData, design = ~condition)
  dds$condition <- relevel(dds$condition, ref = reference)
  
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  dds <- DESeq(dds)
  
  result_col <- resultsNames(dds)[2]
  print(result_col <- resultsNames(dds)[2])
  res <- results(dds, name = result_col)
  res <- data.frame(res)
  res$Gene = rownames(res)
  
  return(res)
}

#Process the pseudobulk DEG for the Endothelial in respect to the cell state and atherosclerosis
res_Cells <- pseudobulk_DESeq2(subset, c("Main_Cell_Types", "Donor"), "Endothelial_1")
res_Disease <- pseudobulk_DESeq2(subset, c("Disease_State", "Donor"), "Control")

# Join the DEG of the cell and disease state
res <- full_join(res_Disease, res_Cells, by = "Gene", suffix = c("_Disease", "_Cells"))
res <- subset(res, res$log2FoldChange_Disease<(20))

# Add colors according to the 
res$padj_Disease[is.na(res$padj_Disease)] <- 0
res$padj_Cells[is.na(res$padj_Cells)] <- 0
res$DEG <- ifelse(res$padj_Disease < 0.1 & res$padj_Cells < 0.1, "Shared",
            ifelse(res$padj_Cells < 0.1, "Cell",
                          ifelse(res$padj_Disease < 0.1, "Disease",NA)))


ggplot(res %>% filter(DEG!="NA"), aes(x=log2FoldChange_Disease, y=log2FoldChange_Cells, color=DEG))+geom_point(size=0.03)+
  theme_minimal()+coord_equal()+theme(axis.text = element_text(colour = "black"), legend.position = "top")+
  guides(colour = guide_legend(override.aes = list(size=5)))+scale_color_brewer(palette = "Dark2")

EC_marker <-  res %>% filter(DEG=="Shared") %>% filter(baseMean_Disease >5) %>% filter(baseMean_Cells >5) %>% 
  filter(log2FoldChange_Disease>0) %>% filter(padj_Disease >0) %>% arrange(padj_Disease*padj_Cells) %>% slice_head(n=20)

r <- res %>% mutate(plotname = as.character(Gene)) %>%
  mutate(plotname = ifelse(plotname %in% c(EC_marker$Gene, "LRRC75A", "PLCG2", "FOXC2", "ITGA6"), plotname, ""))

ggplot(r %>% filter(DEG!="NA"), aes(x=log2FoldChange_Disease, y=log2FoldChange_Cells, color=DEG))+geom_point(size=0.03)+
  theme_minimal()+coord_equal()+theme(axis.text = element_text(colour = "black"), legend.position = "top")+
  geom_label_repel(aes(label = plotname), box.padding = 0.5, size=2, max.overlaps = 30)+
  guides(colour = guide_legend(override.aes = list(size=5)))+scale_color_brewer(palette = "Dark2")

# Genes without significant changes are filtered out with filter(DEG!="NA")
ggplot(r %>% filter(DEG!="NA"), aes(x=log2FoldChange_Disease, y=log2FoldChange_Cells, color=DEG, size=baseMean_Cells))+geom_point()+
  theme_minimal()+theme(axis.text = element_text(colour = "black"), legend.position = "top")+
  geom_label_repel(aes(label = plotname), size = 2, box.padding = 0.5, size=2, max.overlaps = 50)+
  guides(colour = guide_legend(override.aes = list(size=5)))+scale_color_manual(values=c("#27496D", "#0C7B93", "#00A8CC"))+xlim((-7),7)+
  labs(x="Log2 FoldChange - Atherosclerosis vs Control", y="Log2 FoldChange - EC2 vs EC1", size="baseMean")


selected_pathway <- pathways$REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION
s <- res_Cells[res_Cells$Gene %in% selected_pathway, ]

EnhancedVolcano(res_Cells,lab = res_Cells$Gene, x = 'log2FoldChange', y = 'pvalue')

EnhancedVolcano(s, lab = rownames(s), x = 'log2FoldChange',
                y = 'pvalue',
                xlab = bquote(~Log[2]~ 'fold change'),
                pointSize = 3.0,
                labSize = 3.0,
                pCutoff = 0.05,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = TRUE,
                col = c('black', "#00A8CC", "#27496D", "#0C7B93"),
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black') + coord_flip()

Pericytes <- subset(Atlas, subset = (Main_Cell_Types == "Pericytes"))
res_Peri <- pseudobulk_DESeq2(Pericytes, c("Disease_State", "Donor"), "Control")
selected_pathway <- pathways$REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_TRANSPORT_AND_UPTAKE_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS
s <- res_Peri[res_Peri$Gene %in% selected_pathway, ]
EnhancedVolcano(res_Peri,lab = res_Peri$Gene, x = 'log2FoldChange', y = 'pvalue')
EnhancedVolcano(s, lab = rownames(s), x = 'log2FoldChange',
                y = 'pvalue', ylim = c(0,8),
                xlab = bquote(~Log[2]~ 'fold change'),
                pointSize = 3.0,
                labSize = 3.0,
                pCutoff = 0.05,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = TRUE,
                col = c('black', "#E25E3E", "#FF9B50", "#C63D2F"),
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black') + coord_flip()
