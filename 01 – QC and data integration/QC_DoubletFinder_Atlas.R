library(Seurat)
library(purrr)
library(dplyr)
library(ggplot2)
library(scater)
library(scDblFinder)
library(SeuratDisk)

#Load the different datasets
data_set <- list.files(pattern = "*raw.rds") %>% set_names() %>% map(readRDS)

pdf("QC_Atherosclerosis_samples.pdf")
for( i in seq_along(data_set)){
  data_set[[i]][["percent.mt"]] <- PercentageFeatureSet(data_set[[i]], pattern = "^MT-")
  print(VlnPlot(data_set[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by="Study", ncol = 3))
  data_set[[i]] <- subset(data_set[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < quantile(data_set[[i]]$nFeature_RNA, .98) & percent.mt < 20)
  

  #UMAP plots for each data set
  current_data <- NormalizeData(data_set[[i]])
  current_data <- FindVariableFeatures(current_data, selection.method = "vst", nfeatures = 2000)
  current_data <- ScaleData(current_data)
  current_data <- RunPCA(current_data, features = VariableFeatures(object = current_data))
  current_data <- RunUMAP(current_data, dims = 1:30)
  
  sce <- as.SingleCellExperiment(current_data)
  scores <- computeDoubletDensity(sce)
  current_data <- AddMetaData(current_data, scores, col.name = "Doublets")
  data_set[[i]] <- AddMetaData(data_set[[i]], scores, col.name = "Doublets")
  
  print(FeaturePlot(current_data, features = c("Doublets")))
  
  print(DimPlot(current_data, reduction = "umap",  group.by = "Donor"))+print(DimPlot(current_data, reduction = "umap",  group.by = "Donor"))+labs(title=sub("_raw.*", "", names(data_set)[i]), color="Donor")
  }

dev.off()

#Convert in a combined dataset for Seurat and Scanpy based integration method
merged <- merge(data_set[[1]], y = data_set[c(-1)], add.cell.ids = sub("_raw.*", "", names(data_set)), project = "integration")
merged@meta.data$orig.ident <- NULL
SaveH5Seurat(merged, filename = "raw_Atherosclerosis_atlas.h5Seurat", overwrite = T)
Convert("raw_Atherosclerosis_atlas.h5Seurat", dest = "h5ad", overwrite = T)
