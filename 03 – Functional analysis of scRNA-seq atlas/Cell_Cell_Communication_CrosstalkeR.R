#load libraries
library(Seurat)
library(dplyr)
library(patchwork)
library(zellkonverter)
library(SeuratDisk)
require(SingleCellExperiment)
require(liana)
require(tibble)
require(purrr)
require(reticulate)
library(magrittr)
library(SeuratObject)
require(tidyr)
library(tidyverse)
library(harmony)
require(biomaRt)
library(magrittr)
library(CrossTalkeR)
library(SeuratObject)
library(OmnipathR)
library(rmarkdown)
library(tinytex)
#read data
adata_Seurat = readRDS('/data/Kramann_lab/From_Tore/To_Sidrah/Human_Atherosclerosis.rds')
adata_Seurat <- subset(adata_Seurat, subset = (Main_Cell_Types != "Osteoblastic cells"))

#extract data counts and metadata on basis of conditin
data = as.SingleCellExperiment(adata_Seurat)
gp <-assays(data)[["counts"]][,colData(data)$Disease_State=='Control']
pp <-assays(data)[["counts"]][,colData(data)$Disease_State=='Atherosclerosis']

metadata_gp <-tibble::tibble(rownames(colData(data))[colData(data)$Disease_State=='Control'], as.character(colData(data)$Main_Cell_Types[colData(data)$Disease_State =='Control']))
metadata_pp<- tibble::tibble(rownames(colData(data))[colData(data)$Disease_State=='Atherosclerosis'],as.character(colData(data)$Main_Cell_Types[colData(data)$Disease_State=='Atherosclerosis']))


# create seurat object
good = CreateSeuratObject(counts = gp)
Idents(good) = metadata_gp$`as.character(...)`

poor = CreateSeuratObject(counts = pp)
Idents(poor) = metadata_pp$`as.character(...)`

#run liana for both conditions

good = NormalizeData(good, normalization.method = "LogNormalize", scale.factor = 10000)
liana_test <-  liana_wrap(good)

poor = NormalizeData(poor, normalization.method = "LogNormalize", scale.factor = 10000)
liana_test2 <- liana_wrap(poor)

#extract cellphonedb results for crosstalkR comparative analysis
control <-liana_test$cellphonedb%>%
  filter(pvalue<=0.05) %>%
  mutate(gene_A=gsub('_','*',ligand.complex)) %>%
  mutate(gene_B=gsub('_','*',receptor.complex)) %>%
  mutate(source=gsub("_","",source)) %>%
  mutate(target=gsub("_","",target)) %>%
  mutate(type_gene_A = rep('Ligand',length(.data$ligand))) %>%
  mutate(type_gene_B = rep('Receptor',length(.data$receptor))) 

acm <-liana_test2$cellphonedb %>% 
  filter(pvalue<=0.05) %>%
  mutate(gene_A=gsub('_','*',ligand.complex)) %>%
  mutate(gene_B=gsub('_','*',receptor.complex)) %>%
  mutate(source=gsub("_","",source)) %>%
  mutate(target=gsub("_","",target)) %>%
  mutate(type_gene_A = rep('Ligand',length(.data$ligand))) %>%
  mutate(type_gene_B = rep('Receptor',length(.data$receptor))) 

#create list for the selected results
l = list()
l[['Control']] = control
l[['Athero']] = acm

#Generate crosstalkR html report of comparative analysis
data_dis <- generate_report(lrpaths = l,
                            genes=NULL,
                            out_path='~/out/path/',
                            threshold=0,
                            out_file = 'result.html',
                            output_fmt = "html_document",
                            report = T,
                            sel_columns = c('source','target','gene_A','gene_B','type_gene_A','type_gene_B','lr.mean'))
