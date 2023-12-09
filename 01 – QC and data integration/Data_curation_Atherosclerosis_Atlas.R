library(Seurat)
library(SeuratDisk)
library(dplyr)

setwd("/home/uv525372/data_integration_harmony/Human/")

#Alsaigh 2020: https://www.biorxiv.org/content/10.1101/2020.03.03.968123v1.full
Alsaigh <- Read10X("./Alsaigh_2020/")
Alsaigh <- CreateSeuratObject(counts = Alsaigh, project = "Alsaigh", names.field = 2, names.delim = "-")

meta <- Alsaigh@meta.data
meta <- meta %>% mutate(Tissue = case_when(orig.ident == 1 ~ "carotid artery (PA)",
                                           orig.ident == 2 ~ "carotid artery (AC)",
                                           orig.ident == 3 ~ "carotid artery (PA)",
                                           orig.ident == 4 ~ "carotid artery (AC)",
                                           orig.ident == 5 ~ "carotid artery (PA)",
                                           orig.ident == 6 ~ "carotid artery (AC)"))
Alsaigh@meta.data$Tissue <- meta[row.names(Alsaigh@meta.data),]$Tissue
Alsaigh@meta.data$Disease = "type VII calcified plaques endarterectomy"

meta <- meta %>% mutate(Donor = case_when(orig.ident == 1 ~ "P1",
                                          orig.ident == 2 ~ "P1",
                                          orig.ident == 3 ~ "P2",
                                          orig.ident == 4 ~ "P2",
                                          orig.ident == 5 ~ "P3",
                                          orig.ident == 6 ~ "P3"))
Alsaigh@meta.data$Donor <- meta[row.names(Alsaigh@meta.data),]$Donor
Alsaigh@meta.data$Study = "Alsaigh 2020"
Alsaigh@meta.data$Assay = "Cells"

saveRDS(Alsaigh, "./pre_processing/Alsaigh_raw.rds")


#Hu 2019: https://www.ahajournals.org/doi/full/10.1161/ATVBAHA.120.315373
Hu <- Read10X("./Hu_2021/")
Hu <- CreateSeuratObject(counts = Hu, project = "Hu", names.field = 2, names.delim = "-")

meta <- Hu@meta.data
meta$t <- substr(sub(".*_", "", Hu@meta.data$orig.ident), 1, 2)
meta <- meta %>% mutate(Tissue = case_when(t == "AO" ~ "Aorta",
                                           t == "PA" ~ "Pulmonary artery",
                                           t == "CA" ~ "Coronary artery"))
Hu@meta.data$Tissue <- meta[row.names(Hu@meta.data),]$Tissue
Hu@meta.data$Disease = "Control"
Hu@meta.data$Donor = substr(Hu@meta.data$orig.ident, 1, 2)
Hu@meta.data$Study = "Hu 2021"
Hu@meta.data$Assay = "Cells"

saveRDS(Hu, "./pre_processing/Hu_raw.rds")


#Pan 2020:https://www.ahajournals.org/doi/10.1161/CIRCULATIONAHA.120.048378#d1e2347
Pan <- readRDS("./Pan_2020/Pan_unprocessed_data.rds")
Pan@meta.data$Tissue = "carotid artery"
Pan@meta.data$Disease = "atherosclerotic plaques"
Pan@meta.data$Donor = Pan$orig.ident
Pan@meta.data$Study = "Pan 2020"
Pan@meta.data$Assay = "Cells"
saveRDS(Pan, "./pre_processing/Pan_raw.rds")


#Wirka 2019: https://www.nature.com/articles/s41591-019-0512-5#data-availability
#Cellranger with same reference genome
data <- Read10X("/data/uv525372/Atherosclerosis_Atlas/Wirka_data/rca1_ref3/outs/filtered_feature_bc_matrix/")
#Pa1 <- CreateSeuratObject(counts = data, project = "Wirka", names.field = 2, names.delim = "\\.")
Pa1 <- CreateSeuratObject(counts = data, project = "Wirka", names.field = 2, names.delim = "-", min.cells = 3, min.features = 200)
Pa1@meta.data$Donor = "Pa1"

data <- Read10X("/data/uv525372/Atherosclerosis_Atlas/Wirka_data/rca2_ref3/outs/filtered_feature_bc_matrix/")
Pa2 <- CreateSeuratObject(counts = data, project = "Wirka", names.field = 2, names.delim = "-", min.cells = 3, min.features = 200)
Pa2@meta.data$Donor = "Pa2"

data <- Read10X("/data/uv525372/Atherosclerosis_Atlas/Wirka_data/rca3_ref3/outs/filtered_feature_bc_matrix/")
Pa3 <- CreateSeuratObject(counts = data, project = "Wirka", names.field = 2, names.delim = "-", min.cells = 3, min.features = 200)
Pa3@meta.data$Donor = "Pa3"

data <- Read10X("/data/uv525372/Atherosclerosis_Atlas/Wirka_data/rca4_ref3/outs/filtered_feature_bc_matrix/")
Pa4 <- CreateSeuratObject(counts = data, project = "Wirka", names.field = 2, names.delim = "-", min.cells = 3, min.features = 200)
Pa4@meta.data$Donor = "Pa4"

Wirka <- merge(Pa1, y = c(Pa2, Pa3, Pa4), add.cell.ids = c("Pa1", "Pa2", "Pa3", "Pa4"), project = "Wirka")
Wirka@meta.data$Tissue = "Coronary artery"
Wirka@meta.data$Disease = "atherosclerotic lesions"
Wirka@meta.data$Study = "Wirka 2019"
Wirka@meta.data$Assay = "Cells"
Wirka@meta.data$Artery = "Coronary artery"
Wirka@meta.data$State = "Atherosclerosis"

meta <- Wirka@meta.data
meta <- meta %>% mutate(Sex = case_when(Donor == "Pa1" ~ "male",
                                           Donor == "Pa2" ~ "male",
                                           Donor == "Pa3" ~ "male",
                                           Donor == "Pa4" ~ "female",))
Wirka@meta.data$Sex <- meta[row.names(Wirka@meta.data),]$Sex
Wirka@meta.data$orig.ident <- NULL
saveRDS(Wirka, "./pre_processing/Wirka_raw.rds")