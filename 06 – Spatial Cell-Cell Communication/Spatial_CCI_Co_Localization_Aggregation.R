library(tidyverse)
library(Seurat)
library(mistyR)
library(Matrix)
library(liana)
library(openxlsx)
library(RColorBrewer)
source("/data/uv525372/Atherosclerosis_Atlas/Visium/misty_utilities.R")

setwd("/data/uv525372/Atherosclerosis_Atlas/")

#Obtain the LIANA Ligand Receptor Interactions
#https://saezlab.github.io/liana/reference/get_curated_omni.html#details
#liana <- get_curated_omni(curated_resources = c("CellPhoneDB", "CellChatDB", "ICELLNET", "connectomeDB2020", "CellTalkDB"))

liana <- read.table("LIANA_db.csv", sep = ",", header = TRUE)
liana$LR <- paste0(pathways$source_genesymbol, "_", pathways$target_genesymbol, sep="")

# Functions to assign the highest expression to a cell type and evaluate the specficity
calculateSpecificityScore <- function(Atlas, features, group.by) {
  # Average expression from Seurat for the defined group
  AverageExpression(Atlas, features = features, group.by = group.by, slot = "counts") -> tmp
  avg_exp <- data.frame(t(tmp$RNA))
  
  # Calculate specificity score based on ratio of max to second max value
  specificity_scores <- sapply(avg_exp, function(col) {
    sorted_col <- sort(col, decreasing = TRUE)
    max_val <- sorted_col[1]  # Maximum value
    second_max_val <- sorted_col[2]  # Second maximum value
    specificity <- max_val / second_max_val  # Ratio of max to second max value
    return(specificity)
  })
  
  # Find row names corresponding to maximum values for each feature
  max_row_names <- apply(avg_exp, 2, function(col) {
    sorted_col <- sort(col, decreasing = TRUE)
    selected_indices <- which(col == sorted_col[1])
    rownames(avg_exp)[selected_indices]
  })
  
  return(list(Specificity_Scores = specificity_scores, Max_Row_Names = max_row_names))
}

assignSpecificityAndCellType <- function(Atlas, data, Interactors, group.by) {
  features <- data[[Interactors]]
  scores <- calculateSpecificityScore(Atlas, features = features, group.by = group.by)
  max_row_names <- scores$Max_Row_Names
  specificity_scores <- scores$Specificity_Scores
  
  cell_types <- max_row_names[match(features, names(max_row_names))]
  specificities <- specificity_scores[match(features, names(specificity_scores))]
  
  # Add the cell type and specificity score to the the tibble
  data[[paste(Interactors, group.by, sep="_")]] <- as.character(cell_types)
  data[[paste(Interactors, group.by, "Specificity", sep="_")]] <- specificities
  
  return(data)
}

# Aggregate the results and add the sample background
Sample_Disease <- read.table("c2l_main.csv", sep = ",", header = TRUE)
Sample_Disease <- Sample_Disease %>% select(sample, Disease) %>% distinct()

misty_out_folder <- "./results/misty/Spatial/"
misty_outs <- list.files(misty_out_folder, full.names = F)
misty_outs <- set_names(misty_outs)

# misty_outs[c(1,2)] -> misty_outs
misty_res <- collect_results(paste0(misty_out_folder, misty_outs))

sample_importances <- misty_res$importances %>% 
  mutate(sample = sub("^/data/uv525372/Atherosclerosis_Atlas/results/misty/Spatial/", "", sample) %>%
           gsub("_Spatial", "", .)) %>% left_join(Sample_Disease, by = "sample")

sample_importances$LR <- paste0(sample_importances$Predictor, "_", sample_importances$Target, sep="")


#Filtering for interaction from the liana database with the ligand as Predictor
sample_importances <- sample_importances %>%
  filter(LR %in% liana$LR)

#Excluding one sample
sample_importances <- sample_importances %>% filter(sample!="FW106014")

#Summarizing the co-localization results
summarized_interactions_group <- sample_importances %>%
  group_by(view, Predictor, Target, Disease) %>%
  summarize(mean_importance = mean(Importance), sd_importance=sd(Importance)) %>%
  ungroup()

# Re-adjust the mislabled genes by misty
summarized_interactions_group$Predictor <- gsub("\\.", "-", summarized_interactions_group$Predictor)
summarized_interactions_group$Target <- gsub("\\.", "-", summarized_interactions_group$Target)

#Load the scRNA-seq atlas
Atlas <- readRDS("Human_Atherosclerosis_Atlas.rds")
Atlas <- subset(Atlas, subset = (Main_Cell_Types!="cDC1" | Main_Cell_Types=="pDC" | Main_Cell_Types=="Osteoblastic cells"))

# Different specificity levels
Specificity_Assignment <- list(
  list(column = "Predictor", group = "Main_Cell_Types"),
  list(column = "Predictor", group = "Subcluster_Annotation"),
  list(column = "Target", group = "Main_Cell_Types"),
  list(column = "Target", group = "Subcluster_Annotation")
)


for (level in Specificity_Assignment ) {
  summarized_interactions_group <- assignSpecificityAndCellType(
    Atlas,
    summarized_interactions_group,
    level$column,
    level$group
  )
}

wb <- createWorkbook()
addWorksheet(wb, "Ligand_Predictor")
writeData(wb, sheet = "Ligand_Predictor", misty)
saveWorkbook(wb, file = "Supplementary Table 4 - Co_localization_Spatial_CCI_Ligand_Predictor.xlsx", overwrite = TRUE)
saveRDS(summarized_interactions_group, "Spatial_Misty_Ligand_Predictor.rds")

