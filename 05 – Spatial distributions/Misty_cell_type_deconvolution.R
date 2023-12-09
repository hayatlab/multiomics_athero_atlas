# Copyright (c) [2023] [Tore Bleckwehl]

library(tidyverse)
library(Seurat)
library(mistyR)
library(Matrix)
source("/data/uv525372/Atherosclerosis_Atlas/Visium/misty_utilities.R")


setwd("/data/uv525372/Atherosclerosis_Atlas/")

#' Run MISTy for spatial interactions between slides
#' 

future::plan(future::multisession)

# Pipeline definition:
run_colocalization <- function(slide, 
                               assay, 
                               useful_features, 
                               out_label, 
                               misty_out_alias = "./results/misty/main_") {
  
  # Define assay of each view ---------------
  view_assays <- list("main" = assay,
                      "juxta" = assay)
  # Define features of each view ------------
  view_features <- list("main" = useful_features, 
                        "juxta" = useful_features)
  # Define spatial context of each view -----
  view_types <- list("main" = "intra", 
                     "juxta" = "juxta")
  # Define additional parameters (l in case of paraview,
  # n of neighbors in case of juxta) --------
  view_params <- list("main" = NULL, 
                      "juxta" = 6)
  
  misty_out <- paste0(misty_out_alias, 
                      out_label, "_", assay)
  
  run_misty_seurat(visium.slide = slide,
                   view.assays = view_assays,
                   view.features = view_features,
                   view.types = view_types,
                   view.params = view_params,
                   spot.ids = NULL,
                   out.alias = misty_out)
  
  return(misty_out)
}


# Main ------------------------------------------------------------------------

# Getting sample annotations --------------------------------------------------

#sample_dict <- read.table("./markers/visium_annotations_ext.txt", 
#                          sep = "\t", header = T)

#slide_files_folder <- "./processed_visium/objects/"
#slide_files <- list.files(slide_files_folder)
#slide_ids <- gsub("[.]rds", "", slide_files)

# Cell2location proportions - complete
#assay_label <- "c2l_props"

#misty_outs <- map(slide_files, function(slide_file){
 # print(slide_file)
  # Read spatial transcriptomics data
  #slide_id <- gsub("[.]rds", "", slide_file)
  #slide <- readRDS(paste0(slide_files_folder, slide_file))
 
  #useful_features <- useful_features[! useful_features %in% "prolif"]

folder <- "/data/uv525372/Atherosclerosis_Atlas/Visium/rawdata/"
slide_files <- list.files(folder)
slide_files <- slide_files[4:12]
#importances <- data.frame()

map(slide_files, function(slide_file){
  print(slide_file)

  slide_id <- slide_file
  print(slide_id)
  slide <- Load10X_Spatial(
    data.dir= paste(folder, slide_id, sep=""),
    filename = "filtered_feature_bc_matrix.h5")

  
  #selected_features <- rownames(slide)[Matrix::rowSums(slide) > 1000]
  #slide <- subset(slide, features = selected_features)
  #normalize??
  

  # Read and process the CSV file with the cell2loc results
  read_c2l <- function(file_path, slide_id) {
    cell2loc <- read.table(file_path, sep = ",", header = TRUE)
    row.names(cell2loc) <- cell2loc$spot_id
    cell2loc$spot_id <- NULL
    cell2loc$Lymphatic.Endothelial <- NULL
    assay <- subset(cell2loc, cell2loc$sample == slide_id)
    rownames(assay) <- sub(paste("^", slide_id, "_", sep = ""), "", rownames(assay))
    assay$sample <- NULL
    assay$Disease <- NULL
    assay <- Matrix(as.matrix(assay), sparse = TRUE)
    return(assay)}
  
  # Spot deconvolutions for the Main Cell types
  c2lmain_assay <- read_c2l("/data/uv525372/Atherosclerosis_Atlas/c2l_main.csv", slide_id)
  slide[["c2lmain"]] <- CreateAssayObject(data = t(c2lmain_assay))

 
  # Spot deconvolutions for the Subcluster
  c2lsub_assay <- read_c2l("/data/uv525372/Atherosclerosis_Atlas/c2l_sub.csv", slide_id)
  slide[["c2lsub"]] <- CreateAssayObject(data = t(c2lsub_assay))

  # Define the assays to be used
  assays <- c("c2lmain", "c2lsub")
  
  # Loop over each assay
  for (assay in assays) {
    # Set the default assay for the slide
    DefaultAssay(slide) <- assay
    
    # Get the useful features (row names)
    useful_features <- rownames(slide)
    
    # Run colocalization analysis
    mout <- run_colocalization(slide = slide,
                               useful_features = useful_features,
                               out_label = slide_id,
                               assay = assay,
                               misty_out_alias = paste("./results/misty/", assay, "/", sep=""))
    
    # Collect results
    misty_res_slide <- collect_results(mout)
    #importances <- rbind(importances, misty_res_slide$importances)
    #saveRDS(importances, file = "./results/misty/importances.rds")
    
    # Save results to an RDS file
    saveRDS(misty_res_slide, file = paste0("./results/misty/", slide_id, "_", assay, "/misty_res_slide_", slide_id, ".rds"))
    
    # Create a plot folder
    plot_folder <- paste0(mout, "/plots")
    system(paste0("mkdir ", plot_folder))
    
    # Create a PDF file for summary plots
    pdf(file = paste0(plot_folder, "/", slide_id, "_", "summary_plots_", assay, ".pdf"))
    
    # Plot improvement stats and view contributions
    mistyR::plot_improvement_stats(misty_res_slide)
    mistyR::plot_view_contributions(misty_res_slide)
    
    # Plot interaction heatmap and communities for intra and juxta_5
    mistyR::plot_interaction_heatmap(misty_res_slide, "intra", cutoff = 0)
    mistyR::plot_interaction_communities(misty_res_slide, "intra", cutoff = 0.5)
    mistyR::plot_interaction_heatmap(misty_res_slide, "juxta_6", cutoff = 0)
    mistyR::plot_interaction_communities(misty_res_slide, "juxta_6", cutoff = 0.5)
    
    dev.off()
  }
  
})
  
