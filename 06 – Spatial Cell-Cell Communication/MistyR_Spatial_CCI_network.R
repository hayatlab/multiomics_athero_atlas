library(tidyverse)
library(Seurat)
library(mistyR)
library(Matrix)
source("/data/uv525372/Atherosclerosis_Atlas/Visium/misty_utilities.R")


setwd("/data/uv525372/Atherosclerosis_Atlas/")

#' Run MISTy for spatial interactions between slides
#'

#LIANA based interactions
cci_genes <- read.table("cci_genes.csv", sep=",", h=T)


future::plan(future::multisession)

# Misty pipeline definition:
run_colocalization <- function(slide, 
                               assay, 
                               useful_features, 
                               out_label, 
                               misty_out_alias = "./results/misty/cci_") {
  
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

folder <- "/data/uv525372/Atherosclerosis_Atlas/Visium/rawdata/"
slide_files <- list.files(folder)

#consider the slides with decent UMI counts
slide_files <- slide_files[4:12]


map(slide_files, function(slide_file){
  print(slide_file)
  
  slide_id <- slide_file
  print(slide_id)
  slide <- Load10X_Spatial(
    data.dir= paste(folder, slide_id, sep=""),
    filename = "filtered_feature_bc_matrix.h5")
  
  selected_features <- rownames(slide)[rownames(slide) %in% cci_genes$x]
  #selected_features <- rownames(slide)[rownames(slide) %in% cci_genes[1:100,]]
  slide <- subset(slide, features = selected_features)
  selected_features <- rownames(slide)[Matrix::rowSums(slide) > 250]
  slide <- subset(slide, features = selected_features)

  
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
  
  assays <-DefaultAssay(slide) 
  
  # Loop over each assay
  for (assay in assays) {
    # Set the default assay for the slide
    DefaultAssay(slide) <- assay
    
    # Get the useful features (row names)
    useful_features <- rownames(slide)
    
    print(useful_features)
    # Run colocalization analysis
    mout <- run_colocalization(slide = slide,
                               useful_features = useful_features,
                               out_label = slide_id,
                               assay = assay,
                               misty_out_alias = paste("./results/misty/", assay, "/", sep=""))
    
    # Collect results
    misty_res_slide <- collect_results(mout)

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
