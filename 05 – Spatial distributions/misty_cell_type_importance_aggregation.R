library(tidyverse)
library(Seurat)
library(mistyR)
library(Matrix)
source("/data/uv525372/Atherosclerosis_Atlas/Visium/misty_utilities.R")


setwd("/data/uv525372/Atherosclerosis_Atlas/")

Sample_Disease <- read.table("c2l_main.csv", sep = ",", header = TRUE)
Sample_Disease <- Sample_Disease %>% select(sample, Disease) %>% distinct()

misty_out_folder <- "./results/misty/c2lmain/"
misty_outs <- list.files(misty_out_folder, full.names = F)
misty_outs <- set_names(misty_outs) #%>% gsub("_c2lsub", "", .)

misty_res <- collect_results(paste0(misty_out_folder, misty_outs))

sample_importances <- misty_res$importances %>% 
  mutate(sample = sub("^/data/uv525372/Atherosclerosis_Atlas/results/misty/c2lmain/", "", sample) %>%
           gsub("_c2lmain", "", .)) %>% left_join(Sample_Disease, by = "sample")

#Mean importance for each sample of the intra View
summarized_interactions_group <- sample_importances %>%
  group_by(view, Predictor, Target, sample) %>%
  summarize(median_importance = mean(Importance)) %>%
  ungroup() %>%
  group_by(view)

summarized_interactions_sample <- summarized_interactions_group %>%
  filter(view =="intra") %>%
  ggplot(aes(x = Target, y = Predictor, fill = median_importance)) +
  geom_tile() + theme(axis.text = element_text(size = 5))+
  scale_fill_gradient2(high = "blue", 
                       midpoint = 0,
                       low = "white",
                       na.value = "grey") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_wrap(view ~ sample, scales = "free")

#Mean importance for Disease_State and View
summarized_interactions_group <- sample_importances %>%
  group_by(view, Predictor, Target, Disease) %>%
  #group_by(Predictor, Target, Disease) %>%
  summarize(mean_importance = mean(Importance)) %>%
  ungroup() 

summarized_interactions_group$Disease <- factor(summarized_interactions_group$Disease, levels = c("Control", "Intermediate\nlesion", "Atheroma"))
summarized_interactions <- summarized_interactions_group %>% filter(Disease != "Intermediate\nlesion")  %>%
  ggplot(aes(x = Target, y = Predictor, fill = mean_importance)) +
  geom_tile() +
  scale_fill_gradient2(high = "#0080bf", 
                       midpoint = 0,
                       low = "white",
                       na.value = "gray") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title = element_text(size=12)) +
  facet_grid(Disease ~ view)+coord_equal()+labs(fill="mean\nimportance")


# Comparative mean importance
Comparative <- summarized_interactions_group %>%
  filter(Disease!="Intermediate\nlesion") %>%
  group_by(Predictor, Target) %>%
  mutate(ratio = mean_importance - lead(mean_importance)) %>%
  filter(!is.na(ratio))

Comparative_Importance_plot <- ggplot(Comparative, aes(x = Target, y = Predictor, size=mean_importance, color = ratio)) +
  geom_point() + theme(panel.background = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title = element_text(size=12)) +
  scale_color_distiller(palette = "RdBu")+
  coord_equal()+labs(color="Importance\nAtheroma\nvs.\nControl", size="Mean\nImportance\nAtheroma")

# R2 distributions
R2_data <- misty_res$improvements %>%
  dplyr::filter(measure == "multi.R2") %>%
  mutate(sample = sub("^/data/uv525372/Atherosclerosis_Atlas/results/misty/c2lmain/", "", sample) %>%
           gsub("_c2lmain", "", .)) %>% left_join(Sample_Disease, by = "sample")

cell_order <- R2_data %>% 
  group_by(target) %>%
  summarize(med_value = median(value)) %>%
  arrange(-med_value) %>%
  pull(target)

cells_R2_tile <- ggplot(R2_data, aes(x = factor(target,
                                                levels = cell_order), 
                                     y = sample, fill = value)) +
  geom_tile() +
  coord_equal() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_fill_gradient(low = "black", high = "yellow")

cells_R2_box <- ggplot(R2_data, aes(x = factor(target,
                                               levels = cell_order), y = value)) +
  geom_boxplot() +
  geom_point(aes(color = sample)) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size =12),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  ylab("Explained variance") +
  xlab("")



#Summary plots
pdf("Misty_importance_Main_Cell_Types.pdf", width = 6, height = 6)

plot(summarized_interactions_sample)
plot(summarized_interactions)
plot(Comparative_Importance_plot)
mistyR::plot_interaction_communities(misty_res, "intra", cutoff = 0.5)
mistyR::plot_interaction_communities(misty_res, "juxta_6", cutoff = 0.5)
plot(cells_R2_box)
plot(cells_R2_tile)

dev.off()

