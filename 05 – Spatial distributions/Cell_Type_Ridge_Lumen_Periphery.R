library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# Function to generate the visualizations
generate_plot <- function(df, cell_names) {
  df <- df %>%
    select(spot_id, Disease, sample, lumina_to_periphery, !!!cell_names)
  
  new_df <- df %>%
    filter(!sample %in% c("FW106014")) %>%
    filter(!Disease %in% c("Fibroatheroma")) %>%
    pivot_longer(cols = all_of(cell_names), names_to = "Name", values_to = "Value") %>%
    group_by(Name) %>%
    mutate(scaled_value = scale(Value, center = FALSE, scale = max(Value) - min(Value))) %>%
    ungroup() %>%
    filter(Value > 0.4)
  
  new_df$Name <- factor(new_df$Name, levels = cell_names)
  
  p1 <- ggplot(new_df, aes(y = Name, x = lumina_to_periphery, fill = ..x..)) +
    geom_density_ridges_gradient(scale = 0.9) +
    scale_fill_distiller(palette = "YlOrRd") +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 12)
    ) +
    labs(x = "Lumina to Periphery", y = "")
  
  new_df$Disease <- factor(new_df$Disease, levels = c("Control", "Atheroma", "Intermediate \n lesion"))
  
  p2 <- ggplot(new_df, aes(y = Name, x = lumina_to_periphery, color = Disease, size = scaled_value)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7)) +
    theme_classic() +
    theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12)) +
    labs(x = "Lumina to Periphery", y = "", size = "Abundance") +
    scale_color_manual(values = c('#52006A', '#CD113B', '#FF7600'))
  
  ggarrange(p1, p2, widths = c(1, 2))
}


# Load cell2location cell type abundancy and lumen to periphery score
visium <- read_csv("./Atherosclerosis/c2l_meta_main_Atherosclerosis.csv")
visium_sub <- read_csv("./Atherosclerosis/c2l_meta_sub_Atherosclerosis.csv")

# Add luminal score to Main_Cell_Type spot deconvolution
tmp <- visium_sub %>% select(spot_id, subsample, lumina_to_periphery)
visium <- inner_join(visium, tmp, by="spot_id")


# Visium sample stats by disease state
sample_stats <- visium_sub %>% group_by(subsample, Disease, sample) %>% summarise(mean_value = sum(in_tissue, na.rm = TRUE)) %>% filter(mean_value>150) %>% group_by(sample) %>% add_count(sample, name = "tissue_count")
sample_stats$Disease <- factor(sample_stats$Disease, levels = c("Control", "Intermediate \n lesion", "Atheroma", "Fibroatheroma"))
ggplot(sample_stats, aes(y=sample, x=tissue_count, fill=Disease))+ geom_bar(position="dodge", stat="identity")+facet_grid(Disease~., scales = "free")+
  theme_classic()+theme(
    legend.position="none",
    axis.text = element_text(size=12, color = "black"), axis.text.y = element_blank(), strip.text.y = element_text(size = 18), axis.title = element_text(size=21))+labs(x="Arteries per slide", y="10X Visium Slides from human coronary arteries")+
  scale_fill_manual(values = c('#52006A','#FF7600','#CD113B',"gray"))


# Cell types for each plot
cell_types_Main <- c("B cells", "Fibroblast", "Macrophage", "Pericytes", "T cells", "VSMC", "Lymphatic Endothelial", "Mast cells", 'Neuronal cells', 'Plasma cells', 'Endothelial_1', 'Endothelial_2')
cell_types_Sub1 <- c('EC, Adventitia_A', 'EC, Adventitia_B', 'EC, Lumen', 'Fibroblast', 'Fibroblast (COL9A3+)', 'Fibromyocytes', 'Myofibroblast', 'Pericytes (APOE+)', 'Pericytes (MYH11+)', 'VSMC')
cell_types_Sub2 <- c('Macrophages, Resident', 'Macrophages, Foamy', 'Macrophages (SPP1+)', 'Macrophages, M4', 'Monocyte-derived DC')
cell_types_Sub3 <- c('B cells', 'Plasma cells', 'Natural Killer (XCL1+)', 'T cells (CM, CD4+)', 'T cells, Proliferating', 'Natural Killer (GNLY+)', 'T cells (EM, CD4+)', 'T cells (CD8+)')

# Generate the plots
plot_Main_Cell <- generate_plot(visium, cell_types_Main)
plot_Mural <- generate_plot(visium_sub, cell_types_Sub1)
plot_Macrophages <- generate_plot(visium_sub, cell_types_Sub2)
plot_Lympocytes <- generate_plot(visium_sub, cell_types_Sub3)


df <- visium %>% select(spot_id, Disease, sample, lumina_to_periphery, `B cells`:`Neuronal cells`, `Plasma cells`, `Endothelial_1`, `Endothelial_2`)
new_df <- df %>% filter(!sample =="FW106014") %>% filter(!Disease =="Fibroatheroma") %>%
  pivot_longer(cols = c(`B cells`:`Endothelial_2`), names_to = "Name", values_to = "Value") %>% filter(Value >0.6)
new_df <- new_df %>%
  group_by(Name) %>%
  mutate(scaled_value = scale(Value, center = FALSE, scale = max(Value) - min(Value))) %>%
  ungroup()

p1 <- ggplot(new_df, aes(y=Name, x=lumina_to_periphery, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 0.9,) +
  scale_fill_distiller(palette = "YlOrRd")+theme_classic()+theme(
    legend.position="none",
    axis.text = element_text(size=12, color = "black"), axis.title = element_text(size=12))+labs(x="Lumina to Periphery", y="")

new_df$Disease <- factor(new_df$Disease, levels = c("Control", "Atheroma", "Intermediate\nlesion"))
p2 <- ggplot(new_df, aes(y=Name, x=lumina_to_periphery, color=Disease, size=scaled_value))+geom_point(position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.7))+
  theme_classic()+theme(axis.text = element_text(size=12, color = "black"), axis.title = element_text(size=12))+labs(x="Lumina to Periphery", y="", size="Abundance")+
  scale_color_manual(values = c('#52006A','#CD113B', '#FF7600'))
ggarrange(p1,p2, widths = c(1,2))

