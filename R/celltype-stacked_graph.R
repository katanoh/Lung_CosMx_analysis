library(Seurat)
library(ggplot2)
library(viridis)

plan("multisession", workers = 10)
DiscretePalette(32, palette = NULL, shuffle = FALSE)

# Set working directory
results_folder <- './results'

# Load the Seurat object
seu.obj <- readRDS("./data/seurat_object.Rds")

# Change the name of "e" in nb_clus_4 to "B.cell.lineage"
seu.obj@meta.data$nb_clus_4[which(seu.obj@meta.data$nb_clus_4 == "e")] <- "B.cell.lineage"
# Change the name of "monocyte" in nb_clus_4 to "neutrophil"
seu.obj@meta.data$nb_clus_4[which(seu.obj@meta.data$nb_clus_4 == "monocyte")] <- "neutrophil"

# Set idents to nb_clus_4
Idents(seu.obj) <- seu.obj$nb_clus_4

# Rename tissue values and set order
tissue_map <- c(
  "18_76_NTM_MAC_pos" = "NTM1",
  "18_80_NTM_MAC_pos" = "NTM2",
  "Lung_NTM_MAC_18-81" = "NTM3",
  "18_77_Bronchiectasis_MAC_neg" = "NTM4",
  "18_72_TB" = "MTB1", 
  "Lung_TB_18-78" = "MTB2", 
  "Lung_TB_18-79" = "MTB3"
)

# Apply tissue map and set factor levels to ensure order
seu.obj@meta.data$tissue <- factor(seu.obj@meta.data$tissue, levels = names(tissue_map), labels = tissue_map)

# Create a color palette based on the unique values of nb_clus_4
colors <- rainbow(length(unique(seu.obj@meta.data$nb_clus_4)))

# Shuffle the colors randomly
shuffled_colors <- sample(colors)

# Create a bar plot
p <- ggplot(seu.obj@meta.data, aes(x = nb_clus_4)) + 
  geom_bar(aes(fill = nb_clus_4, y = after_stat(count)), color = "black") +
  geom_text(aes(label = after_stat(count), y = after_stat(count)), stat = "count", 
            position = position_dodge(0.9), vjust = -0.3) +
  scale_fill_manual(values = shuffled_colors) +
  facet_wrap(~tissue, scales = "free_y", ncol = 4) +  # Facet by tissue with new labels and order, 4 columns
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14, 
                                   face = "bold", color = "black"),  # Adjust x-axis labels
        legend.position = "none",
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.line = element_line(color = "black"),
        strip.text = element_text(size = 18,  face = "bold",margin = margin(b = 10)),
        strip.background = element_rect(fill = "white"))

# Save the bar plot to the specified folder
ggsave(filename = paste0(results_folder, "/bar1a.png"), plot = p, width = 20, height = 10, bg = "white")


# Reorder tissue levels
order_levels <- c("18_76_NTM_MAC_pos","18_80_NTM_MAC_pos",  
                  "Lung_NTM_MAC_18-81","18_77_Bronchiectasis_MAC_neg",
                  "18_72_TB", "Lung_TB_18-78","Lung_TB_18-79")
seu.obj@meta.data$tissue <- factor(seu.obj@meta.data$tissue, levels = order_levels)

# Create a color palette based on the unique values of nb_clus_4
colors <- rainbow(length(unique(seu.obj@meta.data$nb_clus_4)))

# Shuffle the colors randomly
shuffled_colors <- sample(colors)

# Create a bar plot
p <- ggplot(seu.obj@meta.data, aes(x = nb_clus_4)) + 
  geom_bar(aes(fill = nb_clus_4, y = after_stat(count)), color = "black") +
  geom_text(aes(label = after_stat(count), y = after_stat(count)), 
            stat = "count", position = position_dodge(0.9), vjust = -0.3,size = 3) +
  scale_fill_manual(values = shuffled_colors) +
  facet_wrap(~tissue, scales = "free_y", ncol = 4) +  # Facet by tissue instead of slide_ID_numeric
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, 
                                   size = 10, face = "bold", color = "black"),  # Adjust x-axis labels
        legend.position = "none",
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.line = element_line(color = "black"),
        strip.text = element_text(size = 18, face = "bold",margin = margin(b = 10)),
        strip.background = element_rect(fill = "white"))

# Save the bar plot to the specified folder
ggsave(filename = paste0(results_folder, "/bar2.png"), plot = p, width = 15, height = 9, bg = "white")


library(tidyr)
library(dplyr)

# Group by tissue and nb_clus_4 to get counts
cell_counts <- seu.obj@meta.data %>%
  group_by(tissue, nb_clus_4) %>%
  summarise(count = n()) %>%
  ungroup()

# Calculate the total counts for each tissue
total_counts <- cell_counts %>%
  group_by(tissue) %>%
  summarise(total = sum(count))

# Join the total counts to the cell counts and calculate percentages
cell_counts <- cell_counts %>%
  left_join(total_counts, by = "tissue") %>%
  mutate(percentage = (count / total) * 100) %>%
  select(tissue, nb_clus_4, percentage)

# Spread the data to have tissues as columns and clusters as rows
cell_percentage_table <- cell_counts %>%
  spread(key = tissue, value = percentage, fill = 0)

# Calculate overall percentages
overall_counts <- seu.obj@meta.data %>%
  group_by(nb_clus_4) %>%
  summarise(count = n())

total_count <- sum(overall_counts$count)

overall_percentages <- overall_counts %>%
  mutate(overall_percentage = (count / total_count) * 100) %>%
  select(nb_clus_4, overall_percentage)

# Join the overall percentages to the table
cell_percentage_table <- cell_percentage_table %>%
  left_join(overall_percentages, by = c("nb_clus_4" = "nb_clus_4"))

# Write to CSV
write.csv(cell_percentage_table, file = paste0(results_folder, "/cluster_tissue_percentages_with_overall.csv"), row.names = FALSE)

# Group by tissue and nb_clus_4 to get counts
cell_counts <- seu.obj@meta.data %>%
  group_by(tissue, nb_clus_4) %>%
  summarise(count = n()) %>%
  ungroup()

# Spread the data to have tissues as columns and clusters as rows
cell_count_table <- cell_counts %>%
  spread(key = tissue, value = count, fill = 0)

# Calculate overall counts
overall_counts <- seu.obj@meta.data %>%
  group_by(nb_clus_4) %>%
  summarise(overall_count = n())

# Join the overall counts to the table
cell_count_table <- cell_count_table %>%
  left_join(overall_counts, by = c("nb_clus_4" = "nb_clus_4"))

# Write to CSV
write.csv(cell_count_table, file = paste0(results_folder, "/cluster_tissue_counts_with_overall.csv"), row.names = FALSE)



# Create a stacked bar plot
library(dplyr)

# Summarize the data
df_summary <- seu.obj@meta.data %>%
  group_by(tissue, nb_clus_4) %>%
  summarize(count = n()) %>%
  mutate(prop = count / sum(count))

# Create the plot
p <- ggplot(seu.obj@meta.data, aes(x = nb_clus_4)) + 
  geom_bar(aes(fill = nb_clus_4, y = ..count..), color = "black") +
  geom_text(aes(label = ..count.., y = ..count..), stat = "count", position = position_dodge(0.9), vjust = -0.3) +
  scale_fill_manual(values = shuffled_colors) +
  facet_wrap(~tissue, scales = "free_y", ncol = 4) +  # Adjust ncol parameter to 4
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14, face = "bold", color = "black"),  # Adjust x-axis labels
        legend.position = "none",
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.line = element_line(color = "black"),
        strip.text = element_text(size = 16, margin = margin(b = 10)),
        strip.background = element_rect(fill = "white"))
p
# Save the bar plot to the specified folder
ggsave(filename = paste0(results_folder, "/bar1.png"), plot = p, width = 25, height = 11, bg = "white")


# Data aggregation
df_summary <- seu.obj@meta.data %>%
  filter(tissue == "18_72_TB") %>%  # Select only where tissue is "18_72_TB"
  group_by(tissue, nb_clus_4, fov) %>%
  summarize(count = n()) %>%
  mutate(prop = count / sum(count))

# Plot drawing
p <- ggplot(df_summary, aes(x = fov, y = prop, fill = nb_clus_4)) +  
  geom_bar(stat = "identity", position = "fill", color = "black") +
  scale_fill_manual(values = shuffled_colors, guide = guide_legend(ncol = 1, title = NULL)) +  
  scale_x_discrete(name = "FOV", limits = unique(df_summary$fov)) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, face = "bold", color = "black"),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.line = element_line(color = "black"),
        strip.text = element_text(size = 16, margin = margin(b = 10)),
        strip.background = element_rect(fill = "white"))
p

# Save the plot to the specified folder
ggsave(filename = paste0(results_folder, "/stacked_bar1_18_72_TB_per_fov_horizontal.png"),
       plot = p, width = 8, height = 9, bg = "white")



library(gridExtra)

# 1. Get all unique tissue values
unique_tissues <- unique(seu.obj@meta.data$tissue)

# 2. Create and save plots for each tissue in a list
plot_list <- list()
for (current_tissue in unique_tissues) {
  df_summary <- seu.obj@meta.data %>%
    filter(tissue == current_tissue) %>%
    group_by(tissue, nb_clus_4, fov) %>%
    summarize(count = n()) %>%
    mutate(prop = count / sum(count))
  
  p <- ggplot(df_summary, aes(x = fov, y = prop, fill = nb_clus_4)) +
    geom_bar(stat = "identity", position = "fill", color = "black") +
    scale_fill_manual(values = shuffled_colors, guide = guide_legend(ncol = 1, title = NULL)) +
    scale_x_discrete(name = "Fov", limits = unique(df_summary$fov)) +  # <-- Modification here
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, face = "bold", color = "black"),
          panel.background = element_blank(),
          plot.background = element_blank(),
          axis.line = element_line(color = "black"),
          strip.text = element_text(size = 16, margin = margin(b = 10)),
          strip.background = element_rect(fill = "white"),
          plot.title = element_text(face = "bold"))
  
  plot_list[[current_tissue]] <- p
}

# Render the plots
grid.arrange(grobs = plot_list, ncol = 2)

library(cowplot)

# Add tissue name as title for each plot
for (current_tissue in unique_tissues) {
  plot_list[[current_tissue]] <- plot_list[[current_tissue]] + labs(title = current_tissue)
}

# Combine the plots with legend on the left
combined_plot <- plot_grid(plotlist = plot_list, ncol = 2, align = 'hv', rel_widths = c(1, 1, 1, 1))

# Save the combined plot
file_name <- paste0(results_folder, "/combined_stacked_bar_per_fov_horizontal_with_titles.png")
ggsave(file_name, combined_plot, width = 25, height = 32, bg = "white")


#Stacked bar graph with macrophage highlighted.

# Get all unique values of nb_clus_4 across all tissues
all_clusters <- unique(seu.obj@meta.data$nb_clus_4)

# Assign colors from shuffled_colors for each of these values
complete_colors <- shuffled_colors[all_clusters]

# For missing values, assign random colors
missing_colors <- setdiff(all_clusters, names(shuffled_colors))
complete_colors[missing_colors] <- scales::hue_pal()(length(missing_colors))

# Highlighted clusters
highlighted_clusters <- c("macrophage_a", "macrophage_b","macrophage_c","macrophage_d")

# Modify colors to make non-highlighted clusters translucent
modified_colors_vec <- sapply(names(complete_colors), function(cluster_name) {
  if (cluster_name %in% highlighted_clusters) {
    return(complete_colors[cluster_name])
  } else {
    return(scales::alpha(complete_colors[cluster_name], alpha = 0.1))
  }
})

# Convert to named vector
modified_colors <- setNames(as.vector(modified_colors_vec), names(complete_colors))

# Check the modified colors
print(modified_colors)



# Recreate the plots with the modified colors
plot_list <- list()
for (current_tissue in unique_tissues) {
  df_summary <- seu.obj@meta.data %>%
    filter(tissue == current_tissue) %>%
    group_by(tissue, nb_clus_4, fov) %>%
    summarize(count = n()) %>%
    mutate(prop = count / sum(count))
  
  p <- ggplot(df_summary, aes(x = fov, y = prop, fill = nb_clus_4)) +
    geom_bar(stat = "identity", position = "fill", color = "black") +
    scale_fill_manual(values = modified_colors, guide = guide_legend(ncol = 1, title = NULL)) +
    scale_x_discrete(name = "Fov", limits = unique(df_summary$fov)) +
    labs(title = current_tissue) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, face = "bold", color = "black"),
          panel.background = element_blank(),
          plot.background = element_blank(),
          axis.line = element_line(color = "black"),
          strip.text = element_text(size = 16, margin = margin(b = 10)),
          strip.background = element_rect(fill = "white"),
          plot.title = element_text(face = "bold")) 
  
  plot_list[[current_tissue]] <- p
}

# Combine the plots with legend on the left
combined_plot <- plot_grid(plotlist = plot_list, ncol = 2, align = 'hv', rel_widths = c(1, 1, 1, 1))

# Save the combined plot
file_name <- paste0(results_folder, "/combined_stacked_bar_per_fov_horizontal_with_titles_translucent.png")
ggsave(file_name, combined_plot, width = 25, height = 32, bg = "white")


#Stacked bar graph with T and macrophage highlighted.
# Get all unique values of nb_clus_4 across all tissues
all_clusters <- unique(seu.obj@meta.data$nb_clus_4)

# Assign colors from shuffled_colors for each of these values
complete_colors <- shuffled_colors[all_clusters]

# For missing values, assign random colors
missing_colors <- setdiff(all_clusters, names(shuffled_colors))
complete_colors[missing_colors] <- scales::hue_pal()(length(missing_colors))

# Highlighted clusters
highlighted_clusters <- c("CD4+.T.cell","CD8+.cytotoxic.T.cell",
                          "macrophage_a", "macrophage_b","macrophage_c","macrophage_d", 
                          "B.cell","B.cell.linage","plasma.cell")

# Modify colors to make non-highlighted clusters translucent
modified_colors_vec <- sapply(names(complete_colors), function(cluster_name) {
  if (cluster_name %in% highlighted_clusters) {
    return(complete_colors[cluster_name])
  } else {
    return(scales::alpha(complete_colors[cluster_name], alpha = 0.1))
  }
})

# Convert to named vector
modified_colors <- setNames(as.vector(modified_colors_vec), names(complete_colors))

# Check the modified colors
print(modified_colors)



# Recreate the plots with the modified colors
plot_list <- list()
for (current_tissue in unique_tissues) {
  df_summary <- seu.obj@meta.data %>%
    filter(tissue == current_tissue) %>%
    group_by(tissue, nb_clus_4, fov) %>%
    summarize(count = n()) %>%
    mutate(prop = count / sum(count))
  
  p <- ggplot(df_summary, aes(x = fov, y = prop, fill = nb_clus_4)) +
    geom_bar(stat = "identity", position = "fill", color = "black") +
    scale_fill_manual(values = modified_colors, guide = guide_legend(ncol = 1, title = NULL)) +
    scale_x_discrete(name = "Fov", limits = unique(df_summary$fov)) +
    labs(title = current_tissue) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, face = "bold", color = "black"),
          panel.background = element_blank(),
          plot.background = element_blank(),
          axis.line = element_line(color = "black"),
          strip.text = element_text(size = 16, margin = margin(b = 10)),
          strip.background = element_rect(fill = "white"),
          plot.title = element_text(face = "bold")) 
  
  plot_list[[current_tissue]] <- p
}

# Combine the plots with legend on the left
combined_plot <- plot_grid(plotlist = plot_list, ncol = 2, align = 'hv', rel_widths = c(1, 1, 1, 1))

# Save the combined plot
file_name <- paste0(results_folder, "/combined_stacked_bar_per_fov_horizontal_with_titles_translucentTcell.png")
ggsave(file_name, combined_plot, width = 25, height = 32, bg = "white")




# Load necessary libraries
library(ggplot2)
library(gridExtra)

# Function to create dot plot
create_dot_plot <- function(current_tissue, modified_colors, seu.obj) {
  df_summary <- seu.obj@meta.data %>%
    filter(tissue == current_tissue) %>%
    group_by(fov, nb_clus_4) %>%
    summarize(count = n()) %>%
    arrange(desc(count))
  
  p <- ggplot(df_summary, aes(x = reorder(fov, -count), y = count, color = nb_clus_4)) +
    geom_point(size = 3) +
    scale_color_manual(values = modified_colors) +
    labs(title = current_tissue, x = "Fov", y = "Cell Count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, face = "bold", color = "black"),
          axis.line = element_line(color = "black"),
          plot.title = element_text(face = "bold"))
  
  return(p)
}

# Create dot plots for each tissue
dot_plot_list <- list()
for (current_tissue in unique_tissues) {
  dot_plot <- create_dot_plot(current_tissue, modified_colors, seu.obj)
  dot_plot_list[[current_tissue]] <- dot_plot
}

# Combine dot plots
combined_dot_plot <- grid.arrange(grobs = dot_plot_list, ncol = 2)

# Save the combined dot plot
file_name_dot_plot <- "mac-CD8-dotplot.png"
ggsave(file_name_dot_plot, combined_dot_plot, width = 16, height = 16, bg = "white")
