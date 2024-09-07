library(Seurat)
library(ggplot2)

# Set working directory
results_folder <- "./results"  
seurat_object_path <- "./data/seurat_object.Rds"

# Load Seurat object
seu.obj <- readRDS(seurat_object_path)

# Update metadata
seu.obj@meta.data$nb_clus_4[which(seu.obj@meta.data$nb_clus_4 == "monocyte")] <- "neutrophil"
seu.obj@meta.data$nb_clus_4[which(seu.obj@meta.data$nb_clus_4 == "e")] <- "B.cell.linage"

# Normalize data
seu.obj <- NormalizeData(seu.obj)

# Scale data
seu.obj <- ScaleData(seu.obj)

# Identify variable features
seu.obj <- FindVariableFeatures(seu.obj)

# Limit the number of variable features (using the top 26 variable features here)
top_features <- head(VariableFeatures(seu.obj), 26)

# Specify clusters to use (using all clusters here)
clusters_to_use <- unique(seu.obj@meta.data$nb_clus_4)

# Create a heatmap of feature expression
DoHeatmap(seu.obj, features = top_features, group.by = "nb_clus_4", size = 8, angle = 90) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu"))) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_text(size = 12),  
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),  # Set panel background to white
        plot.background = element_rect(fill = "white", color = NA),   # Set plot background to white
        strip.background = element_rect(fill = "white", color = NA)) +  # Set strip background to white
  guides(color = "none", fill = guide_colorbar(title = "Expression"))  # Remove identity legend and show color scale legend

# Save plot
ggsave(filename = paste0(results_folder, "/heatmap_single_cell_feature_expression.png"), 
       plot = last_plot(), width = 17, height = 12)

DoHeatmap(seu.obj, features = top_features, group.by = "nb_clus_4") +
  guides(color = "none") +
  theme(axis.text.y = element_text(size = 12))

ggsave(filename = paste0(results_folder, "/heatmap_single_cell_feature_expression1.png"), 
       plot = last_plot(), width = 15, height = 10)

# Create a clustered heatmap

# Reorder by clustering
seu.obj <- ScaleData(seu.obj, features = top_features)
seu.obj <- RunPCA(seu.obj, features = top_features)
seu.obj <- FindNeighbors(seu.obj, dims = 1:10)
seu.obj <- FindClusters(seu.obj, resolution = 0.5)

# Extract reordered data for heatmap
heatmap_data <- seu.obj@assays$RNA@scale.data[top_features,]
ordered_genes <- hclust(dist(heatmap_data))$order
ordered_cells <- hclust(dist(t(heatmap_data)))$order

# Create a heatmap with reordered genes and cells
DoHeatmap(seu.obj, features = top_features[ordered_genes], cells = Cells(seu.obj)[ordered_cells], group.by = "nb_clus_4", size = 8, angle = 90) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1),  # Set text size for cell names to 8 and angle to vertical
        axis.text.y = element_text(size = 15),  # Set text size for gene names to 15
        axis.ticks.x = element_blank()) +  # Remove x-axis ticks
  guides(fill = guide_colorbar(title = "Expression"))  # Remove identity legend and show color scale legend

# Save clustered heatmap
ggsave(filename = paste0(results_folder, "/heatmap_single_cell_feature_expression2.png"), 
       plot = last_plot(), width = 15, height = 12)
