library(Seurat)
library(ggplot2)
library(gridExtra)
library(future)
library(future.apply)
library(ggrepel)
library(patchwork)

plan("multisession", workers = 10)
DiscretePalette(32, palette = NULL, shuffle = FALSE)

# Set working directory
results_folder = './results'

seu.obj <- readRDS(".data/seurat_object.Rds")

# Change the name of "e" in nb_clus_4 to "B.cell.lineage".
seu.obj@meta.data$nb_clus_4[which(seu.obj@meta.data$nb_clus_4 == "e")] <- "B.cell.lineage"
# Change the name of "monocyte" in nb_clus_4 to "neutrophil".
seu.obj@meta.data$nb_clus_4[which(seu.obj@meta.data$nb_clus_4 == "monocyte")] <- "neutrophil"

# Set idents to nb_clus_4
Idents(seu.obj) <- seu.obj$nb_clus_4

# The following steps insert genes into the RNA slot
# Retrieve the layers from the Nanostring assay
nanostring_assay <- seu.obj[["Nanostring"]]
layer_names <- names(nanostring_assay@layers)

# Retrieve data for each layer and store in a list
layer_list <- lapply(layer_names, function(layer_name) {
  GetAssayData(nanostring_assay, layer = layer_name)
})

# Combine the layers along the rows
combined_counts <- do.call(cbind, layer_list)

# Display the first few rows of the combined data
print(dim(combined_counts))
print(head(combined_counts))

# Create temporary gene expression data
gene_names <- rownames(combined_counts)
cell_names <- colnames(combined_counts)

# Create a new assay with the temporary data
seu.obj[["RNA"]] <- CreateAssayObject(counts = combined_counts)

# Set the default assay to RNA
DefaultAssay(seu.obj) <- "RNA"

# Confirm that the data has been added successfully
print(head(GetAssayData(seu.obj, slot = "counts")))


# Start UMAP
DimPlot(seu.obj, group.by = "nb_clus_4", raster = FALSE, label = TRUE, label.size = 4)
ggsave(file = "UMAP1.png", path = results_folder, dpi = 300, width = 10, height = 5)

# Create a UMAP plot where labels do not overlap
# Create UMAP plot
umap_plot <- DimPlot(seu.obj, group.by = "nb_clus_4", raster = FALSE, label = FALSE) +
  ggtitle(NULL) +
  theme_minimal()

# Get UMAP plot data
umap_data <- umap_plot$data

# Calculate the median coordinates for each cluster
library(dplyr)
cluster_centers <- umap_data %>%
  group_by(nb_clus_4) %>%
  summarize(umap_1 = median(umap_1), umap_2 = median(umap_2))

# Add cluster name labels to the UMAP plot
umap_plot <- umap_plot + 
  geom_text_repel(data = cluster_centers, aes(x = umap_1, y = umap_2, label = nb_clus_4), 
                  size = 2, box.padding = 0.5, point.padding = 0.5, max.overlaps = Inf)

# Display the plot
print(umap_plot)

# Save the plot
ggsave(file = "UMAP1a.png", plot = umap_plot, path = results_folder, dpi = 300, width = 9, height = 5)
ggsave(file = "UMAP1b.png", plot = umap_plot, path = results_folder, dpi = 300, width = 8, height = 3.5)
ggsave(file = "UMAP1c.png", plot = umap_plot, path = results_folder, dpi = 400, width = 8, height = 4)

DimPlot(seu.obj, group.by = "nb_clus_4", raster = FALSE, label = FALSE) # Set label to FALSE
ggsave(filename = "UMAP2.png", path = results_folder, dpi = 300, width = 15, height = 10)


# Create FeaturePlot
FeaturePlot(seu.obj, features = c("SPP1", "FN1", "CD44", "MARCO"), 
            raster = FALSE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
ggsave(file = "UMAP_spp1_FN1_CD44_MARCO.png", path = results_folder,
       dpi = 300, width = 10, height = 10)

FeaturePlot(seu.obj, features = c("SPP1", "FN1", "CD44", "IL1RN"), 
            raster = FALSE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
ggsave(file = "UMAP_spp1_FN1_CD44_IL1RN.png", path = results_folder,
       dpi = 300, width = 10, height = 10)

FeaturePlot(seu.obj, features = c("SPP1", "CD68", "CD14", "MARCO"), 
            raster = FALSE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
ggsave(file = "UMAP_spp1_CD68_CD14_MARCO.png", path = results_folder,
       dpi = 300, width = 10, height = 10)


FeaturePlot(seu.obj, features = "SPP1", raster = FALSE, pt.size = 0.2, min.cutoff = 0, max.cutoff = 10)
ggsave(file = "UMAP_spp1.png", path = results_folder,
       dpi = 300, width = 10, height = 10)

FeaturePlot(seu.obj, features = "CD44", raster = FALSE, pt.size = 0.2, min.cutoff = 0, max.cutoff = 10)
ggsave(file = "UMAP_CD44.png", path = results_folder,
       dpi = 300, width = 10, height = 10)


markers1.to.plot <- c("CCL19", "SCGB3A1", "LAMP3", "IGHM", "IGHG2", "VWF", "IL7R", "CCL5",
                      "TUBB4B", "ITGAX", "IL1RN", "COL6A3", "CCL2", "COL3A1",
                      "CCL21", "FCER1G", "LYZ", "CD14", "SPP1", "TPSB2", "CPA3", "S100A9",
                      "MYH11", "GNLY", "IGHA1", "JCHAIN", "CTLA4")
# Create FeaturePlot and remove scale
plots <- lapply(markers1.to.plot, function(marker) {
  FeaturePlot(seu.obj, features = marker, raster = FALSE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5) + 
    theme(legend.position = "none")
})

# Arrange plots in 5 columns
combined_plot <- wrap_plots(plots, ncol = 6)

# Save the plot
ggsave(filename = "UMAP_markers1.png", plot = combined_plot, 
       path = results_folder, dpi = 300, width = 15, height = 10)


# UMAP by sample
# Select FOVs for batch processing
sample_list <- c("NTM1", "NTM2", "NTM3", "NTM4", "MTB1", "MTB2", "MTB3")
slide <- c("18_76_NTM_MAC_pos", "18_80_NTM_MAC_pos", 
           "Lung_NTM_MAC_18-81", "18_77_Bronchiectasis_MAC_neg",
           "18_72_TB", "Lung_TB_18-78", "Lung_TB_18-79")

# Create a named vector to map `tissue` to `Sample-ID`
sample_mapping <- c("18_76_NTM_MAC_pos" = "NTM1",
                    "18_80_NTM_MAC_pos" = "NTM2",  "Lung_NTM_MAC_18-81" = "NTM3",
                    "18_77_Bronchiectasis_MAC_neg" = "NTM4", "18_72_TB" = "MTB1", 
                    "Lung_TB_18-78" = "MTB2", 
                    "Lung_TB_18-79" = "MTB3")

plot_list <- list()

for (i in sample_list) {
  # Find the corresponding `tissue` value
  tissue <- names(sample_mapping[sample_mapping == i])
  
  # Create a data frame with UMAP coordinates and subset information
  umap_data <- as.data.frame(seu.obj@reductions$umap@cell.embeddings)
  umap_data$subset <- seu.obj@meta.data$tissue == tissue
  
  # Create UMAP plot with subset cells highlighted in blue
  plot_umap <- ggplot() +
    geom_point(data = umap_data, aes(x = UMAP_1, y = UMAP_2, color = subset), size = 0.1) +
    scale_color_manual(values = c("gray", "blue")) +
    ggtitle(i) +
    theme_minimal()
  
  plot_list[[i]] <- plot_umap
}

# Arrange plots in 3 columns
combined_plot <- grid.arrange(grobs = plot_list, ncol = 3)

# Save the plot
ggsave(filename = "UMAP_by_sample_v2.png", plot = combined_plot, 
       path = results_folder, dpi = 300, width = 20, height = 15)
