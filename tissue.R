library(Seurat)
library(ggplot2)

plan("multisession", workers = 10)

DiscretePalette(32, palette = NULL, shuffle = FALSE)

# set working directory
results_folder = '/R/nanostring7/work4'

# set object
seu.obj <- readRDS("/R/nanostring7/6 Analysis/Data objects/seurat_object.Rds")


# change the name of "monocyte" in nb_clus_4 to "neutrophil".
seu.obj@meta.data$nb_clus_4[which(seu.obj@meta.data$nb_clus_4 == "monocyte")] <- "neutrophil"
seu.obj@meta.data$nb_clus_4[which(seu.obj@meta.data$nb_clus_4 == "e")] <- "B.cell.linage"

# set idents to nb_clus_4
Idents(seu.obj) <- seu.obj$nb_clus_4


# Continuous processing: selecting FOV 1 Run5609_18_72_TB Background White
numbers <- c(1,2)
for (i in numbers)  {
  
  # Select FOV
  use_fov <- i # FOV desired
  use_slide_image <- 'Run5609.18.72.TB' # Slide desired, as named in images
  use_slide_metadata <- 'Run5609_18_72_TB' # Slide desired, as named in the metadata column ‘Run_Tissue_name’
  
  # First get cells in your FOV
  cells_of_interest <- seu.obj$id[(seu.obj$fov == use_fov) & (seu.obj$Run_Tissue_name == use_slide_metadata)]
  # Then find spatial boundaries of the rectangle containing the centroids of these cells
  centroid_data <- seu.obj@images[[use_slide_image]]$centroids
  zoom_fov <- apply(centroid_data@coords[centroid_data@cells %in% cells_of_interest,], 2, range)
  
  # Cell types + gene expression (SPP1, MARCO)
  ImageDimPlot(seu.obj,
               fov = use_slide_image,
               cols = "glasbey",
               alpha = 0.4,
               border.color = "gray",
               molecules = c("SPP1","MARCO"),
               mols.cols = c("green","red"),
               mols.size = 0.5,
               nmols = 10000,
               size = 0.01,
               border.size = 0.2)+ 
    xlim(zoom_fov[, 2]) +
    ylim(zoom_fov[, 1])+
    theme(panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"),
          legend.background = element_rect(fill = "white"),
          legend.key = element_rect(fill = "white"),
          legend.text = element_text(color = "black", face = "bold"),
          legend.title = element_text(color = "black", face = "bold"))
  filename <- paste0("18_72_TBfov", i, "c_white.png")
  ggsave(filename, bg="white", path = results_folder, dpi = 300, width = 20, height = 12)
  
}

# Continuous processing: selecting FOV 1 Run5609_18_72_TB
numbers <- 1:25
for (i in numbers)  {
  
  # Select FOV
  use_fov <- i # FOV desired
  use_slide_image <- 'Run5609.18.72.TB' # Slide desired, as named in images
  use_slide_metadata <- 'Run5609_18_72_TB' # Slide desired, as named in the metadata column ‘Run_Tissue_name’
  
  # First get cells in your FOV
  cells_of_interest <- seu.obj$id[(seu.obj$fov == use_fov) & (seu.obj$Run_Tissue_name == use_slide_metadata)]
  # Then find spatial boundaries of the rectangle containing the centroids of these cells
  centroid_data <- seu.obj@images[[use_slide_image]]$centroids
  zoom_fov <- apply(centroid_data@coords[centroid_data@cells %in% cells_of_interest,], 2, range)
  
  # Cell types + gene expression (SPP1, MARCO)
  ImageDimPlot(seu.obj,
               fov = use_slide_image,
               cols = "glasbey",
               alpha = 0.4,
               border.color = "gray",
               molecules = c("SPP1","MARCO"),
               mols.cols = c("green","red"),
               mols.size = 0.5,
               nmols = 10000,
               size = 0.01,
               border.size = 0.2)+ 
    xlim(zoom_fov[, 2]) +
    ylim(zoom_fov[, 1])
  filename <- paste0("18_72_TBfov", i, "e.png")
  ggsave(filename, path = results_folder, dpi = 300, width = 20, height = 12)
  
}

# Continuous processing: selecting FOV 2 Run5612_18_76_NTM_MAC_pos
numbers <- c(1,2,3,4,6,7,8,9,10,11,14,15,16,17,19,20,21,22,23,24,25)
for (i in numbers)  {
  
  # Select FOV
  use_fov <- i # FOV desired
  use_slide_image <- 'Run5612.18.76.NTM.MAC.pos' # Slide desired, as named in images
  use_slide_metadata <- 'Run5612_18_76_NTM_MAC_pos' # Slide desired, as named in the metadata column ‘Run_Tissue_name’
  
  # First get cells in your FOV
  cells_of_interest <- seu.obj$id[(seu.obj$fov == use_fov) & (seu.obj$Run_Tissue_name == use_slide_metadata)]
  # Then find spatial boundaries of the rectangle containing the centroids of these cells
  centroid_data <- seu.obj@images[[use_slide_image]]$centroids
  zoom_fov <- apply(centroid_data@coords[centroid_data@cells %in% cells_of_interest,], 2, range)
  
  # Cell types + gene expression (SPP1, MARCO)
  ImageDimPlot(seu.obj,
               fov = use_slide_image,
               cols = "glasbey",
               alpha = 0.4,
               border.color = "gray",
               molecules = c("SPP1","MARCO"),
               mols.cols = c("green","red"),
               mols.size = 0.5,
               nmols = 10000,
               size = 0.01,
               border.size = 0.2)+ 
    xlim(zoom_fov[, 2]) +
    ylim(zoom_fov[, 1])
  filename <- paste0("18_76_NTMfov", i, "e.png")
  ggsave(filename, path = results_folder, dpi = 300, width = 20, height = 12)
  
}

# Continuous processing: selecting FOV 3 Run5688_18_77_Bronchiectasis_MAC_neg
numbers <- c(1,2,3,4,6,7,9,17,19,20,21,22,24,25)
for (i in numbers)  {
  
  # Select FOV
  use_fov <- i # FOV desired
  use_slide_image <- 'Run5688.18.77.Bronchiectasis.MAC.neg' # Slide desired, as named in images
  use_slide_metadata <- 'Run5688_18_77_Bronchiectasis_MAC_neg' # Slide desired, as named in the metadata column ‘Run_Tissue_name’
  
  # First get cells in your FOV
  cells_of_interest <- seu.obj$id[(seu.obj$fov == use_fov) & (seu.obj$Run_Tissue_name == use_slide_metadata)]
  # Then find spatial boundaries of the rectangle containing the centroids of these cells
  centroid_data <- seu.obj@images[[use_slide_image]]$centroids
  zoom_fov <- apply(centroid_data@coords[centroid_data@cells %in% cells_of_interest,], 2, range)
  
  # Cell types + gene expression (SPP1, MARCO)
  ImageDimPlot(seu.obj,
               fov = use_slide_image,
               cols = "glasbey",
               alpha = 0.4,
               border.color = "gray",
               molecules = c("SPP1","MARCO"),
               mols.cols = c("green","red"),
               mols.size = 0.5,
               nmols = 10000,
               size = 0.01,
               border.size = 0.2)+ 
    xlim(zoom_fov[, 2]) +
    ylim(zoom_fov[, 1])
  filename <- paste0("18_77_Bronchiectasisfov", i, "e.png")
  ggsave(filename, path = results_folder, dpi = 300, width = 20, height = 12)
  
}

# Continuous processing: selecting FOV 4 NTM_MAC_18-80 (SPP1, MARCO)
numbers <- c(1,2,3,5,6,7,8,10,12,13,14,15,16,17,18,19,20,21,22,23,24,25)
for (i in numbers)  {
  
  # Select FOV
  use_fov <- i # FOV desired
  use_slide_image <- 'Run5689.18.80.NTM.MAC' # Slide desired, as named in images
  use_slide_metadata <- 'Run5689_18_80_NTM_MAC' # Slide desired, as named in the metadata column ‘Run_Tissue_name’
  
  # First get cells in your FOV
  cells_of_interest <- seu.obj$id[(seu.obj$fov == use_fov) & (seu.obj$Run_Tissue_name == use_slide_metadata)]
  # Then find spatial boundaries of the rectangle containing the centroids of these cells
  centroid_data <- seu.obj@images[[use_slide_image]]$centroids
  zoom_fov <- apply(centroid_data@coords[centroid_data@cells %in% cells_of_interest,], 2, range)
  
  # Cell types + gene expression (SPP1, MARCO)
  ImageDimPlot(seu.obj,
               fov = use_slide_image,
               cols = "glasbey",
               alpha = 0.4,
               border.color = "gray",
               molecules = c("SPP1","MARCO"),
               mols.cols = c("green","red"),
               mols.size = 0.5,
               nmols = 10000,
               size = 0.01,
               border.size = 0.2)+ 
    xlim(zoom_fov[, 2]) +
    ylim(zoom_fov[, 1])
  filename <- paste0("18_80_NTM_MACfov", i, "e.png")
  ggsave(filename, path = results_folder, dpi = 300, width = 20, height = 12)
  
}


# Continuous processing: selecting FOV 5 R6021_Lung_TB_18-78 (SPP1, MARCO)
numbers <- 1:40
for (i in numbers)  {
  
  # Select FOV.
  use_fov <- i # FOV desired
  use_slide_image <- 'R6021.Lung.TB.18.78' # Slide desired, as named in images
  use_slide_metadata <- 'R6021_Lung_TB_18-78' # Slide desired, as named in the metadata column ‘Run_Tissue_name’
  
  # First get cells in your FOV
  cells_of_interest <- seu.obj$id[(seu.obj$fov == use_fov) & (seu.obj$Run_Tissue_name == use_slide_metadata)]
  # Then find spatial boundaries of the rectangle containing the centroids of these cells
  centroid_data <- seu.obj@images[[use_slide_image]]$centroids
  zoom_fov <- apply(centroid_data@coords[centroid_data@cells %in% cells_of_interest,], 2, range)
  
  # Cell types + gene expression (SPP1, MARCO)
  ImageDimPlot(seu.obj,
               fov = use_slide_image,
               cols = "glasbey",
               alpha = 0.4,
               border.color = "gray",
               molecules = c("SPP1","MARCO"),
               mols.cols = c("green","red"),
               mols.size = 0.5,
               nmols = 10000,
               size = 0.01,
               border.size = 0.2)+ 
    xlim(zoom_fov[, 2]) +
    ylim(zoom_fov[, 1])
  filename <- paste0("18_78_TBfov", i, "e.png")
  ggsave(filename, path = results_folder, dpi = 300, width = 10, height = 6)
  
}

# Continuous processing: selecting FOV 6 R6021_Lung_TB_18-79
numbers <- 1:40
for (i in numbers)  {
  
  # Select FOV.
  use_fov <- i # FOV desired
  use_slide_image <- 'R6021.Lung.TB.18.79' # Slide desired, as named in images
  use_slide_metadata <- 'R6021_Lung_TB_18-79' # Slide desired, as named in the metadata column ‘Run_Tissue_name’
  
  # First get cells in your FOV
  cells_of_interest <- seu.obj$id[(seu.obj$fov == use_fov) & (seu.obj$Run_Tissue_name == use_slide_metadata)]
  # Then find spatial boundaries of the rectangle containing the centroids of these cells
  centroid_data <- seu.obj@images[[use_slide_image]]$centroids
  zoom_fov <- apply(centroid_data@coords[centroid_data@cells %in% cells_of_interest,], 2, range)
  
  # Cell types + gene expression (SPP1, MARCO)
  ImageDimPlot(seu.obj,
               fov = use_slide_image,
               cols = "glasbey",
               alpha = 0.4,
               border.color = "gray",
               molecules = c("SPP1","MARCO"),
               mols.cols = c("green","red"),
               mols.size = 0.5,
               nmols = 10000,
               size = 0.01,
               border.size = 0.2)+ 
    xlim(zoom_fov[, 2]) +
    ylim(zoom_fov[, 1])
  filename <- paste0("18_79_TBfov", i, "e.png")
  ggsave(filename, path = results_folder, dpi = 300, width = 10, height = 6)
  
}

# Continuous processing: selecting FOV 7 NTM_MAC_18-81
numbers <- 1:40
for (i in numbers)  {
  
  # Select FOV.
  use_fov <- i # FOV desired
  use_slide_image <- 'R6021.Lung.NTM.MAC.18.81' # Slide desired, as named in images
  use_slide_metadata <- 'R6021_Lung_NTM_MAC_18-81' # Slide desired, as named in the metadata column ‘Run_Tissue_name’
  
  # First get cells in your FOV
  cells_of_interest <- seu.obj$id[(seu.obj$fov == use_fov) & (seu.obj$Run_Tissue_name == use_slide_metadata)]
  # Then find spatial boundaries of the rectangle containing the centroids of these cells
  centroid_data <- seu.obj@images[[use_slide_image]]$centroids
  zoom_fov <- apply(centroid_data@coords[centroid_data@cells %in% cells_of_interest,], 2, range)
  
  # Cell types + gene expression (SPP1, MARCO)
  ImageDimPlot(seu.obj,
               fov = use_slide_image,
               cols = "glasbey",
               alpha = 0.4,
               border.color = "gray",
               molecules = c("SPP1","MARCO"),
               mols.cols = c("green","red"),
               mols.size = 0.5,
               nmols = 10000,
               size = 0.01,
               border.size = 0.2)+  
    xlim(zoom_fov[, 2]) +
    ylim(zoom_fov[, 1])
  filename <- paste0("18_81_MACfov", i, "e.png")
  ggsave(filename, path = results_folder, dpi = 300, width = 10, height = 6)
  
}
