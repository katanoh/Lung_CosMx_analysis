# Load necessary packages
library(Seurat)
library(CellChat)
library(patchwork)

# set working directory
results_folder = '/R/nanostring7/work5cellchat/all6'

seu.obj <- readRDS("/R/nanostring7/6 Analysis/Data objects/seurat_object.Rds")

# change the name of "monocyte" in nb_clus_4 to "neutrophil"
seu.obj@meta.data$nb_clus_4[which(seu.obj@meta.data$nb_clus_4 == "monocyte")] <- "neutrophil"
seu.obj@meta.data$nb_clus_4[which(seu.obj@meta.data$nb_clus_4 == "e")] <- "B.cell.linage"

# Retrieve metadata
meta <- seu.obj@meta.data

# Retrieve data from each layer
data_matrices <- lapply(paste0("counts.", 1:7), function(layer) {
  GetAssayData(seu.obj[['Nanostring']], layer = layer)
})

# Check data matrices
lapply(data_matrices, dim)

# Integrate data
data_matrix <- do.call(cbind, data_matrices)

# Check consistency of metadata
all(rownames(meta) %in% colnames(data_matrix))
meta <- meta[match(colnames(data_matrix), rownames(meta)), ]

# CosMx
spatial.locs = GetTissueCoordinates(seu.obj, scale = NULL, cols = c("x", "y"))
spatial.locs <- spatial.locs[, c("x", "y")]
conversion.factor = 0.18

set.seed(123)  # Set seed for reproducibility
sampled_locs <- spatial.locs[sample(1:nrow(spatial.locs), size = 1000), ]
d <- computeCellDistance(sampled_locs)

spot.size = min(d) * conversion.factor # converting the distance in Pixels to Micrometers
spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size / 2)
meta$samples <- seu.obj$tissue  # Example: use `orig.ident` from the Seurat object

# Create CellChat object
cellchat <- createCellChat(object = data_matrix, meta = meta, group.by = "nb_clus_4",
                           datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors)

# Set intercellular communication database in the CellChat object
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB

# Preprocessing to identify intercellular communication networks
cellchat <- subsetData(cellchat)

# Calculate intercellular communication networks
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Compute communication network
cellchat <- computeCommunProb(cellchat, 
                              distance.use = TRUE, 
                              contact.range = 50,  # Communication within a range of 50 micrometers
                              scale.distance = 20,  # Increased scaling of distance
                              population.size = TRUE)

cellchat <- filterCommunication(cellchat, min.cells = 10)

# Calculate communication network pathways
cellchat <- computeCommunProbPathway(cellchat)

# Aggregate the network
cellchat <- aggregateNet(cellchat)

# Visualize communication network
groupSize <- as.numeric(table(cellchat@idents))

# Visualize the number of interactions
png(filename = file.path(results_folder, "Number_of_interactions.png"), 
    width = 7, height = 7, units = "in",  res = 600)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, vertex.label.cex = 1.0,
                 weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
dev.off()

# Visualize interaction strength
png(filename = file.path(results_folder, "Interaction_weights.png"), 
    width = 10, height = 10, units = "in", res = 300)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, 
                 vertex.label.cex = 1.0, label.edge = FALSE, title.name = "Interaction weights/strength")
dev.off()

png(filename = file.path(results_folder, "Interaction_weights_sq.png"), 
    width = 10, height = 10, units = "in", res = 300)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, 
                 vertex.label.cex = 1.0, label.edge = FALSE, shape = "square",
                 title.name = "Interaction weights/strength")
dev.off()

# Pathway-level analysis
pathways.show <- c("SPP1")  # Display the "SPP1" pathway as an example

# Circle plot
png(filename = file.path(results_folder, "Circle_plot_SPP1.png"), 
    width = 10, height = 10, units = "in", res = 300)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle",vertex.label.cex = 1.0)
dev.off()

# Check data for heatmap plot
if (length(unique(cellchat@netP[[pathways.show]])) > 1) {
  png(filename = file.path(results_folder, "Heatmap_SPP1.png"))
  netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
  dev.off()
} else {
  cat("Error: Not enough distinct break values for heatmap.\n")
}

# Communication between each cell
png(filename = file.path(results_folder, "Contribution_SPP1.png"), 
    width = 3, height = 3, units = "in", res = 600)
netAnalysis_contribution(cellchat, signaling = pathways.show,font.size = 10,font.size.title =10)
dev.off()

# Overall heatmap
png(filename = file.path(results_folder, "Overall_Heatmap.png"), 
    width = 10, height = 10, units = "in", res = 300)
netVisual_heatmap(cellchat)
dev.off()

# Network characteristic analysis
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Visualize hub cells
png(filename = file.path(results_folder, "SignalingRole_Network_SPP1.png"), 
    width = 10, height = 10, units = "in", res = 300)
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show,
                                  width = 11,
                                  height = 4,
                                  font.size = 10,font.size.title = 14)
dev.off()

# Role analysis
png(filename = file.path(results_folder, "SignalingRole_Scatter.png"), 
    width = 7, height = 7, units = "in", res = 300)
netAnalysis_signalingRole_scatter(cellchat,label.size = 6,font.size = 20)
dev.off()

# Hub cells bubble plot
png(filename = file.path(results_folder, "Bubble_Plot_SPP1.png"), 
    width = 10, height = 10, units = "in", res = 300)
netVisual_bubble(cellchat, sources.use = NULL, targets.use = NULL, signaling = pathways.show, remove.isolate = FALSE)
dev.off()

# Detailed pathway analysis
png(filename = file.path(results_folder, "Aggregate_Plot_SPP1.png"), 
    width = 7, height = 7, units = "in", res = 600)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle",vertex.label.cex = 1.0)
dev.off()

png(filename = file.path(results_folder, "Aggregate_Plot_chord_SPP1.png"),
    height = 900, width = 900)
netVisual_aggregate(cellchat, signaling = pathways.show, vertex.label.cex = 1.5, layout = "chord")
dev.off()

# Check data for heatmap plot
if (length(unique(cellchat@netP[[pathways.show]]$prob)) > 1) {
  png(filename = file.path(results_folder, "Heatmap_Detailed_SPP1.png"))
  netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
  dev.off()
} else {
  cat("Error: Not enough distinct break values for detailed heatmap.\n")
}

png(filename = file.path(results_folder, "Contribution_Detailed_SPP1.png"), 
    width = 10, height = 10, units = "in", res = 300)
netAnalysis_contribution(cellchat, signaling = pathways.show,font.size = 20,font.size.title = 20)
dev.off()

# Extract list of pathways from cellchat@netP$pathways and store in pathways_list
pathways_list <- cellchat@netP$pathways

# Check the extracted pathway names
print(pathways_list)

# Loop through each pathway
for (pathway in pathways_list) {
  
  pathways.show <- c(pathway)
  
  # Circle plot
  png(filename = file.path(results_folder, paste0("Circle_plot_", pathway, ".png")), 
      width = 10, height = 10, units = "in", res = 300)
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", vertex.label.cex = 1.0)
  dev.off()
  
  # Check data for heatmap plot
  if (length(unique(cellchat@netP[[pathways.show]])) > 1) {
    png(filename = file.path(results_folder

# Data validation for heatmap plot
if (length(unique(cellchat@netP[[pathways.show]])) > 1) {
    png(filename = file.path(results_folder, paste0("Heatmap_", pathway, ".png")))
    netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
    dev.off()
} else {
    cat(paste0("Error: Not enough distinct break values for heatmap for ", pathway, ".\n"))
}

# Intercellular communication
png(filename = file.path(results_folder, paste0("Contribution_", pathway, ".png")), 
    width = 3, height = 3, units = "in", res = 600)
netAnalysis_contribution(cellchat, signaling = pathways.show, font.size = 10, font.size.title = 10)
dev.off()

# Hub cell visualization
png(filename = file.path(results_folder, paste0("SignalingRole_Network_", pathway, ".png")), 
    width = 10, height = 10, units = "in", res = 300)
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show,
                                  width = 11, height = 4,
                                  font.size = 10, font.size.title = 14)
dev.off()

# Hub cell bubble plot
png(filename = file.path(results_folder, paste0("Bubble_Plot_", pathway, ".png")), 
    width = 10, height = 10, units = "in", res = 300)
netVisual_bubble(cellchat, sources.use = NULL, targets.use = NULL, signaling = pathways.show, remove.isolate = FALSE)
dev.off()

# Detailed pathway analysis
png(filename = file.path(results_folder, paste0("Aggregate_Plot_", pathway, ".png")), 
    width = 7, height = 7, units = "in", res = 600)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", vertex.label.cex = 1.0)
dev.off()

png(filename = file.path(results_folder, paste0("Aggregate_Plot_chord_", pathway, ".png")),
    height = 900, width = 900)
netVisual_aggregate(cellchat, signaling = pathways.show, vertex.label.cex = 1.5, layout = "chord")
dev.off()

# Data validation for detailed heatmap plot
if (length(unique(cellchat@netP[[pathways.show]]$prob)) > 1) {
    png(filename = file.path(results_folder, paste0("Heatmap_Detailed_", pathway, ".png")))
    netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
    dev.off()
} else {
    cat(paste0("Error: Not enough distinct break values for detailed heatmap for ", pathway, ".\n"))
}

# Detailed contribution analysis
png(filename = file.path(results_folder, paste0("Contribution_Detailed_", pathway, ".png")), 
    width = 10, height = 10, units = "in", res = 300)
netAnalysis_contribution(cellchat, signaling = pathways.show, font.size = 20, font.size.title = 20)
dev.off()

# All pathways
pathways_list <- cellchat@netP$pathways
print(pathways_list)

# Intercellular communication (all pathways)
png(filename = file.path(results_folder, "Contribution.png"), 
    width = 4, height = 3, units = "in", res = 600)
netAnalysis_contribution(cellchat, signaling=pathways_list, font.size = 10,font.size.title =10)
dev.off()

saveRDS(cellchat, file = file.path(results_folder, "cellchat_nanostring7.rds"))

library(NMF)
library(ggalluvial)

# Outgoing signaling
png(filename = file.path(results_folder, "OutgoingPattern.png"), 
    width = 7, height = 4, units = "in", res = 300)
selectK(cellchat, pattern = "outgoing")
dev.off()

# Outgoing pattern heatmap
png(filename = file.path(results_folder, "OutgoingPatternheatmap.png"), 
    width = 10, height = 6, units = "in", res = 300)
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
dev.off()

# Outgoing river plot
png(filename = file.path(results_folder, "outgoing_riverplot.png"), 
    width = 6, height = 6, units = "in", res = 600)
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()

# Outgoing dot plot
png(filename = file.path(results_folder, "outgoing_dotplot.png"), 
    width = 6, height = 6, units = "in", res = 600)
netAnalysis_dot(cellchat, dot.size = c(1, 6), font.size = 14, font.size.title = 14, 
                pattern = "outgoing")
dev.off()

# Incoming signaling
png(filename = file.path(results_folder, "IncomingPattern.png"), 
    width = 7, height = 4, units = "in", res = 300)
selectK(cellchat, pattern = "incoming")
dev.off()

# Incoming pattern heatmap
png(filename = file.path(results_folder, "IncomingPatternheatmap.png"), 
    width = 10, height = 6, units = "in", res = 300)
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
dev.off()

# Incoming river plot
png(filename = file.path(results_folder, "incoming_riverplot.png"), 
    width = 6, height = 3, units = "in", res = 600)
netAnalysis_river(cellchat, pattern = "incoming")
dev.off()

# Incoming dot plot
png(filename = file.path(results_folder, "incoming_dotplot.png"), 
    width = 5, height = 3, units = "in", res = 600)
netAnalysis_dot(cellchat, pattern = "incoming")
dev.off()

library(reticulate)
reticulate::py_install("umap-learn")

# Identify signaling groups based on functional similarity
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")

library(future)

# Run clustering
cellchat <- netClustering(cellchat, type = "functional")

umap <- reticulate::import("umap")
print(umap)

# Classification learning of signaling networks for a single dataset
# Visualization in 2D-space
png(filename = file.path(results_folder, "func_signaling_groups.png"), 
    width = 5, height = 3, units = "in", res = 600)
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
dev.off()

# Identify signaling groups based on structural similarity
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")

# Visualization in 2D-space
png(filename = file.path(results_folder, "struct_signaling_groups.png"), 
    width = 5, height = 3, units = "in", res = 600)
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
dev.off()

png(filename = file.path(results_folder, "struct_signaling_groups_dots.png"), 
    width = 5, height = 3, units = "in", res = 600)
netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)
dev.off()

sample_list <- unique(cellchat@meta$sample)

# SPP1-FN1 spatial plot for specific samples
png(filename = file.path(results_folder, "18_72_TB_SPP1-FN1-spatial.png"), 
    width = 50, height = 30, units = "in", res = 300)
spatialFeaturePlot(cellchat, features = c("SPP1","FN1"), sample.use = "18_72_TB", 
                   point.size = 1, color.heatmap = "Reds", direction = 1)
dev.off()

# CD44-SPP1 ligand-receptor pair plot
png(filename = file.path(results_folder, "18_72_TB_CD44_SPP1-spatial.png"), 
    width = 90, height = 40, units = "in", res = 300)
spatialFeaturePlot(cellchat, pairLR.use = "SPP1_CD44", sample.use = "18_72_TB", 
                   point.size = 0.1, do.binary = TRUE, cutoff = 0.05, 
                   enriched.only = F, color.heatmap = "Reds", direction = 1)
dev.off()

# Loop through each sample
for (sample_name in sample_list) {
  tryCatch({
    # Set the file name for each sample
    file_name <- paste0("CD44-SPP1_", sample_name, "-spatial.png")
    
    # Open the PNG device
    png(filename = file.path(results_folder, file_name), 
        width = 90, height = 40, units = "in", res = 300)
    
    # Execute the plot
    spatialFeaturePlot(cellchat, pairLR.use = "SPP1_CD44", sample.use = sample_name, 
                       point.size = 0.1, do.binary = TRUE, cutoff = 0.05, 
                       enriched.only = FALSE, color.heatmap = "Reds", direction = 1)
    
    # Close the PNG device
    dev.off()
    
    # Success message
    message(paste("Successfully created plot for", sample_name))
    
    # Free up memory
    gc()
    
  }, error = function(e) {
    # Display the error message and skip
    message(paste("Error in sample:", sample_name, "-", e$message))
  })
}




png(filename = file.path(results_folder, "18_72_TB_incoming-spatial.png"), 
    width = 90, height = 40, units = "in", res = 300)
netVisual_aggregate(cellchat, signaling = pathways.show, sample.use = "18_72_TB", 
                    layout = "spatial", edge.width.max = 2, alpha.image = 0.2, 
                    vertex.weight = "incoming", vertex.size.max = 6, 
                    vertex.label.cex = 0)
dev.off()



