# Libraries
library(tictoc)
library(igraph)
library(pheatmap)
library(scales)
library(Matrix)
library(viridis)
library(InSituType)
library(SeuratObject)
library(SeuratDisk)
library(InSituCor)

# Load Giotto object
load(".data/complete_giotto_object_PostAnalysis.RData")

#change the name of "monocyte" in nb_clus_4 to "neutrophil".
gem@cell_metadata$rna$nb_clus_4[which(gem@cell_metadata$rna$nb_clus_4 == "monocyte")] <- "neutrophil"
gem@cell_metadata$rna$nb_clus_4[which(gem@cell_metadata$rna$nb_clus_4 == "e")] <- "B.cell.linage"

annot <- as.data.frame(gem@cell_metadata$rna)
rownames(annot) <- annot$cell_ID
clust <- gem@cell_metadata$rna$nb_clus_4

## Normalized counts, cells x genes matrix
normcounts <- t(gem@expression$rna$normalized)
raw <- t(gem@expression$rna$raw)



# Calculate column means for the dgCMatrix
negmean <- colMeans(gem@expression$negprobes$raw)

# Convert the result to a data frame if needed
negmean_df <- as.data.frame(negmean)

# First, add the negative counts mean to metadata
gem@cell_metadata$rna$negmean <- colMeans(gem@expression$negprobes$raw)

# Next, extract a data frame of variables to condition on
conditionon <- as.data.frame(gem@cell_metadata$rna)[, c("fov", "totalcounts", "negmean", "nb_clus_4")]

# Add cell IDs as row names
rownames(conditionon) <- gem@cell_ID

## XY positions
xy <- gem@spatial_locs$viz[,1:2]

## Cell types
celltype <- gem@cell_metadata$rna$nb_clus_4

### Run InSituCor with default values
res <- insitucor(
  counts = normcounts,
  conditionon = conditionon,
  celltype = celltype,
  xy = xy,
  k = 50,
  tissue = gem@cell_metadata$rna$tissue
)

# set working directory
results_folder = './results'

# Check the class of 'raw'
class(raw)

# If it's a vector, convert it to a matrix (assuming a single row)
if (is.vector(raw)) {
  raw <- matrix(raw, nrow = 1)
}

# normalize:
norm = sweep(raw, 1, rowSums(raw), "/") * mean(rowSums(raw))

annot <- cbind(annot, negmean_df)

#### run whole SPARC pipeline at once to test runtime ---------------------------

tic()
res <- insitucor(counts = norm, 
                 conditionon = cbind(annot[, c("totalcounts", "negmean")], clust),
                 celltype = clust,
                 neighbors = NULL,
                 xy = xy,
                 k = 50)

toc()  #2.4945

#### run SPARC step by step ---------------------------------

step1 <- calcSpatialCor(counts = norm, 
                        conditionon = cbind(annot[, c("totalcounts", "negmean")], clust),
                        neighbors = NULL,
                        xy = xy,
                        k = 50,
                        radius = NULL,
                        tissue = NULL, # example of using tissue in neighbor definition
                        roundcortozero = 0.1,
                        max_cells = 10000,
                        verbose = FALSE)

saveRDS(step1, file = "/R/nanostring7/insitucor/processed_data/step1.RDS")


# derive Modules
modules <- defineModules(condcor = step1$condcor,
                         env = step1$env,
                         min_module_size = 2,
                         max_module_size = 20,
                         resolution = 0.02,
                         corthresh = 0.1, 
                         min_module_cor = 0.1,
                         gene_weighting_rule = "inverse_sqrt")
saveRDS(modules, file = "/R/nanostring7/insitucor/processed_data/modules.RDS")


# get scores:
scores <- scoreModules(counts = norm,
                       weights = modules$weights,
                       neighbors = step1$neighbors)
saveRDS(scores, file = "/R/nanostring7/insitucor/processed_data/scores.RDS")


#### attribution analysis -----------------------------------------------

tic()
attribution <- cellTypeAttribution(
  modulescores = scores$scores_env,
  weights = modules$weights,
  counts = norm,
  celltype = clust,
  neighbors = step1$neighbors,
  nsub = 1000,
  verbose = TRUE)
toc()
saveRDS(attribution, file = "/R/nanostring7/insitucor/processed_data/attribution.RDS")


#### plots for manuscript --------------------------------------------

# cell type map:
cols = InSituType::colorCellTypes(freqs = table(clust), palette = "brewers")

png("/R/nanostring7/insitucor/results/fig 1a - cell type map.png", width = 6, height = 4, units= "in", res = 2400)
par(mar = c(0,0,0,0))
plot(xy, pch = 16, col = cols[clust], cex = 0.1, asp = 1, ylim = c(-7,0),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(c(2.1,3.1), c(-0.2,-0.2))
text(2.6,-0.3,"1 mm", cex = 0.8)
dev.off()

svg("/R/nanostring7/insitucor/results/fig 1a - cell type legend.svg", width = 3, height = 3)
par(mar = c(0,0,0,0))
frame()
o = order(names(cols), decreasing = F)
legend("center", pch = 16, col = cols[o], legend = gsub(" naive", "", gsub(" memory", "", names(cols)[o])), cex = 0.55, ncol = 2)
dev.off()

# cartoon of neighbors:
ind = 10000

# Extract columns from data.table as vectors
x_coord = xy[[1]]  # First column
y_coord = xy[[2]]  # Second column
x_ind = xy[[1]][ind]  # First column of the ind-th row
y_ind = xy[[2]][ind]  # Second column of the ind-th row

# Calculate 'use'
use = ((x_coord > x_ind - 0.1) & (x_coord < x_ind + 0.1)) & 
  ((y_coord > y_ind - 0.1) & (y_coord < y_ind + 0.1))
ind.neighbors = which(step1$neighbors[ind, ] > 1e-6)

# Create SVG for nearest neighbors plot
svg("/R/nanostring7/insitucor/results/fig 1b - nearest neighbors.svg", width = 5, height = 5)
par(mar = c(0, 0, 0, 0))

# Plot points based on the 'use' condition
plot(xy[use, ], pch = 16, col = cols[clust[use]], cex = 0.8, asp = 1,
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
points(xy[ind, 1], xy[ind, 2], cex = 1, pch = 16)

# Plot nearby points
# Extract x and y coordinates as vectors
x_coords_neighbors = xy[ind.neighbors, 1, with = FALSE][[1]]
y_coords_neighbors = xy[ind.neighbors, 2, with = FALSE][[1]]

points(x_coords_neighbors, y_coords_neighbors, cex = 1, pch = 1, col = "blue")

# Close the SVG file
dev.off()


# heatmap of environment matrix:
set.seed(0)
subr = sample(1:nrow(step1$env), 500)
subc = order(colMeans(step1$env), decreasing = T)[1:100]

mat = step1$env[subr, subc]
mat = sweep(mat,2,apply(mat,2,max), "/")
p1 = pheatmap(mat)
dev.off()
png("/R/nanostring7/insitucor/results/fig 1c - env matrix.png", width = 3, height = 3, units ="in",res = 400)
pheatmap(mat[p1$tree_row$order, p1$tree_col$order], 
         cluster_rows = F, cluster_cols = F,
         show_rownames = F, show_colnames = F, legend = FALSE,
         col = viridis_pal(option = "B")(100))
dev.off()

# heatmap of cormat
rawcor = cor(step1$env)
saveRDS(rawcor, file = "/R/nanostring7/insitucor/processed_data/rawcor.RDS")


set.seed(0)
inds = sample(1:nrow(rawcor), 300)

hc1 = hclust(dist(rawcor[inds, inds]))
png("/R/nanostring7/insitucor/results/fig 1e - raw cormat.png", width = 3, height = 3, units ="in",res = 400)
pheatmap(rawcor[inds, inds][hc1$order, hc1$order], cluster_rows = F, cluster_cols = F,
         col = colorRampPalette(c("darkblue",'blue', "white","red","darkred"))(100),
         breaks = seq(-0.6,0.6,length.out = 101),
         #col = colorRampPalette(c('blue', "white","red"))(100),
         #breaks = seq(-0.4,0.4,length.out = 101),
         show_rownames = F, show_colnames = F, legend = FALSE)
dev.off()

# heatmap of cond cor:
hc2 = hclust(dist(step1$condcor[inds, inds]))
png("/R/nanostring7/insitucor/results/fig 1f - conditional cormat.png", width = 3, height = 3, units ="in",res = 400)
pheatmap(step1$condcor[inds, inds][hc2$order, hc2$order], cluster_rows = F, cluster_cols = F,
         col = colorRampPalette(c("darkblue",'blue', "white","red","darkred"))(100),
         breaks = seq(-0.6,0.6,length.out = 101),
         show_rownames = F, show_colnames = F, legend = FALSE)
dev.off()

# heatmap of env confounders
condmat = InSituCor:::build_conditional_matrix(cbind(annot[, c("totalcounts", "negmean")], clust))
colnames(condmat)[1:2] = c("total counts", "negprobe counts")
condmat[, "total counts"] = condmat[, "total counts"] / mean(condmat[, "total counts"])

mat = condmat[1:500, order(colMeans(condmat[1:500, ]), decreasing = T)[1:20]]
p1 = pheatmap(mat)
dev.off()
png("/R/nanostring7/insitucor/results/fig 1d - conditioning matrix.png", width = 3, height = 3, units ="in",res = 400)
pheatmap(mat[p1$tree_row$order, p1$tree_col$order], cluster_rows = F, cluster_cols = F,
         #fontsize_col = 6,
         show_colnames = F,
         col = colorRampPalette(c("grey80","darkviolet"))(100),
         show_rownames = F, legend = FALSE)
dev.off()


## show examples of marker genes losing conditional cor:

# find marker genes:
profiles = InSituType::Estep(counts = raw, clust = clust, neg = annot$neg)
markers = list()
for (cell in unique(clust)) {
  markers[[cell]] = rownames(profiles)[
    order(profiles[, cell] / apply(profiles[, setdiff(colnames(profiles), cell)], 1, max) * sqrt(profiles[, cell]),
          decreasing = T)[1:4]]
}


# get non-rounded conditional cor of only marker genes:
markerstep1 <- calcSpatialCor(counts = norm[, unlist(markers)], 
                              conditionon = cbind(annot[, c("totalcounts", "negmean")], clust),
                              neighbors = NULL,
                              xy = xy,
                              k = 50,
                              radius = NULL,
                              tissue = NULL, # example of using tissue in neighbor definition
                              roundcortozero = 1e-4,
                              max_cells = 10000,
                              verbose = FALSE)

unroundedstep1 <- calcSpatialCor(counts = norm, 
                                 conditionon = cbind(annot[, c("totalcounts", "negmean")], clust),
                                 neighbors = NULL,
                                 xy = xy,
                                 k = 50,
                                 radius = NULL,
                                 tissue = NULL, # example of using tissue in neighbor definition
                                 roundcortozero = 1e-2,
                                 max_cells = 10000,
                                 verbose = FALSE)
saveRDS(unroundedstep1, file = "/R/nanostring7/insitucor/processed_data/unroundedstep1.RDS")



# look at cors:
gpairs = list(c("MS4A1", "CD19"),   #, "BLK", "BCL2"
              c("CD68", "CD163"),
              c("CD3E", "CD3D"), 
              c("COL5A1", "COL5A2"))

for (i in 1:length(gpairs)) {
  print(c(rawcor[gpairs[[i]][1], gpairs[[i]][2]],
          unroundedstep1$condcor[gpairs[[i]][1], gpairs[[i]][2]]))
}
for (i in 1:length(gpairs)) {
  print(rawcor[gpairs[[i]], gpairs[[i]]])
}

for (i in 1:length(gpairs)) {
  print(step1$condcor[gpairs[[i]], gpairs[[i]]])
}
#次はなぜかできない。
for (i in 1:length(gpairs)) {
  print(markerstep1$condcor[gpairs[[i]], gpairs[[i]]])
}


## summary of cors before and after:

# summarize how many cors are lost after conditioning:
top5kthresh = quantile(rawcor[upper.tri(rawcor)], 1 - 5000/sum(upper.tri(rawcor)))
top5k = (rawcor > top5kthresh) & upper.tri(rawcor)
sum(step1$condcor[top5k] > 0.2)
sum(step1$condcor[top5k] <= 0.2)
sum(top5k)


# plot cors before and after:
ynudge = c(0.035,-0.05,0.05,0.035)
# Confirm that there are no NAs or NaNs, then extract the upper triangular portion
valid_rawcor <- rawcor[upper.tri(rawcor)]
valid_condcor <- unroundedstep1$condcor[upper.tri(step1$condcor)]

# Extract only the matching data points
valid_index <- !is.na(valid_rawcor) & !is.na(valid_condcor)
valid_rawcor <- valid_rawcor[valid_index]
valid_condcor <- valid_condcor[valid_index]

# Perform the plot with matching data points in the matrix
png("/R/nanostring7/insitucor/results/fig 1g - scatterplot of cors before and after.png", 
    width = 4, height = 5, units = "in", res = 300)
par(mar = c(5, 4.5, 0, 0))
plot(valid_rawcor, valid_condcor, 
     pch = 16, cex = 0.5,
     xlab = "Correlation of gene pairs\nin environment matrix", 
     ylab = "Correlation conditional on confounders",
     cex.lab = 1.1,
     xlim = c(0, 1), ylim = c(-0.1, 1),
     col = alpha("dodgerblue4", 0.25))

# Close the PNG file
dev.off()


abline(0,1)
abline(h = 0, col = 1)
abline(v = 0, col = 1)
for (i in 1:length(gpairs)) {
  points(rawcor[gpairs[[i]][1], gpairs[[i]][2]],
         step1$condcor[gpairs[[i]][1], gpairs[[i]][2]],
         pch = 16, col = "darkred")
  text(rawcor[gpairs[[i]][1], gpairs[[i]][2]],
       step1$condcor[gpairs[[i]][1], gpairs[[i]][2]] + ynudge[i],
       labels = paste0(gpairs[[i]][1], "/", gpairs[[i]][2]),
       col = "darkred", cex = 0.75)
}
rect(top5kthresh, 0.205, 1, 1, border = "orange")
rect(top5kthresh, -0.1, 1, 0.195, border = "darkviolet")
text(0.85,-0.05,paste0(sum(step1$condcor[top5k] <= 0.2), " pairs"), col = "darkviolet")
text(0.85,0.85,paste0(sum(step1$condcor[top5k] > 0.2), " pairs"), col = "orange")
dev.off()



### spatial plots -----------------------
name = "FN1_SPP1_MARCO_3"
png(paste0("/R/nanostring7/insitucor/results/fig 1h - ", name, "env scores.png"), width = 6, height = 4, units= "in", res = 2400)
par(mar = c(0,0,0,0))
plot(xy, pch =16, asp = 1, cex = 0.2, ylim = c(-7,0),
     col = viridis_pal(option = "B")(101)[
       1 + round(100 * pmin(scores$scores_env[, name] / quantile(scores$scores_env[, name], 0.95), 1))],
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(c(2.1,3.1), c(-0.2,-0.2))
text(2.6,-0.3,"1 mm", cex = 0.8)
#rect(0.2, 2.5, 1, 3.5, border = "cyan", lwd = 2)
rect(1.4, 2.5, 3, 4.5, border = "dodgerblue1", lwd = 2)
dev.off()

png(paste0("/R/nanostring7/insitucor/results/fig 1i - ", name, "env scores zoom1.png"), width = 8, height = 5, units= "in", res = 300)
par(mar = c(0,0,0,0))
plot(xy, pch =16, asp = 1, 
     xlim = c(0.2,5), ylim = c(-3.5,-0.2),
     cex = 0.3,
     col = viridis_pal(option = "B")(101)[
       1 + round(100 * pmin(scores$scores_env[, name] / quantile(scores$scores_env[, name], 0.95), 1))],     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
dev.off()

png(paste0("/R/nanostring7/insitucor/results/fig 1j - ", name, "sc scores zoom.png"), 
    width = 7, height = 5, units= "in", res = 300)
par(mar = c(0,0,0,0))
plot(xy, pch =16, asp = 1, 
     ylim = c(-7,0),
     cex = .2 + 0.2 * pmin(scores$scores_sc[, name] / quantile(scores$scores_sc[, name],0.99),1), 
     col = cols[clust],
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
showcells = c("neutrophil", "macrophage_d")
legend("bottomleft", pch = 16, col = cols[showcells], legend = showcells, cex = 0.8)
rect(1.55, 3.55, 1.95, 4.05, border = "dodgerblue1", lwd = 2)
dev.off()

### network diagram ------------------

svg("/R/nanostring7/insitucor/results/fig 1k - network.svg")
par(mar = c(0,0,0,0))
plotCorrelationNetwork(x = step1$condcor,
                       modules = modules$weightsdf, show_gene_names = F)
dev.off()

svg("/R/nanostring7/insitucor/results/fig 1k - network v2.svg", width = 8, height = 8)
par(mar = c(0,0,0,0))
plotCorrelationNetwork(x = step1$condcor,
                       modules = modules$weightsdf, show_gene_names = F, vertex_size = 3.5)
dev.off()

pdf("/R/nanostring7/insitucor/results/fig 1k - network w names.pdf")
par(mar = c(0,0,0,0))
plotCorrelationNetwork(x = step1$condcor,
                       modules = modules$weightsdf, corthresh = 0.2, show_gene_names = T,  
                       vertex_size = 0.1)
dev.off()


### attribution plots --------------

svg("/R/nanostring7/insitucor/results/attribution matrix transposed.svg", height = 8, width = 14)
pheatmap(t(attribution$involvescores),
         border_color = NA, legend = F,
         col = colorRampPalette(c("white","darkblue"))(100),
         breaks = seq(0,0.9,length.out = 101))
#fontsize_col = 6)
dev.off()

p1 = pheatmap(attribution$attributionmats[[name]], main = name,
              col = colorRampPalette(c("white","darkblue"))(100),
              breaks = seq(0,1,length.out=101))

svg(paste0("/R/nanostring7/insitucor/results/attribution for ", name, ".svg"), width = 6, height = 7)
pheatmap(attribution$attributionmats[[name]][p1$tree_row$order, p1$tree_col$order],
         border_color = NA, 
         cluster_rows = F, cluster_cols = F,
         col = colorRampPalette(c("white","darkblue"))(100),
         breaks = seq(0,0.9,length.out = 101),
         fontsize_row = 10, fontsize_col = 8, legend = F)
dev.off()


pdf("/R/nanostring7/insitucor/results/attributionmats.pdf")
for (name in names(attribution$attributionmats)) {
  pheatmap(attribution$attributionmats[[name]],
           main = name,
           col = colorRampPalette(c("white","darkblue"))(100),
           breaks = seq(0,1,length.out=101))
}
dev.off()


load("/R/nanostring7/insitucor/CellChatDB.human.rda") # (from https://github.com/sqjin/CellChat/blob/master/data/CellChatDB.human.rda)
ligands = intersect(unique(CellChatDB.human$interaction$ligand), colnames(raw))
receptors = intersect(unique(CellChatDB.human$interaction$receptor), colnames(raw))
lrpairs = CellChatDB.human$interaction[, c("ligand", "receptor")]
rownames(lrpairs) = CellChatDB.human$interaction$interaction_name
lrpairs = lrpairs[is.element(lrpairs[,1], colnames(raw)) & is.element(lrpairs[,2], colnames(raw)), ]
lrpairs = lrpairs[lrpairs[,1] != lrpairs[,2], ]


step1 = readRDS(file = "/R/nanostring7/insitucor/processed_data/step1.RDS")
modules <- readRDS("/R/nanostring7/insitucor/processed_data/modules.RDS")
scores <- readRDS("/R/nanostring7/insitucor/processed_data/scores.RDS")
attribution = readRDS("/R/nanostring7/insitucor/processed_data/attribution.RDS")
rawcor = readRDS("/R/nanostring7/insitucor/processed_data/rawcor.RDS")
unroundedstep1 = readRDS(file = "/R/nanostring7/insitucor/processed_data/unroundedstep1.RDS")


# get condcor
ligstep1 <- calcSpatialCor(counts = norm[, ligands], 
                           conditionon = cbind(annot[, c("log10totalcounts", "background")], clust),
                           neighbors = NULL,
                           xy = xy,
                           k = 50,
                           radius = NULL,
                           tissue = NULL, # example of using tissue in neighbor definition
                           roundcortozero = 1e-4,
                           max_cells = 10000,
                           verbose = FALSE)
saveRDS(ligstep1, file = "/R/nanostring/insitucor/processed_data/ligstep1.RDS")

# get modules:
lmods <- defineModules(condcor = ligstep1$condcor, 
                       env = ligstep1$env, 
                       min_module_size = 2,
                       max_module_size = 20,
                       resolution = 0.02,
                       corthresh = 0.1, 
                       min_module_cor = 0.1,
                       gene_weighting_rule = "inverse_sqrt")
lmods$modules

lscores =  scoreModules(counts = norm,
                        weights = lmods$weights,
                        neighbors = ligstep1$neighbors)

pheatmap(cor(lscores$scores_env))

# attribution:
ligattribution <- cellTypeAttribution(
  modulescores = lscores$scores_env,
  weights = lmods$weights,
  counts = norm,
  celltype = clust,
  neighbors = ligstep1$neighbors,
  nsub = 1000,
  verbose = TRUE)
saveRDS(ligattribution, file = "/R/nanostring7/insitucor/processed_data/ligattribution.RDS")
ligattribution = readRDS("/R/nanostring7/insitucor/processed_data/ligattribution.RDS")


p1 = pheatmap(ligattribution$involvescores,
              col = colorRampPalette(c("white", "darkblue"))(100),
              breaks = seq(0,0.8,length.out=101))
svg("/R/nanostring7/insitucor/results/fig 2b - attrib heatmap.svg", width = 8, height = 6)
pheatmap(ligattribution$involvescores[p1$tree_row$order, p1$tree_col$order],
         cluster_rows = F, cluster_cols = F,
         col = colorRampPalette(c("white", "darkblue"))(100),
         breaks = seq(0,0.8,length.out=101), legend = FALSE)
dev.off()

svg("/R/nanostring7/insitucor/results/fig 2a - network.svg")
par(mar = c(0,0,0,0))
plotCorrelationNetwork(x = ligstep1$condcor,corthresh = 0.1, 
                       modules = lmods$weightsdf, show_gene_names = T)
dev.off()

# spatial plots of ligand env scores:
for (name in colnames(lscores$scores_env)) {
  
  png(paste0("/R/nanostring7/insitucor/results/ligand modules env scores - ", name, ".png"), 
      width = 12, height = 8, units= "in", res = 500)
  par(mar = c(0,0,2,0))
  plot(xy, pch =8, asp = 1, cex = 0.05, ylim = c(-7,0), 
       main = paste0(lmods$modules[[name]], collapse = ", "), cex.main = 0.65,
       col = viridis_pal(option = "B")(101)[
         1 + round(100 * pmin(lscores$scores_env[, name] / quantile(lscores$scores_env[, name], 0.95), 1))],
       xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  lines(c(2.1,3.1), c(-0.2,-0.2))
  text(2.6,-0.3,"1 mm", cex = 0.8)
  dev.off()
}



#### LR pair analysis -----------------------------------
# get condcor for just the lr genes:

# look at cor values for all LR pairs:
lrcors = lrrawcors = c()
for (name in rownames(lrpairs)) {
  lrcors[name] = step1$condcor[lrpairs[name, 1], lrpairs[name, 2]]
  lrrawcors[name] = rawcor[lrpairs[name, 1], lrpairs[name, 2]]
}

plot(lrrawcors, lrcors, asp = 1, col = 0)
abline(0,1)
abline(h = 0)
text(lrrawcors, jitter(lrcors, amount = 0.106), names(lrcors), cex = 0.5)

# histogram of cor values for all LR pairs, with lines at cor values for selected pairs
svg("/R/nanostring7/insitucor/results/fig2f - LR cor histogram.svg", height = 3, width = 4)
par(mar = c(4,4,0,1))
hist(lrcors, breaks = 50, main = "", xlab = "Conditional correlation",
     ylab = "Number of L-R pairs", col = alpha("dodgerblue4", 0.75))
dev.off()

lrcors[order(lrcors, decreasing = T)[1:20]]

# Convert lrcors to a data frame
lrcors_df <- data.frame(LR_Pair = names(lrcors), Correlation = lrcors)

# Save as a CSV file
write.csv(lrcors_df, file = "/R/nanostring7/insitucor/results/LRlist_all.csv", row.names = FALSE)


# pick 1-2 pairs, and show transcript plots:

# make pseudomodules, then score:
lrmods = list()
for (i in 1:5) {
  name = names(lrcors)[order(lrcors, decreasing = T)[i]]
  lrmods[[name]] = c(lrpairs[name, 1], lrpairs[name, 2])
}

#### correlation network around a LR pair: 
name = "SPP1_CD44"
gl = lrpairs[name, "ligand"]
gr = lrpairs[name, "receptor"]
usegenes = names(which(colSums(step1$condcor[c(gl, gr), ] > 0) > 0))

gr <- igraph::graph_from_data_frame(lrpairs)

gr0 <- InSituCor:::get_coregulation_network(cormat = step1$condcor[usegenes, usegenes], 
                                             corthresh = 0.1)
svg("/R/nanostring7/insitucor/results/fig 2h - network around LR pair.svg", width = 7, height = 7)
par(mar = c(0,0,0,0))
igraph::plot.igraph(gr0, 
                    vertex.label = unlist(list(NA, names(igraph::V(gr0)))[1 + TRUE]), 
                    vertex.size = 1, 
                    vertex.color = "grey",
                    vertex.label.color = c("dodgerblue4", "red")[1 + is.element(usegenes, c(gl, gr))]) 
dev.off()



# make a module score:
gcscore = rowMeans(norm[, setdiff(usegenes, c("SPP1", "CD44"))])
png(paste0("/R/nanostring7/insitucor/results/ligand modules env scores - ", name, ".png"),
    width = 6, height = 4, units= "in", res = 2400)
par(mar = c(0,0,0,0))
plot(xy, pch =16, asp = 1, cex = 0.1, ylim = c(-7,0),
     col = viridis_pal(option = "B")(101)[
       1 + round(100 * pmin(gcscore / quantile(gcscore, 0.995), 1))],
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
dev.off()



png("/R/nanostring7/insitucor/results/fig2g12 - spatial plot of 2 genes.png", 
    width = 6, height = 6, units= "in", res = 1200)
par(mar = c(0,0,0,0))
par(bg = "black")
show = rowSums(raw[, c(gl, gr)]) > 0
o = order(rowSums(norm[, c(gl, gr)]))
plot(xy[o, ], 
     col = c("grey20", "yellow", "cyan", "red")[
       1 + (norm[, gl] > 2) + 2 * (norm[, gr] > 2)
     ][o],
     pch = 16, cex = 0.1 + 0.05 * (apply(norm[, c(gl, gr)], 1, max) > 2), asp = 1,
     ylim = c(-7, 0), #ylim = c(-0.51,4.61),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend(7,1.5, pch = c(16,16,16),  cex = 0.7,  
       col = c("yellow", "cyan", "red"),
       legend = c(paste0("high ", gl), paste0("high ", gr), "high both"), 
       text.col = "white")
lines(c(3, 6), rep(2.7,2), col = "white")
text(2,1,"1 mm", cex = 0.8, col = "white")
dev.off()


# macrophage_d GC polygons:
gcpolys = list()
for (id in setdiff(names(which(table(db) > 100)), 0)) {
  gcpolys[[id]] = xy[clust == "macrophage_d", ][db == id, ][chull(xy[clust == "macrophage_d", ][db == id, ]), ]
}
name = "SPP1_CD44"
gl = lrpairs[name, "ligand"]
gr = lrpairs[name, "receptor"]

png("/R/nanostring7/insitucor/results/fig2g - spatial plot of 2 genes.png", 
    width = 6, height = 6, units= "in", res = 1200)
par(mar = c(0,0,0,0))
par(bg = "black")
show = rowSums(raw[, c(gl, gr)]) > 0
o = order(rowSums(norm[, c(gl, gr)]))
plot(xy[o, ], 
     col = c("grey20", "yellow", "cyan", "red")[
       1 + (norm[, gl] > 2) + 2 * (norm[, gr] > 2)
     ][o],
     pch = 16, cex = 0.1 + 0.05 * (apply(norm[, c(gl, gr)], 1, max) > 2), asp = 1,
     ylim = c(-7, 0), #ylim = c(-0.51,4.61),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
for (i in 1:length(gcpolys)) {
  polygon(gcpolys[[i]], col = NULL, border = "white")
}
legend(7,1.5, pch = c(16,16,16,NA), lty = c(NA,NA,NA,1), cex = 0.7,  
       #col = c(rgb(1,0,0,1), rgb(0,0,1,1), rgb(1,0,1,1), "white"),
       col = c("yellow", "cyan", "red", "white"),
       legend = c(paste0("high ", gl), paste0("high ", gr), "high both", "Lymphoid structure"), 
       text.col = "white")
lines(c(3, 6), rep(2.7,2), col = "white")
text(2,1,"1 mm", cex = 0.8, col = "white")
dev.off()

png("/R/nanostring7/insitucor/results/fig2g2 - spatial plot of 2 genes.png", 
    width = 6, height = 6, units= "in", res = 1200)
par(mar = c(0,0,0,0))
par(bg = "black")
show = rowSums(raw[, c(gl, gr)]) > 0
o = order(rowSums(norm[, c(gl, gr)]))
plot(xy[o, ], 
     col = c("grey20", "yellow", "cyan", "red")[
       1 + (norm[, gl] > 2) + 2 * (norm[, gr] > 2)
     ][o],
     pch = 16, cex = 0.1 + 0.05 * (apply(norm[, c(gl, gr)], 1, max) > 2), asp = 1,
     ylim = c(-7, 0), #ylim = c(-0.51,4.61),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend(7,1.5, pch = c(16,16,16),  cex = 0.7,  
       col = c("yellow", "cyan", "red"),
       legend = c(paste0("high ", gl), paste0("high ", gr), "high both"), 
       text.col = "white")
lines(c(3, 6), rep(2.7,2), col = "white")
text(2,1,"1 mm", cex = 0.8, col = "white")
dev.off()


