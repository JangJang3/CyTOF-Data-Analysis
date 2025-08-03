################################################################
# Set working environment 
################################################################
current_path <- "/mnt/volumeb/jin/baseline_pipeline/analysis/analysis_only_Tcells_test_run"

getwd()
setwd(dir=current_path)


################################################################
# Install and Load Packages 
################################################################

source("/mnt/volumeb/jin")


################################################################
# Create File Path and Folders
################################################################
print("Status:  Create File Paths")


dir_name <- "spread_clusters"

# create data folder where all analysis will be stored
if(!dir.exists(dir_name)) dir.create(dir_name)

# set it as a main directory 
dir <- file.path(getwd(), dir_name)


################################################################
# Add  annotations 
################################################################
print("Status:  Add Annotations")

srat_obj <- readRDS("seurat_object_clustering_leiden.rds")



# Annotations 
# Cluster 1 - CD127+CD28+CD27high TCM CD4 cells
# Cluster 2 -CD27+CD28midCD127low naïve CD4 cells
# Cluster 3 - CD127+CD28+CD27highCD161low TCM CD4 cells 
# Cluster 4 - CD314+Ki67highCD57midHLA-DRmidCD38lowCD56lowCD27low TEMRA CD8 cells
# Cluster 5 - CD314+CD27highCD28midCD127midPDL1low naive CD8 cells
# Cluster 6 - Tregs
# Cluster 7 - PD1midCD57highCD127low TEM CD4 cells
# Cluster 8 - VD2 cells
# Cluster 9 - VD1 cells
# Cluster 10 - CD127+CD161+CD56midCD28mid MAIT cells
# Cluster 11 - CD27+CD28+CD127highTCRVa7.2CD25lowICOSlow TCM CD4 cells 
# Cluster 12 - NKT2 cells 
# Cluster 13 - CD16+PD1+CD57+CD314highCD56low MAIT cells


Idents(srat_obj) <- "res_0.2"

cluster_ids.new <- c(
  "CD127+CD28+CD27high TCM CD4 cells",
  "CD27+CD28midCD127low naïve CD4 cells", 
  "CD127+CD28+CD27highCD161low TCM CD4 cells", 
  "CD314+Ki67highCD57midHLA-DRmidCD38lowCD56lowCD27low TEMRA CD8 cells", 
  "CD314+CD27highCD28midCD127midPDL1low naive CD8 cells",
  "Tregs", 
  "PD1midCD57highCD127low TEM CD4 cells",
  "VD2 cells",
  "VD1 cells",
  "CD127+CD161+CD56midCD28mid MAIT cells",
  "CD27+CD28+CD127highTCRVa7.2CD25lowICOSlow TCM CD4 cells",
  "NKT2 cells",
  "CD16+PD1+CD57+CD314highCD56low MAIT cells"
)
names(cluster_ids.new) <- levels(srat_obj)
srat_obj <- RenameIdents(srat_obj, cluster_ids.new)
saveRDS(srat_obj, file.path(dir, "seurat_object_annotated.rds"))

#number_of_clusters = n_distinct(srat_obj@meta.data["res_0.2"])


################################################################
# Select Colors
################################################################
print("Status:  Select Colors")


# met_palettes <- c("Archambault", "Derain", "Homer1", "Navajo", "Pillement", "Tiepolo", "Austria",
#                   "Egypt", "Homer2", "Klimt", "NewKingdom", "Pissaro", "Troy", "Gauguin", "Lakota",
#                   "Nizami", "Redon", "Greek", "Renoir", "VanGogh2", "Hiroshige", "Isfahan2",
#                   "Signac", "Java", "Thomas")
# colors <- MetBrewer::met.brewer("Signac", number_of_clusters)
# 
# 
# # Your hex color palette
# palette_colors <- c(
#   '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99', '#ff7899', '#337690',
#   '#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99'
# )
# 
# 
# # install.packages("paletteer")
# library(paletteer)
# 
# # Use in a ggplot2 chart:
# scale_colour_paletteer_d("MoMAColors::Klein")
# scale_fill_paletteer_d("MoMAColors::Klein")

colors <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  '#882255', "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)


################################################################
# Spread and Distance
################################################################
print("Status:  Test different min_dist and spreads")


umap_params <- expand.grid(
  min_dist = c(0.01, 0.1, 0.3, 0.5, 1.0),
  spread = c(0.5, 1.0)
)


#' Calculates the numbers of PCS needed to explain 95% of the variance in the data 
#' 
#' @param seurat_obj Seurat object 
#' @return The number of PCs needed to explain 95% variance in the data as an Integer
calculate_npcs <- function(seurat_obj){
  srat_obj_variance <- seurat_obj@reductions$pca@stdev**2
  srat_obj_variance <- cumsum(srat_obj_variance)/sum(srat_obj_variance)
  srat_obj_variance
  
  which(srat_obj_variance > 0.95)[1] #
  
  needed_pcs <- which(srat_obj_variance > 0.95)[1]
  
  print(str_interp("Needed PCs are: ${needed_pcs}"))
  
  return(needed_pcs)
}


umaps <- lapply(seq(nrow(umap_params)), function(i) {
  message(i)
  set.seed(2025)
  srat_obj <- RunUMAP(srat_obj,
                      dims = 1:calculate_npcs(srat_obj),
                      min.dist = umap_params$min_dist[i],
                      spread = umap_params$spread[i])
  
  name <- glue("spread: {umap_params$spread[i]}\nmin_dist: {umap_params$min_dist[i]}")
  p <- DimPlot(srat_obj, reduction = "umap", label = F) + 
    scale_color_manual(values = colors) + 
    theme(legend.position="none") + 
    ggtitle(name)
  
  return(p)
  #plot_list[[name]] <- p
})

# Display all plots in a grid
plots <- grid.arrange(grobs = umaps, ncol = 3)

ggsave(file.path(dir, "UMAPs_spreading.png"),plot = plots,  width = 15, height = 15)
