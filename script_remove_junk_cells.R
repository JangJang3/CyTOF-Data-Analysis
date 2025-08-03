################################################################
# Set working environment 
################################################################
current_path <- "/mnt/volumeb/jin/batch_no_HD2_test"

getwd()
setwd(dir=current_path)

################################################################
# Installing and Loading Packages 
################################################################


source("/mnt/volumeb/jin/installations.R")



################################################################
# Create File Path and Folders
################################################################
print("Status:  Create File Paths")

dir_name <- "remove_junks"

# create data folder where all analysis will be stored
if(!dir.exists(dir_name)) dir.create(dir_name)

# set it as a main directory 
dir <- file.path(getwd(), dir_name)


################################################################
# Remove Junk Cells
################################################################
print("Status:  Read Seurat Objects")

srat_obj <- readRDS(file.path(dir, "seurat_object_clustering_leiden.rds"))

res_values <- srat_obj@meta.data$res_0.1 

DimPlot(srat_obj, reduction = "umap", raster = FALSE, label = TRUE, group.by = "res_0.1") + ggtitle(str_interp("UMAP with Leiden Algorithm and res_0.1"))
ggsave(file.path(dir, "plot_umap_leiden_${res}.png"), width = 4000, height = 4000, units = "px")


ggplot(data.frame(res_values), aes(x = res_values)) +
  geom_bar(stat = "count", fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of res_0.1",
       x = "Clusters",
       y = "Frequency")

ggsave(file.path(dir, "plot_hist_res_0.1.png"), width = 4000, height = 4000, units = "px")





print("Status:  Remove Junk cells")

# 
junk_cells <- subset(srat_obj, subset = res_0.1 == 9)

useable_cells <- subset(srat_obj, subset = res_0.1 != 9)

useable_cells.sample_id <- useable_cells@meta.data$sample_id


ggplot(data.frame(useable_cells.sample_id), aes(x = useable_cells.sample_id)) +
  geom_bar(stat = "count", fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of Samples",
       x = "Clusters",
       y = "Frequency")

ggsave(file.path(dir, "plot_hist_useable_cells.png"), width = 4000, height = 4000, units = "px")



################################################################
# Save useable cells
################################################################
print("Status:  Save useable cells")

saveRDS(useable_cells, file.path(dir, "seurat_object_useable_cells.rds"))


