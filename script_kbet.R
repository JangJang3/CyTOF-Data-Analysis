################################################################
# Set working environment 
################################################################
current_path <- "/mnt/volumeb/jin/batch_all_smaller"

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

dir_name <- "kBET"

# create data folder where all analysis will be stored
if(!dir.exists(dir_name)) dir.create(dir_name)

# set it as a main directory 
dir <- file.path(getwd(), dir_name)



################################################################
# Using kBET 
################################################################

# already downsampled file with 20,000 cells per samples
srat_obj <- readRDS(str_interp("${current.dir}/seurat_object_UMAP.rds"))

dim(srat_obj)

markers <- rownames(srat_obj)
cell_names <- colnames(srat_obj)

#subset_obj <- subset(srat_obj, subset = batch_group %in% c("Batch_1", "Batch_2"))
#exprs_data <- as.matrix(subset_obj@assays$RNA@layers$data) 

exprs_data <- as.matrix(subset_obj@assays$RNA@layers$scale.data)
exprs_data <- t(exprs_data)

rownames(exprs_data) <- cell_names
colnames(exprs_data) <- markers

# extract the batch labels
batch <- srat_obj@meta.data$batch_group
batch <- gsub("Batch_", "", batch)


require('FNN')

# From Github kBET
# data: a matrix (rows: samples, columns: features (genes))
# k0=floor(mean(table(batch))) #neighbourhood size: mean batch size 
# knn <- get.knn(umap_coords, k=k0, algorithm = 'cover_tree')
  
#data: a matrix (rows: cells or other observations, columns: features (genes); will be transposed if necessary)
#batch: vector or factor with batch label of each cell/observation; length has to match the size of the corresponding data dimension  
#now run kBET with pre-defined nearest neighbours.
batch.estimate <- kBET(exprs_data, batch = batch, k0 = 20,  do.pca = TRUE, dim.pca = 31, plot = TRUE)
saveRDS(file.path(dir, batch.estimate, "kBET_results.rds"))

#batch.estimate <- readRDS(file.path(dir, kBET_results.rds"))


# # Convert matrix to a data frame
umap_coords <- Embeddings(srat_obj, reduction = "umap")
umap_coords_df <- as.data.frame(umap_coords)

# Add kBET p-values to your UMAP data
umap_coords_df$kbet_pvalues <- batch.estimate$results$kBET.pvalue.test  # Add kBET p-values to your UMAP data


ggplot(umap_coords_df, aes(x = umap_1, y = umap_2, color = kbet_pvalues)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "green") +  # Color scale for p-values
  theme_minimal() +
  labs(title = "UMAP with kBET p-values", color = "kBET p-value")
ggsave(file.path(dir, "kBET_UMAP.pdf"), width = 10, heigh = 10)


# Create the expected and observed plot
plot.data <- data.frame(class=rep(c('observed', 'expected'), 
                                  each=length(batch.estimate$stats$kBET.observed)), 
                        data =  c(batch.estimate$stats$kBET.observed,
                                  batch.estimate$stats$kBET.expected))

ggplot(plot.data, aes(class,data)) + geom_boxplot() + 
  labs(x='Test', y='Rejection rate',title='kBET test results') +
  theme_bw() +
  scale_y_continuous(limits=c(0,1))

ggsave(file.path(dir, "kBET_result.pdf"), width = 5, height = 5)
