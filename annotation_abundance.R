################################################################
# Set working environment 
################################################################
current_path <- "/mnt/volumeb/jin/baseline_pipeline/analysis/analysis_only_NK_cells_with_CD3/Analysis_remove_CD3_cluster_6/Analysis_remove_CD3_cluster_6_rerun_0.2"

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


dir_name <- "annotation_abundances"

# create data folder where all analysis will be stored
if(!dir.exists(dir_name)) dir.create(dir_name)

# set it as a main directory 
dir <- file.path(getwd(), dir_name)


################################################################
# Load Seurat Object 
################################################################
print("Status:  Load Seurat Object")

srat_obj <- readRDS("seurat_object_clustering_leiden.rds")


################################################################
# Add  Annotations 
################################################################

# Annotations 
# Cluster_1	Naive CD8⁺ T cells
# Cluster_2	TCM-like CD8⁺ T cells
# Cluster_3	CD57+ CD8⁺ T cells
# Cluster_4	CD57 NK-like CD8⁺ T cells
# Cluster_5	Proliferating CD8⁺ T cells
# Cluster_6	Exhausted CD8⁺ T cells
# Cluster_7	CD161+ NK-like CD8⁺ T cells


Idents(srat_obj) <- "res_0.2"

cluster_ids.new <- c(
  "5" = "CD16dimCD56high NK",
  "4" = "Proliferating NK",
  "7" = "helper ILCs",
  "3" = "CD16highCD56dim CD11b+ NK",
  "8" = "CD16highCD56dim exhausted-like NK",
  "1" = "CD16highCD56dim CD57+ NK",
  "2" = "CD16highCD56dim NK",
  "6" = "CD16highCD56dim ICOS+ NK"
)


names(cluster_ids.new) <- levels(srat_obj)
srat_obj <- RenameIdents(srat_obj, cluster_ids.new)
srat_obj$annotation <- Idents(srat_obj)
saveRDS(srat_obj, file.path(dir, "seurat_object_annotated.rds"))
#srat_obj <- readRDS(file.path(dir, "seurat_object_annotated.rds"))


#number_of_clusters = n_distinct(srat_obj@meta.data["res_0.2"])


################################################################
# Add Colors
################################################################
print("Status:  Add Colors")


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
# Create Annotated Plots 
################################################################
print("Status:  Create Annotated Plots")

# UMAP with annotations
dim_plt <-  DimPlot(srat_obj, reduction = "umap", raster = T, label = F) +
  scale_color_manual(values = colors) +
  ggtitle(str_interp("UMAP Clustering"))

ggsave(file.path(dir, "UMAP_annotation_rastered.pdf"),plot = dim_plt,  width = 10, height = 8)

dim_plt <-  DimPlot(srat_obj, reduction = "umap", raster = F, label = F) +
  scale_color_manual(values = colors) +
  ggtitle(str_interp("UMAP Clustering"))

ggsave(file.path(dir, "UMAP_annotation_editable.pdf"),plot = dim_plt,  width = 10, height = 8)


dim_plt_data <- layer_data(dim_plt)
dim_plt_data <- dim_plt_data[, c("colour", "group")]
dim_plt_data$group <- factor(dim_plt_data$group, levels = sort(unique(dim_plt_data$group), decreasing = FALSE))
dim_plt_df <- unique(dim_plt_data)

color_df <- dim_plt_df %>% arrange(group) %>% mutate(group = levels(srat_obj))
my_colors <- setNames(color_df$colour, color_df$group)


#create histogram
ggplot(data.frame(dim_plt_data), aes(x = group, fill = colour)) +
  geom_bar(position = "dodge", color = "black") +
  theme_minimal() +
  geom_text(
    stat = "count",
    aes(
      label = after_stat(count)
    ),
    position = position_dodge(),
    color = "black",
    size = 3,
    vjust = -0.2
  ) +
  labs(title = "Annotated Histogram",
       x = "Clusters",
       y = "Frequency") +
  scale_fill_identity()

file.name.hist <- str_interp("plot_hist_annotated.pdf")
ggsave(file.path(dir, file.name.hist), width = 10, height = 10, bg= "white")


# Extract the Average Expression of the scaled.data
combined_averages <- AverageExpression(srat_obj, return.seurat = TRUE, group.by = res)

xxx= combined_averages@assays$RNA@layers$data
rownames(xxx) = rownames(combined_averages)
colnames(xxx) = colnames(combined_averages)

# Scale per row (features) and rescale to values -1 and 1
xx = t(scale(t(xxx)))
xx=t(apply(xxx, 1, function(x) scales::rescale(x, to = c(-1, 1))))

clusters <- gsub("g", "Cluster_", colnames(xx))
annotation_col <- data.frame(Cluster = clusters)
rownames(annotation_col) <- clusters
colnames(xx) <- clusters

umap_data <- layer_data(dim_plt)

# Get cluster -> color mapping
cluster_color_df <- unique(umap_data[, c("group", "colour")])
cluster_color_df <- cluster_color_df %>% arrange(group) %>% mutate(group = clusters)
cluster_colors <- setNames(cluster_color_df$colour, cluster_color_df$group)

ann_colors <- list(
  Cluster = cluster_colors
)


# Create matrixplot
cairo_pdf(file.path(cluster_plots.dir, str_interp("matrixplot_annotated.pdf")), width = 8 + number_of_clusters * 0.3, height = 12)
pheatmap::pheatmap(xx, cluster_rows = F, cluster_cols = T, scale = "none", main = str_interp("Average Scaled Expression per Feature with Leiden ${res}"),
                   angle_col = c("90"), annotation_colors = ann_colors, annotation_col = annotation_col,
                   color =rev(RColorBrewer::brewer.pal(n =11, name = "RdBu")))
dev.off()

# bubble heatmap 
# https://stackoverflow.com/questions/76631139/how-to-use-plotmeth-in-pheatmap-title
# file.name.dot <- str_interp("plot_bubble_heatmap_${res}.png")
# 
# dot_plt <- DotPlot(srat_obj, features = markers, scale = TRUE, cluster.idents = FALSE, group.by = res) +
#   scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
#   ggtitle(str_interp("Bubble Heatmap using ${res}")) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
# ggsave(filename = file.path(cluster_plots.dir, file.name.dot), plot = dot_plt, device = "png",  width = 4000, height = 4000,  units = "px", bg = "white")


# heatmap
file.name.heatmap <- str_interp("plot_heatmap_annotated_rastered.png")

heatmap_plt <- DoHeatmap(srat_obj, features = markers, size = 3, slot = "scale.data") +
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")), limits = c(-1, 1), oob=scales::squish) +
  ggtitle("Annotated Heatmap")
ggsave(filename = file.path(cluster_plots.dir, file.name.heatmap), plot = heatmap_plt, device = "png",  width = 6000, height = 4000, units = "px")


file.name.heatmap <- str_interp("plot_heatmap_annotated_editable.pdf")

heatmap_plt <- DoHeatmap(srat_obj, features = markers, size = 3, slot = "scale.data") +
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")), limits = c(-1, 1), oob=scales::squish) +
  ggtitle("Annotated Heatmap")
ggsave(filename = file.path(cluster_plots.dir, file.name.heatmap), plot = heatmap_plt, device = "pdf",  width = 15, height = 10)



################################################################
# Create Abundance Plot and Table 
################################################################
print("Status:  Create Abundance Plot and Table")

controls <- c("HD_ID_1", "HD_ID_3", "HD_ID_4", "HD_ID_5", "HD_ID_6")

# Convert to character to avoid level issues
srat_obj@meta.data$timepoint <- as.character(srat_obj@meta.data$timepoint)

# Perform assignment
srat_obj@meta.data$timepoint[srat_obj@meta.data$patient_id %in% controls] <- "Control"

# Convert back to factor
srat_obj@meta.data$timepoint <- as.factor(srat_obj@meta.data$timepoint)

# change PatID to Pat_ID
srat_obj@meta.data$patient_id <- gsub("PatID", "Pat_ID", srat_obj@meta.data$patient_id)


#cairo_pdf("abundance_plot")
abundance_plot <- dittoBarPlot(
  object = srat_obj,
  var = "res_0.2",         
  group.by = "patient_id",          
  split.by = "timepoint",          
  scale = "percent",               
  main = "Abundance Plot", 
  color.panel =  colors,
  var.labels.rename = levels(srat_obj),
  split.adjust = list(scales = "free_x") 
) 

abundance_plot$layers[[1]]$aes_params$width <- 0.7

cairo_pdf(file.path(dir, "NK_cells_abundance_plot_annotated.pdf"), width = 20)
print(abundance_plot)
dev.off()


# Now extract plot information for the abundance table
abundance_plot <- dittoBarPlot(
  object = srat_obj,
  var = "res_0.2",         
  group.by = "patient_id",          
  split.by = "timepoint",          
  scale = "percent",               
  main = "Abundance Plot", 
  color.panel =  colors,
  #var.labels.rename = cluster_ids.new,
  split.adjust = list(scales = "free_x") 
) 

abundance_data <- abundance_plot$data

abundance_data <- abundance_data %>%
  rename(cluster_id = label) %>%
  rename(sample_id = grouping) %>%
  rename(total_count_per_sample = label.count.total.per.facet) %>% 
  mutate(percent = percent *100) %>% 
  arrange(cluster_id)

annotations <- levels(srat_obj$annotation)
annotations <- setNames(seq_along(annotations), annotations)

# Invert the named vector so you can look up names by cluster_id
id_to_annotation <- names(annotations)[match(abundance_data$cluster_id, annotations)]
abundance_data$annotation <- id_to_annotation

# add Cluster prefix
abundance_data$cluster_id <- paste("Cluster_", abundance_data$cluster_id, sep = "")


clinical_data <- readRDS("/mnt/volumeb/jin/survival_analysis/new_clinical_data_filtered_July.rds")
#clinical_data$Staging_6months[is.na(clinical_data$Staging_6months)] <- "Response"


# Separate dataframes 
short_IDs <- clinical_data[clinical_data$survival == "Short term", ]$patient_id
mid_IDs <- clinical_data[clinical_data$survival == "Mid term", ]$patient_id
long_IDs <- clinical_data[clinical_data$survival == "Long term", ]$patient_id

df_short <- abundance_data[abundance_data$sample_id %in% short_IDs,] 
df_mid <- abundance_data[abundance_data$sample_id %in% mid_IDs,] 
df_long <- abundance_data[abundance_data$sample_id %in% long_IDs,] 

df_control <- abundance_data[abundance_data$sample_id %in% controls,] 

# for writing a data.frame or list of data.frames to an xlsx file
write.xlsx(df_control, file.path(dir, 'NK_cells_abundance_data_control.xlsx'))
write.xlsx(df_short, file.path(dir, 'NK_cells_abundance_data_short_term.xlsx'))
write.xlsx(df_mid, file.path(dir, 'NK_cells_abundance_data_mid_term.xlsx'))
write.xlsx(df_long, file.path(dir, 'NK_cells_abundance_data_long_term.xlsx'))

# ################################################################
# # Spread and Distance
# ################################################################
# print("Status:  Test different min_dist and spreads")
# 
# 
# umap_params <- expand.grid(
#   min_dist = c(0.01, 0.1, 0.3, 0.5, 1.0),
#   spread = c(0.5, 1.0)
# )
# 
# 
# # # How many PCs explain 95% of the variance?
# calculate_npcs <- function(seurat_obj){
#   srat_obj_variance <- seurat_obj@reductions$pca@stdev**2
#   srat_obj_variance <- cumsum(srat_obj_variance)/sum(srat_obj_variance)
#   srat_obj_variance
#   
#   which(srat_obj_variance > 0.95)[1] #
#   
#   needed_pcs <- which(srat_obj_variance > 0.95)[1]
#   
#   print(str_interp("Needed PCs are: ${needed_pcs}"))
#   
#   return(needed_pcs)
# }
# 
# 
# umaps <- lapply(seq(nrow(umap_params)), function(i) {
#   message(i)
#   set.seed(2025)
#   srat_obj <- RunUMAP(srat_obj,
#                       dims = 1:calculate_npcs(srat_obj),
#                       min.dist = umap_params$min_dist[i],
#                       spread = umap_params$spread[i])
#   
#   name <- glue("spread: {umap_params$spread[i]}\nmin_dist: {umap_params$min_dist[i]}")
#   p <- DimPlot(srat_obj, reduction = "umap", label = F) + 
#     scale_color_manual(values = colors) + 
#     theme(legend.position="none") + 
#     ggtitle(name)
#   
#   return(p)
#   #plot_list[[name]] <- p
# })
# 
# # Display all plots in a grid
# plots <- grid.arrange(grobs = umaps, ncol = 3)
# 
# ggsave(file.path(dir, "UMAPs_spreading.png"),plot = plots,  width = 15, height = 15)
