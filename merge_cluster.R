################################################################
# Set working environment 
################################################################
current_path <- "/mnt/volumeb/jin/baseline_pipeline/analysis/analysis_only_Tcells_test_run/Analysis_CD8+_T_cells_rerun"

getwd()
setwd(dir=current_path)


################################################################
# Install and Load Packages 
################################################################

source("/mnt/volumeb/jin/installations.R")


################################################################
# Create File Path and Folders
################################################################
print("Status:  Create File Paths")


dir_name <- "Merged_clusters"

# create data folder where all analysis will be stored
if(!dir.exists(dir_name)) dir.create(dir_name)

# set it as a main directory 
dir <- file.path(getwd(), dir_name)


################################################################
# Merge Clusters 
################################################################
print("Status:  Create File Paths")

srat_obj <- readRDS("seurat_object_clustering_leiden.rds")

# Merge clusters
# Cluster 5, 6 at res 0.3

res <- srat_obj@meta.data$res_0.3

res[res == 6] <- 5

# rename the rest 
res[res == 7] <- 6
res[res == 8] <- 7

srat_obj <- AddMetaData(object = srat_obj  , metadata = res, col.name = "merged_clusters")
saveRDS(srat_obj, file.path(dir, "seurat_object_leiden_merged_clusters.rds"))



################################################################
# Create Plots 
################################################################
print("Status:  Create Plots")

# Save the output in a new directory 
cluster_plots.dir <- file.path(dir, "cluster_plots")
if(!dir.exists(cluster_plots.dir)) dir.create(cluster_plots.dir)


markers <- rownames(srat_obj)

resolutions <- c("merged_clusters")


# save normalized data in data slot (important)
srat_obj[["RNA"]]$data  <- as(object = srat_obj[["RNA"]]$counts, Class = "dgCMatrix")

for(res in resolutions){
  
  number_of_clusters = n_distinct(srat_obj@meta.data[res])
  
  # # UMAP split by clusters 
  dim_plt_split <- DimPlot(srat_obj, reduction = "umap", label = TRUE, cols = DiscretePalette_scCustomize(num_colors = number_of_clusters, palette = "varibow"),
                           group.by = res, split.by = res, ncol = 5) + ggtitle(str_interp("UMAP splitted by clusters"))
  file.name.split <- str_interp("plot_umap_split_leiden_${res}.pdf")
  ggsave(file.path(cluster_plots.dir, file.name.split), plot = dim_plt_split, width = 10, height = 10)
  
  
  # UMAP clustering (editable & rastered)
  dim_plt <- DimPlot(srat_obj, reduction = "umap",
                     cols = DiscretePalette_scCustomize(num_colors = number_of_clusters, palette = "varibow"),
                     label = T, group.by = res) + ggtitle(str_interp("UMAP"))
  file.name.dim <- str_interp("plot_umap_leiden_${res}_rastered.pdf")
  ggsave(file.path(cluster_plots.dir, file.name.dim), width = 10, height = 10)
  
  
  dim_plt <- DimPlot(srat_obj, reduction = "umap",
                     cols = DiscretePalette_scCustomize(num_colors = number_of_clusters, palette = "varibow"),
                     label = T, group.by = res, raster = F) + ggtitle(str_interp("UMAP with Leiden Algorithm and ${res}"))
  file.name.dim <- str_interp("plot_umap_leiden_${res}_editable.pdf")
  ggsave(file.path(cluster_plots.dir, file.name.dim), plot = dim_plt, width = 10, height = 10)
  
  
  # Extract color information from DimPlot
  dim_plt_data <- layer_data(dim_plt)
  dim_plt_data <- dim_plt_data[, c("colour", "group")]
  dim_plt_data$group <- factor(dim_plt_data$group, levels = sort(unique(dim_plt_data$group), decreasing = FALSE))
  
  # Create Histogram 
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
    labs(title = str_interp("Histogram of ${res}"),
         x = "Clusters",
         y = "Frequency") +
    scale_fill_identity() # keep seurat color code
  
  file.name.hist <- str_interp("plot_hist_${res}.pdf")
  ggsave(file.path(cluster_plots.dir, file.name.hist), width = 10, height = 10, bg= "white")
  
  
  # Extract the average expression of the scaled.data
  combined_averages <- AverageExpression(srat_obj, return.seurat = TRUE, group.by = res)
  
  xxx= combined_averages@assays$RNA@layers$data
  rownames(xxx) = rownames(combined_averages)
  colnames(xxx) = colnames(combined_averages)
  
  # scale per row (features) and rescale to values between -1 and 1
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
  cairo_pdf(file.path(cluster_plots.dir, str_interp("matrixplot_${res}.pdf")), width = 8 + number_of_clusters * 0.3, height = 12)
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
  
  
  # Create Complex Heatmap (each line one expression as editable & rastered)
  file.name.heatmap <- str_interp("plot_heatmap_${res}_rastered.png")
  
  heatmap_plt <- DoHeatmap(srat_obj, features = markers, size = 3, slot = "scale.data", group.by = res) +
    scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")), limits = c(-1, 1), oob=scales::squish) +
    ggtitle(str_interp("Heatmap using ${res}"))
  ggsave(filename = file.path(cluster_plots.dir, file.name.heatmap), plot = heatmap_plt, device = "png",  width = 6000, height = 4000, units = "px")
  
  
  file.name.heatmap <- str_interp("plot_heatmap_${res}_editable.pdf")
  
  heatmap_plt <- DoHeatmap(srat_obj, features = markers, size = 3, slot = "scale.data", group.by = res) +
    scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")), limits = c(-1, 1), oob=scales::squish) +
    ggtitle(str_interp("Heatmap using ${res}"))
  ggsave(filename = file.path(cluster_plots.dir, file.name.heatmap), plot = heatmap_plt, device = "pdf",  width = 15, height = 10)
  
}



