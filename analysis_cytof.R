################################################################
# Set working environment 
################################################################
current_path <- "/mnt/volumeb/jin/baseline_pipeline"

getwd()
setwd(dir=current_path)

################################################################
# Installing and Loading
################################################################

source("/mnt/volumeb/jin")

################################################################
# Create File Path and Folders
################################################################
print("Status:  Create File Paths")

dir_name <- "analysis"

# create data folder where all analysis will be stored
if(!dir.exists(dir_name)) dir.create(dir_name)

# set it as a main directory 
dir <- file.path(getwd(), dir_name)
current.dir <- file.path(getwd(), dir_name)


################################################################
# For more than one Batch: Open files 
################################################################
print("Status:  Read Files")


#' Reads the FCS files 
#' 
#' @param dir directory where the data is stored
#' @return a flowSet
#' 
read_fcs <- function(dir = "/mnt/volumeb/jin/Preprocessed_data"){
  
  file.dir <- list.files(dir, pattern="\\.fcs$", recursive = TRUE)
  
  # Check for HD_2_Batch_7 (different colnames, so no merging possible) and remove all other Benny HD_donors except HD_1_Batch_5
  exclude_pattern <- grep("HD_1_Batch_[1-5]|HD_[2]|HD_1_Batch_7", file.dir, invert = TRUE)
  file.dir <- file.dir[exclude_pattern]
  
  fs <- read.flowSet(file.dir, path=dir, transformation = FALSE)
  
  return(fs)
}

#fs <- read_fcs()
#saveRDS(fs, file = file.path(dir, "fs.rds"))
#fs <- readRDS(file.path(dir, "fs.rds"))

# Names of the files (currently 54 samples)
#sampleNames(fs)

# Column names: abbreviation/dyes/chemicals used
#colnames(fs)

# marker expression data for first flow set
#exprs(fs[[1]])

# expression data for the first flow set and its first marker
#exprs(fs[[1]])[, 1]


# ################################################################
# # Create Metadata file 
# ################################################################
print("Status:  Create Meta File")

# create a metadata frame with
# file name
# condition healthy or HCC
# sample id   Pat_ID_1.0.1
# timepoint T0
# patient id Pat_ID_1
# batch group


#' help function to capture information from sample ID
#' 
#' @param sample_ids list of all sample_ids
#' @return list with batch group, condidtions, patient ids, and batch number
extract_info <- function(sample_ids){
  conditions <- c()
  timepoints <- c()
  patient_ids <- c()
  batch_groups <- c()
  
  for(id in sample_ids){
    if (grepl("Batch_", id)){
      batch_number <- sub(".*Batch_(\\d+).*", "\\1", id)
      batch_groups <- c(batch_groups, paste("Batch_", batch_number, sep = ""))
      conditions <- c(conditions, "Healthy")
      timepoints <- c(timepoints, "T0")
      patient_id <- sub(".*HD_(\\d+).*", "\\1", id)
      patient_ids <- c(patient_ids, paste("HD_ID_", patient_id, sep = ""))
    } else {
      # extract the last number after the last dot (e.g. 1.0.1, extract 1)
      batch_groups <- c(batch_groups, sub(".*(\\d+)\\.*$", "Batch_\\1", id))
      conditions <- c(conditions, "HCC")
      timepoints <- c(timepoints, sub(".*(\\d+)\\.(\\d+)\\.(\\d+)\\.*", "T\\2", id)) # captures the middle number
      patient_ids <- c(patient_ids, sub("^(.*?)(\\..*)?$", "\\1", id)) # captures everything before and including the number before the first dot
    }
  }
  
  return(list(conditions, timepoints, patient_ids, batch_groups))
}



#' creates a metadata file from the FCS files
#' @param fcs.files 
#' @param save_file boolean
#' @param file_dir 
#' @return md file 

create_md <- function(fcs.files, save_file = FALSE, file_dir = dir){
  file_name <- sampleNames(fcs.files)
  sample_id <- sub("^Norm_(.*)_CC_gated\\.fcs$", "\\1", sampleNames(fcs.files)) # returns as 'HD_PBMC_Batch_1'
  info <- extract_info(sample_id)
  condition <- info[[1]]
  patient_id <-  info[[3]]
  timepoint <- info[[2]]
  batch_group <- info[[4]]
  
  md <- data.frame(file_name, sample_id, condition, timepoint, patient_id, batch_group)
  
  if(save_file){
    write.csv(md, file = file.path(file_dir, "metadata_all.csv"), row.names = FALSE)
  }
  
  return(md)
}

# metatable as md
md <- create_md(fs, save_file = TRUE)

# load md file

#md <- read.csv(file.path(dir, "metadata_all.csv"))


# ################################################################
# #  Create Panel file 
# ################################################################
print("Status:  Create Panel")

# Create a panel dataframe that contains the marker names, the antigen names and the marker class information
# You can either save the file or not
# Input: reference file with the marker names you want to use, flowset 

create_panel <- function(marker_xml, fs, save_file = FALSE, file_dir = dir){
  panel <- pData(parameters(fs[[2]]))
  panel <- panel %>%
    select(name, desc) 
  
  colnames(panel) <- c("fcs_colname", "marker_name") # rename columes
  
  # remove "Di" and rearrange the order of the elements, add new column for linking (thus we wont change original column)
  panel <- panel %>%
    mutate(fcs_link = str_remove(fcs_colname, "Di")) %>%  # Remove "Di"
    mutate(fcs_link = str_replace(fcs_link, "(\\D+)(\\d+)", "\\2\\1"))  # Rearrange so number comes before letters
  
  # ensure case-insensitivity by converting both columns to lowercase
  marker_xml <- marker_xml %>%
    mutate(`marker_link`= tolower(str_trim(`Isotype/fluorescein`)))
  
  # ensure case-insensitivity by converting to lowercase 
  panel <- panel %>%
    mutate(fcs_link = tolower(str_trim(fcs_link)))
  
  panel <- panel %>%
    left_join(marker_xml, by = c("fcs_link" = "marker_link")) %>%
    filter( !is.na(Marker)) %>%
    select(-marker_name, -fcs_link, -`Isotype/fluorescein`) # Drop the columns that are no longer needed
  
  colnames(panel)[2] <- c("antigen")
  
  # functional markers: FoxP3 and Ki-67
  # the rest are lineage markers
  
  lineage_markers <- c("Ki-67", "FoxP3")
  
  # differentiate between marker class 
  panel <- panel %>%
    mutate(marker_class = ifelse(antigen %in% lineage_markers, "state", "type"))
  
  # remove CD45, because every other marker is CD45+
  panel <- panel %>%
    filter(antigen != "CD45")
  
  if(save_file){
    write.csv(panel,file = file.path(file_dir,"panel_all.csv"), row.names = FALSE)
  }
  
  return(panel)
}

# Open marker information from excel sheet
marker_xml <- openxlsx::read.xlsx("/mnt/volumeb/jin/2024_06_10_HCC Study CYTOF Panel.xlsx", rows = 1:33, cols = 1:2, sheet = 1)

#head(marker_xml)

panel <- create_panel(marker_xml, fs, save.file = TRUE)

# read from panel file
#panel <- read.csv(file.path(dir, "panel_all.csv"))

# spot check that all panel columns are in the flowSet object
#all(panel$fcs_colname %in% colnames(fs))


# ################################################################
# # Pre-processing 
# ################################################################
print("Status:  Preprocessing")

transform_data <- function(fcs_file, md, save_file = FALSE, file_dir = dir){
  
  md$condition <- factor(md$condition, levels = c("Healthy", "HCC"))
  md$timepoint <- factor(md$timepoint, levels = c("T0", "T1"))
  md$batch_group <- factor(md$batch_group, levels = c("Batch_1", "Batch_2", "Batch_3", "Batch_4", "Batch_5", "Batch_6", "Batch_7"))
  md$patient_id <- factor(md$patient_id, levels = unique(md$patient_id))
  md$sample_id <- factor(md$sample_id,
                         levels = md$sample_id[order(md$condition)])
  # specify md_cols
  factors <- list(factors = c("patient_id", "sample_id", "condition", "timepoint", "batch_group"))
  
  # construct SingleCellExperiment
  sce <- prepData(fcs_file, panel, md, features = panel$fcs_colname, md_cols = factors, transform = TRUE, cofactor = 5) # saves only the data
  
  if(save_file){
    saveRDS(sce, file = file.path(file_dir, "sce_all.rds"))
  }
  return(sce)
}

sce <- transform_data(fs, md, save.file = TRUE)

# load sce
#sce <- readRDS(file.path(dir, "sce_all.rds"))

# view number of cells per sample
#table(sce$sample_id)

#n_cells(sce) # alternative methode for number of cells per sample

plot_count_cells <- function(data, file_dir = dir, show = FALSE){
  plt <- plotCounts(data, group_by = "sample_id", color_by = "condition") + labs(fill = "Condition")
  
  file_name <- str_interp("${file_dir}/plot_count.png")
  ggsave(filename = file_name, plot = plt, width = 5000, height = 4000, units = "px")
  
  if(show){
    print(plt)
  }
}

plot_count_cells(sce)

# number of cells after preprocessing
#dim(sce)

# view non-mass channels
#names(int_colData(sce))



# ################################################################
# Use only the useable cell from batch 7 
# ################################################################

# batch7.useable.cells <- readRDS("/mnt/volumeb/jin/other_workflow/based_paper/analysis_with_transformation/remove_junk_cells/seurat_object_useable_cells.rds")
# batch7.useable.cells.indices <- colnames(batch7.useable.cells)


# ################################################################
# Downsampling/subsampling
# ################################################################
# Goal: Reduce the number of cells but keep biological significance 
# Downsample 100 000 cells from each sample

print("Status:  Subsampling")

subsampling_data <- function(data, number = 100000, save.files = FALSE, file_dir = dir, file.name = "subsample_data",  seed = 2025){
  set.seed(seed)
  
  # dim(assay(sce, "exprs")) is 31 x #cells 
  
  expr_matrix <- assay(data, "exprs")
  
  # returns metadata of each cell
  meta_data <- colData(data)
  #dim(meta_data) is  #cells x 5
  
  # add missing rownames for better information extraction
  rownames(meta_data) <- 1:nrow(colData(data))
  colnames(expr_matrix) <- rownames(meta_data)
  
  
  # filter the metadata by condition
  healthy_cells <- meta_data[meta_data$condition == "Healthy", ]
  hcc_cells <- meta_data[meta_data$condition == "HCC", ]
  
  # filter the names of the PatIDs by condition
  healthy_IDs <- md$sample_id[grepl("HD_", md$sample_id)]
  hcc_IDs <- md$sample_id[grepl("Pat_ID", md$sample_id)]
  
  
  # only include those cells that are useable (no junk cells)
  # meaning only keep cells that are included by batch7_useable_cells_indices
  
  healthy_cells.batch7 <- healthy_cells[healthy_cells$batch_group == "Batch_7", ]
  healthy_cells.batch7 <- healthy_cells.batch7[rownames(healthy_cells.batch7) %in% batch7.useable.cells.indices, ]
  
  hcc_cells.batch7 <- hcc_cells[hcc_cells$batch_group == "Batch_7", ]
  hcc_cells.batch7 <- hcc_cells.batch7[rownames(hcc_cells.batch7) %in% batch7.useable.cells.indices, ]
  
  
  
  # group the data
  # select from each group 100,000 cells
  # but check if the cell amount exceed 100,000 cells
  # if yes, then sample 100, 000 cells
  # if not, take everything (without random sampling)
  # data.frame(healthy_cells) %>% group_by(sample_id) %>% slice_sample(n= 100000)
  # works but check if the sample is randomized
  
  healthy_cells <- data.frame(healthy_cells) %>%  
    rownames_to_column(var="row_IDs") %>% 
    filter(batch_group != "Batch_7") %>% 
    group_by(sample_id) %>% 
    slice_sample(n= number)
  
  hcc_cells <- data.frame(hcc_cells) %>%  
    rownames_to_column(var="row_IDs") %>% 
    filter(batch_group != "Batch_7") %>%
    group_by(sample_id) %>% 
    slice_sample(n= number)
  
  # for batch 7 only 
  healthy_cells.batch7 <- data.frame(healthy_cells.batch7) %>%  
    rownames_to_column(var="row_IDs") %>% 
    group_by(sample_id) %>% 
    slice_sample(n= number)
  
  hcc_cells.batch7 <- data.frame(hcc_cells.batch7) %>%  
    rownames_to_column(var="row_IDs") %>% 
    group_by(sample_id) %>% 
    slice_sample(n= number)
  
  # extract the exprs info for the selected cells
  indices_sample <- c(healthy_cells$row_IDs, healthy_cells.batch7$row_IDs, hcc_cells$row_IDs, hcc_cells.batch7$row_IDs)
  indices_sample.ordered <- as.character(sort(as.numeric(indices_sample)))
  
  # Seurat’s convention is that the cells are columns and the features are rows.
  # Dimension would be 32 5,313,657 with markers as rows and cells as columns
  expr_matrix <- expr_matrix[, colnames(expr_matrix) %in% indices_sample.ordered ]
  
  # To add meta data your row names are cells and the corresponding meta data information should be the row values/column names
  dim(meta_data) # 18729297 x 5 is correct, no transpose is needed
  meta_data <- meta_data[colnames(expr_matrix), ]
  
  if(save_files){
    #file.name.indices <- str_interp("${file.dir}/subsample_indices.rds")
    #saveRDS(indices_sample.ordered, file = file.name.indices)
    
    #file.name.meta <- str_interp("${file.dir}/subsample_meta_data.rds")
    #saveRDS(meta_data, file = file.name.meta)
    
    #file.name.exprs <- str_interp("${file.dir}/subsample_exprs.rds")
    #saveRDS(expr_matrix, file = file.name.exprs)
    
    subsample.data <- list(indices_sample.ordered, meta_data, expr_matrix)
    file.name <- str_interp("${file.name}.rds")
    saveRDS(subsample.data, file = file.path(file_dir, file.name))
    
  }
  
  return(subsample.data)
}

#subsample.data <- subsampling_data(sce, number = 80000, save.files = TRUE, file.name = "subsample_data_80000")
#subsample.data <- subsampling_data(sce, number = 90000, save.files = TRUE, file.name = "subsample_data_90000")
#subsample.data <- subsampling_data(sce, number = 100000, save.files = TRUE, file.name = "subsample_data_100000")

# subsample.data <- subsampling_data(sce, number = 100000, save.files = TRUE)
# subsample.data <- readRDS(file.path(dir, "subsample_data_90000.rds"))

# ################################################################
# Downsampling/subsampling
# ################################################################
# Goal: Reduce the number of cells but keep biological significance 
# Downsample 100 000 cells from each sample

print("Status:  Subsampling")

subsampling_data <- function(data, number = 100000, save_files = FALSE, file_dir = dir, file.name = "subsample_data", seed = 2025){
  set.seed(seed)
  
  # dim(assay(sce, "exprs")) is 31 x #cells 
  
  count_matrix <- assay(data, "counts")
  expr_matrix <- assay(data, "exprs")
  
  # returns metadata of each cell
  meta_data <- colData(data)
  #dim(meta_data) is  #cells x 5
  
  # add missing rownames for better information extraction
  rownames(meta_data) <- 1:nrow(colData(data))
  colnames(count_matrix) <- rownames(meta_data)
  colnames(expr_matrix) <- rownames(meta_data)
  
  
  # filter the subsamples with the correct conditions
  healthy_cells <- meta_data[meta_data$condition == "Healthy", ]
  hcc_cells <- meta_data[meta_data$condition == "HCC", ]
  
  # subsample the data, meaning for each sample filter 100.000 cells
  healthy_IDs <- md$sample_id[grepl("HD_", md$sample_id)]
  hcc_IDs <- md$sample_id[grepl("Pat_ID", md$sample_id)]
  
  
  # group the data
  # select from each group 100,000 cells
  # but check if the cell amount exceed 100,000 cells
  # if yes, then sample 100, 000 cells
  # if not, take everything (without random sampling)
  # data.frame(healthy_cells) %>% group_by(sample_id) %>% slice_sample(n= 100000)
  # works but check if the sample is randomized
  
  healthy_cells <- data.frame(healthy_cells) %>%  rownames_to_column(var="row_IDs") %>% group_by(sample_id) %>% slice_sample(n= number)
  hcc_cells <- data.frame(hcc_cells) %>%  rownames_to_column(var="row_IDs") %>% group_by(sample_id) %>% slice_sample(n= number)
  
  # extract the exprs info for the selected cells
  indices_sample <- c(healthy_cells$row_IDs, hcc_cells$row_IDs)
  indices_sample.ordered <- as.character(sort(as.numeric(indices_sample)))
  
  # Seurat’s convention is that the cells are columns and the features are rows.
  # Dimension would be 32 5,313,657 with markers as rows and cells as columns
  count_matrix <- count_matrix[, colnames(count_matrix) %in% indices_sample.ordered ]
  expr_matrix <- expr_matrix[, colnames(expr_matrix) %in% indices_sample.ordered ]
  
  # To add meta data your row names are cells and the corresponding meta data information should be the row values/column names
  dim(meta_data) # 18729297 x 5 is correct, no transpose is needed
  meta_data <- meta_data[colnames(expr_matrix), ]
  
  
  indices_sample.ordered <- rownames(meta_data)
  
  # short test
  sum(colnames(expr_matrix) == rownames(meta_data))
  
  if(save_files){
    #file.name.indices <- str_interp("${file.dir}/subsample_indices.rds")
    #saveRDS(indices_sample.ordered, file = file.name.indices)
    
    #file.name.meta <- str_interp("${file.dir}/subsample_meta_data.rds")
    #saveRDS(meta_data, file = file.name.meta)
    
    #file.name.exprs <- str_interp("${file.dir}/subsample_exprs.rds")
    #saveRDS(expr_matrix, file = file.name.exprs)
    
    subsample.data <- list(indices_sample.ordered, meta_data, count_matrix, expr_matrix)
    file.name <- str_interp("${file.name}.rds")
    saveRDS(subsample.data, file = file.path(file_dir, file.name))
    
  }
  
  return(subsample.data)
}

#subsample.data <- subsampling_data(sce, number = 80000, save.files = TRUE, file.name = "subsample_data_80000")
subsample.data <- subsampling_data(sce, number = 90000, save.files = TRUE, file.name = "subsample_data_90000")
#subsample.data <- subsampling_data(sce, number = 100000, save.files = TRUE, file.name = "subsample_data_100000")



# ################################################################
# Create Seurat Object and Scaling 
# ################################################################
print("Status:  Create Seurat Object")

initiate_seurat_object <- function(subsample.data, save.file = FALSE, file.dir = current.dir){
  count_matrix <-subsample.data[[3]]
  expr_matrix <- subsample.data[[4]]
  meta_data   <- subsample.data[[2]]
  
  seurat_obj  <- Seurat::CreateSeuratObject(counts = count_matrix, data = expr_matrix, meta.data = meta_data)
  
  # change Anti-Fitc and subset cell specific markers
  rownames(seurat_obj)[28] <- "TCR Vδ1 "
  
  # only cell-specific 
  marker_ref <- openxlsx::read.xlsx(file.path("/mnt/volumeb/jin/baseline_pipeline", "cell_specific_markers.xlsx"), rows = 1:16, cols = 1:2, sheet = 1)
  
  rownames(seurat_obj) %in% marker_ref$Marker
  
  seurat_obj <- subset(seurat_obj, features = marker_ref$Marker)
  
  # Add meta data information
  seurat_obj  <- AddMetaData(object = seurat_obj  , metadata = meta_data$patient_id, col.name = "patient_id")
  seurat_obj  <- AddMetaData(object = seurat_obj  , metadata = meta_data$sample_id, col.name = "sample_id")
  seurat_obj  <- AddMetaData(object = seurat_obj  , metadata = meta_data$condition, col.name = "condition")
  seurat_obj  <- AddMetaData(object = seurat_obj  , metadata = meta_data$timepoint, col.name = "timepoint")
  seurat_obj  <- AddMetaData(object = seurat_obj  , metadata = meta_data$batch_group, col.name = "batch_group")
  
  
  seurat_obj  <- Seurat::FindVariableFeatures(seurat_obj , nfeatures = nrow(seurat_obj))
  
  all.markers <- rownames(seurat_obj)
  
  # scale data using z-scoring, score each feature individual (not overall)
  seurat_obj  <- ScaleData(object = seurat_obj  , do.scale = TRUE, do.center = TRUE, layer = "counts", features = all.markers) 
  
  if(save.file){
    saveRDS(seurat_obj, file = file.path(dir, "seurat_object_init.rds"))
  }
  return(seurat_obj)
}

#srat_obj <- initiate_seurat_object(subsample.data, save.file = TRUE)

# load seurat object
#srat_obj <- readRDS(file.path(dir, "seurat_object_init.rds"))


################################################################
# Create RidgePlot 
################################################################
print("Status:  Create ridge plot")

#Create ridge plot with quantile cutoff of q90

create_ridge_plot <- function(seurat_obj, file.dir = current.dir, assay = "scale.data", name = "plot_ridge", group.plt = "condition", width = 8000, height = 4000){
  markers <- rownames(seurat_obj)
  rownames(seurat_obj@assays$RNA@layers$scale.data) <- rownames(seurat_obj@assays$RNA)
  
  plot_lst <- list()
  
  if(group.plt == "condition"){
    
    for (marker in markers){
      q90 <- quantile(seurat_obj@assays$RNA@layers$scale.data[marker, ], 0.90, na.rm = TRUE)
      p <- RidgePlot(seurat_obj, features = marker, group.by = group.plt, layer = assay, cols = c("#7FDF1F", "#EB2B2B")) &
        coord_cartesian(xlim = c(NA, q90)) &
        theme(legend.position = "right", axis.title.y = element_blank())
      plot_lst[[marker]] <- p
    }
    
    
  } else {
    
    for (marker in markers){
      q90 <- quantile(seurat_obj@assays$RNA@layers$scale.data[marker, ], 0.90, na.rm = TRUE)
      p <- RidgePlot(seurat_obj, features = marker, group.by = group.plt, layer = assay) &
        coord_cartesian(xlim = c(NA, q90)) &
        theme(legend.position = "right", axis.title.y = element_blank())
      plot_lst[[marker]] <- p
    }
    
  }
  
  final_ridge_plot <- wrap_plots(plot_lst, ncol = 4, nrow = 8)
  
  plt.name <- str_interp("${name}.png")
  ggsave(filename = file.path(dir, plt.name), plot = final_ridge_plot, width = width, height = height, units = "px")
}

create_ridge_plot(srat_obj)


################################################################
# Perform PCA 
################################################################
print("Status:  Perform PCA")

# create data folder where all plots will be stored
plot.dir <- file.path(dir, "plots")
if(!dir.exists(plot.dir)) dir.create(plot.dir)

perform_pca <- function(seurat_obj, number_pcs = 15, file.dir = current.dir, save.data = FALSE, name = "seurat_object_pca"){
  seurat_obj <- RunPCA(seurat_obj , npcs = number_pcs)
  
  # Examine and visualize PCA results a few different ways
  #print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 10)
  
  VizDimLoadings(seurat_obj, dims = 1:2, reduction = "pca")
  ggsave(filename = file.path(plot.dir, "plot_vizDimLoading.png"), width = 4000, height = 4000, units = "px")
  
  DimPlot(seurat_obj , reduction = "pca", raster = FALSE)
  ggsave(filename = file.path(plot.dir, "plot_pca.png"), width = 4000, height = 4000, units = "px")
  
  DimPlot(seurat_obj, reduction = "pca", group.by = "batch_group", raster = FALSE)
  ggsave(filename = file.path(plot.dir, "plot_pca_batch.png"), width = 4000, height = 4000, units = "px")
  
  DimPlot(seurat_obj , reduction = "pca", group.by = c("condition", "batch_group"), raster = FALSE)
  ggsave(filename = file.path(plot.dir, "plot_pca_condition_batch.png"), width = 6000, height = 4000, units = "px")
  
  # Alternative to Jack Straw
  ElbowPlot(seurat_obj , ndims = number_pcs)
  ggsave(filename = file.path(plot.dir, "plot_elbow.png"), width = 4000, height = 4000, units = "px")
  
  if(save.data){
    file.name <- str_interp("${name}.rds")
    saveRDS(seurat_obj, file = file.path(dir, file.name))
  }
  
  return(seurat_obj)
}

srat_obj <- perform_pca(srat_obj, save.data = TRUE)
#srat_obj <- readRDS("/mnt/volumeb/jin/batch_all_hd_hcc_smaller/seurat_object_pca.rds")

#srat_obj_no_HD_2 <- perform_pca(srat_obj_no_HD_2, save.data = TRUE, name = "seurat_object_pca_no_HD_2")
#srat_obj_no_HD_2 <- readRDS("/mnt/volumeb/jin/batch_all_hd_hcc/seurat_object_pca_no_HD_2.rds")


# # How many PCs explain 95% of the variance?
calculate_npcs <- function(seurat_obj){
  srat_obj_variance <- seurat_obj@reductions$pca@stdev**2
  srat_obj_variance <- cumsum(srat_obj_variance)/sum(srat_obj_variance)
  srat_obj_variance
  
  which(srat_obj_variance > 0.95)[1] #
  
  needed_pcs <- which(srat_obj_variance > 0.95)[1]
  
  print(str_interp("Needed PCs are: ${needed_pcs}"))
  
  return(needed_pcs)
}

# Suggesting that the majority of signal are in the first X PCs, with little more signal in the next Y PCs.
# Meaning that a total of X PCs explain the variance in the data. 
#npcs <- calculate_npcs(srat_obj)


################################################################
# Perform UMAP
################################################################
print("Status:  Create UMAP")

# create data folder where all plots will be stored
plot.dir <- file.path(dir, "plots")
if(!dir.exists(plot.dir)) dir.create(plot.dir)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
generate_UMAPs <- function(seurat_obj, file.dir = current.dir, save.data = FALSE){
  seurat_obj  <- RunUMAP(seurat_obj, dims = 1:calculate_npcs(seurat_obj))
  
  if(save.data){
    saveRDS(seurat_obj, file.path(dir, "seurat_object_UMAP.rds"))
  }
  
  DimPlot(seurat_obj, reduction = "umap", raster = FALSE)
  ggsave(file.path(plot.dir, "plot_umap.png"), width = 4000, height = 4000, units = "px")
  
  DimPlot(seurat_obj, reduction = "umap", group.by = "batch_group", raster = FALSE) + facet_wrap(~ seurat_obj@meta.data$condition)
  ggsave(file.path(plot.dir, "plot_umap_condition_batch.png"), width = 6000, height = 4000, units = "px")
  
  DimPlot(seurat_obj, reduction = "umap", group.by = "condition", raster = FALSE, split.by = "condition")
  ggsave(file.path(plot.dir, "plot_umap_splitby_condition.png"), width = 6000, height = 4000, units = "px")
  
  
  DimPlot(seurat_obj, reduction = "umap", group.by = "batch_group", raster = FALSE, split.by = "batch_group", ncol = 3)
  ggsave(file.path(plot.dir, "plot_umap_splitby_batch.png"), width = 8000, height = 6000, units = "px")
  
  DimPlot(seurat_obj, reduction = "umap", group.by = "sample_id", raster = FALSE, split.by = "sample_id", ncol = 6)
  ggsave(file.path(plot.dir, "plot_umap_splitby_samples.png"), width = 10000, height = 10000, units = "px")
  
  
  return(seurat_obj)
}

srat_obj <- generate_UMAPs(srat_obj, save.data = TRUE)
#srat_obj <- readRDS(file.path(dir, "seurat_object_UMAP.rds"))


# ################################################################
# # Create Feature Plots
# ################################################################
print("Status:  Create Feature Plots")

# create data folder where all plots will be stored
plot.dir <- file.path(dir, "plots2")
if(!dir.exists(plot.dir)) dir.create(plot.dir)


# adjust the marker scores for the specific markers
adjust_markers <- c(
  "TCR Vδ2 " = 1,
  "CD274 (PD-L1)" = 0.5,
  "CD25 (IL-2R)" = 1,
  "CD223 (LAG-3)" = 0.5,
  "TCRγδ" = 1,
  "TCR Vα7.2 " = 1.5,
  "CD279 (PD-1)" = 1,
  "CD86" = 1.5,
  "FoxP3" = 0.5,
  "CD152 (CTLA-4)" = 0.5,
  "Ki-67" = 1,
  "Anti-FITC" = 1,
  "CD278/ICOS" = 0.5
)

# returns the keys
# names(adjust_markers)

create_feature_plots <- function(seurat_obj, file.dir = current.dir , name = "UMAP_feature_plot", assay = "scale.data",  adjust = FALSE, adjust.markers = adjusted_markers){
  markers <- rownames(seurat_obj)
  
  if(adjust){
    
    feature_lst <- c()
    for(marker in markers){
      # check if we need to adjust the marker
      if(marker %in% names(adjust.markers)){
        # adjust marker with the value from adjust_markers vector
        plot <- FeaturePlot(seurat_obj, features = marker, max.cutoff = adjust.markers[marker], split.by = "condition", slot = assay, raster = FALSE, ncol = 2) & labs(color = "z-score") & 
          theme(legend.position = "right") & 
          scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
      } else {
        # if not, then just do the regular setting
        plot <- FeaturePlot(seurat_obj, features = marker, max.cutoff = "q90", split.by = "condition", slot = assay, raster = FALSE, ncol = 2) & labs(color = "z-score") & 
          theme(legend.position = "right") & 
          scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
      }
      feature_lst[[marker]] <- plot
    }
    
    
    # plot each feature plots in sets of 4
    for(i in c(4, 8, 12, 16)){
      file.name <- str_interp("${name}_${i/4}.png")
      if(i == 32){
        final_feature_plot <- wrap_plots(feature_lst[(i-3):length(markers)], ncol = 1, nrow = 3)
      } else {
        final_feature_plot <- wrap_plots(feature_lst[(i-3):i], ncol = 1, nrow = 4)
      }
      ggsave(filename = file.path(plot.dir, file.name), width = 4000, height = 4000, units = "px")
    }
    
    
  } else {
    for(i in c(4, 8, 12, 16)){
      file.name <- str_interp("${name}_${i/4}_rastered.png")
      plt <- FeaturePlot(seurat_obj, features = markers[(i-3):i], max.cutoff = "q90", slot = assay, raster = T) & 
        scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")), limits = c(-1, 1), oob=scales::squish) & theme(legend.position = "right")
      ggsave(filename = file.path(plot.dir, file.name), width = 4000, height = 4000, units = "px")
      
      # editable version
      file.name <- str_interp("${name}_${i/4}_editable.pdf")
      plt <- FeaturePlot(seurat_obj, features = markers[(i-3):i], max.cutoff = "q90", slot = assay, raster = F) & 
        scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")), limits = c(-1, 1), oob=scales::squish) & theme(legend.position = "right")
      ggsave(filename = file.path(plot.dir, file.name), plot = plt,  width = 10, height = 10)
    }
  }
}

# create feature plots 
create_feature_plots(srat_obj)
# create_feature_plots(srat_obj, adjust = TRUE, adjust.markers = adjust_markers)



################################################################
# Write Data on Disc 
################################################################
print("Status:  On Disc")

compress_seurat_obj <- function(seurat_obj, file.dir = current.dir){
  
  # Embeddings(srat_obj, reduction="pca") %>% is.na() %>% sum
  # Embeddings(srat_obj, reduction="pca") %>% is.nan() %>% sum
  # 
  # Embeddings(srat_obj, reduction="umap") %>% is.na() %>% sum
  # Embeddings(srat_obj, reduction="umap") %>% is.nan() %>% sum
  
  # 
  rownames <- rownames(seurat_obj)
  colnames <- colnames(seurat_obj)
  # 
  # #https://satijalab.org/seurat/articles/seurat5_bpcells_interaction_vignette
  # # reduce memory usage by storing all scale.data as sparse matrix 
  seurat_obj[["RNA"]]@layers$counts <- as(seurat_obj[["RNA"]]@layers$counts,  Class = "dgCMatrix")
  seurat_obj[["RNA"]]@layers$data<- as(seurat_obj[["RNA"]]@layers$data,  Class = "dgCMatrix")
  seurat_obj[["RNA"]]@layers$scale.data<- as(seurat_obj[["RNA"]]@layers$scale.data,  Class = "dgCMatrix")
  
  rownames(seurat_obj) <- rownames
  colnames(seurat_obj) <- colnames
  
  # Write the counts layer to a directory
  write_matrix_dir(mat = seurat_obj[["RNA"]]@layers$counts, dir = file.path(dir, "seurat_exprs"))
  counts.mat <- open_matrix_dir(dir = file.path(dir, "seurat_exprs"))
  
  write_matrix_dir(mat = seurat_obj[["RNA"]]@layers$data, dir = file.path(dir, "seurat_scale_data"))
  data.mat <- open_matrix_dir(dir = file.path(dir, "seurat_data"))
  
  
  write_matrix_dir(mat = seurat_obj[["RNA"]]@layers$scale.data, dir = file.path(dir, "seurat_scale_data"))
  scale.data.mat <- open_matrix_dir(dir = file.path(dir, "seurat_scale_data"))
  
  seurat_obj[["RNA"]]@layers$counts <- counts.mat
  seurat_obj[["RNA"]]@layers$data <- data.mat
  seurat_obj[["RNA"]]@layers$scale.data <- scale.data.mat
  rownames(seurat_obj)
  
  saveRDS(seurat_obj, file.path(dir, "seurat_object_compressed.rds"))
  
  return(seurat_obj)
}

srat_obj <- compress_seurat_obj(srat_obj)
#srat_obj <- readRDS(file.path(dir,"seurat_object_compressed.rds"))


################################################################
# Perform Clustering 
################################################################

# create data folder where all analysis will be stored

# rm(list = setdiff(ls(), c("srat_obj", "current.dir", "dir", "calculate_npcs")))
# gc()
# #plan("multisession", workers = 25)
# options(future.globals.maxSize = 900 * 1024^3) # 900 GB
# 
set.seed(2025)
# # calculate_npcs is 23 (default)
srat_obj <- FindNeighbors(srat_obj, dims = 1:calculate_npcs(srat_obj), reduction = "pca", verbose = TRUE)
saveRDS(srat_obj, file = file.path(dir, "seurat_object_SNN.rds"))

#srat_obj <- readRDS(file.path(dir, "seurat_object_SNN.rds"))

################################################################
# Perform Clustering 
################################################################
print("Status:  Perform Clustering")


perform_clustering  <- function(seurat_obj, npcs = calculate_npcs(seurat_obj), file.dir = current.dir, file.name = "leiden", seed = 123, resolutions = c(0.2, 0.5, 0.8, 1.0), algo = 1){
  # The following is basically PhenoGraph.
  # increase maxSize
  # options(future.globals.maxSize = 1000 * 1024^2) 
  # https://satijalab.org/seurat/archive/v3.0/future_vignette.html 
  
  options(future.globals.maxSize = 900 * 1024^3)  # 700 GB
  
  plan("multisession", workers = 35)
  
  set.seed(seed)
  
  # This will return a list of seurat objects with different resolutions 
  srat_obj_clusters <- future_lapply(resolutions, function(res) {
    tryCatch({
      message(paste("Running clustering for resolution:", res))
      srat_obj_tmp <- FindClusters(seurat_obj, resolution = res, algorithm = algo, random.seed = 123)
      
      # Extract the clustering results and add them to the original Seurat object
      cluster_column_name <- paste0("res_", res)
      clustering_result <- srat_obj_tmp@meta.data$seurat_clusters
      
      return(list(resolution = res, clusters = clustering_result))
    }, error = function(e) {
      message(paste("Error in resolution", res, ":", e$message))
    })
  }, future.seed = 42L)
  
  file.name.list <- str_interp("seurat_object_clustering_${file.name}_list.rds")
  saveRDS(srat_obj_clusters,file.path(dir, file.name.list))
  
  return(srat_obj_clusters)
}

rm(list = setdiff(ls(), c("srat_obj", "dir", "calculate_npcs", "perform_clustering")))
gc()
# default leiden algorithm
# srat_obj_clusters <- perform_clustering(srat_obj, calculate_npcs(srat_obj), file.name = "louvain", algo = 1 , resolutions = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6))

# Default resolution for Leiden is 1.0
#srat_obj_clusters <- perform_clustering(srat_obj, calculate_npcs(srat_obj), algo = 4 , resolutions = c(1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2.0))
srat_obj_clusters <- perform_clustering(srat_obj, calculate_npcs(srat_obj), algo = 4 , resolutions = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1.0, 1.1))

srat_obj_clusters <- readRDS(str_interp("${current.dir}/seurat_object_clustering_leiden_list.rds"))

plan("sequential")
for (result in srat_obj_clusters) {
   if (!is.null(result)) {
     # Add the clustering result to the original srat_obj metadata
     cluster_column_name <- str_interp("res_${result$resolution}")
     srat_obj <- AddMetaData(srat_obj, metadata = result$clusters, col.name = cluster_column_name)
  }
 }

saveRDS(srat_obj, file.path(dir, "seurat_object_clustering_leiden.rds"))


################################################################
# Create Heatmaps
################################################################
print("Status:  Create Heatmaps")

srat_obj <- readRDS(file.path(dir, "seurat_object_clustering_leiden.rds"))

markers <- rownames(srat_obj)

# UMAPs with clustering
resolutions <- colnames(srat_obj@meta.data)[grep("res_", colnames(srat_obj@meta.data))]

plot.dir <- file.path(dir, "plots")
if(!dir.exists(plot.dir)) dir.create(plot.dir)

cluster_plots.dir <- file.path(plot.dir, "cluster_plots")
if(!dir.exists(cluster_plots.dir)) dir.create(cluster_plots.dir)

# save normalized data in data slot (important)
srat_obj[["RNA"]]$data  <- as(object = srat_obj[["RNA"]]$counts, Class = "dgCMatrix")

for(res in resolutions){
  # UMAP with cluster 
  number_of_clusters = n_distinct(srat_obj@meta.data[res])
  
  # # UMAP split by clusters 
  dim_plt_split <- DimPlot(srat_obj, reduction = "umap", label = TRUE, cols = DiscretePalette_scCustomize(num_colors = number_of_clusters, palette = "varibow"),
                           group.by = res, split.by = res, ncol = 5) + ggtitle(str_interp("UMAP with Leiden Algorithm and ${res}"))
  file.name.split <- str_interp("plot_umap_split_leiden_${res}.pdf")
  ggsave(file.path(cluster_plots.dir, file.name.split), plot = dim_plt_split, width = 10, height = 10)
  
  # UMAP clustering (all)
  dim_plt <- DimPlot(srat_obj, reduction = "umap",
                     cols = DiscretePalette_scCustomize(num_colors = number_of_clusters, palette = "varibow"),
                     label = T, group.by = res) + ggtitle(str_interp("UMAP with Leiden Algorithm and ${res}"))
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
    labs(title = str_interp("Histogram of ${res}"),
         x = "Clusters",
         y = "Frequency") +
    scale_fill_identity() # keep seurat color code
  
  file.name.hist <- str_interp("plot_hist_${res}.pdf")
  ggsave(file.path(cluster_plots.dir, file.name.hist), width = 10, height = 10, bg= "white")
  
  
  # extract the average expression of the scaled.data
  combined_averages <- AverageExpression(srat_obj, return.seurat = TRUE, group.by = res)
  
  xxx= combined_averages@assays$RNA@layers$data
  rownames(xxx) = rownames(combined_averages)
  colnames(xxx) = colnames(combined_averages)
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
  
  
  # heatmap
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

