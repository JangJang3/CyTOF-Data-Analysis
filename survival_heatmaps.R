################################################################
# Set working environment 
################################################################
current_path <- "/mnt/volumeb/jin/baseline_pipeline/analysis"

getwd()
setwd(dir=current_path)

################################################################
# Installing
################################################################

source("/mnt/volumeb/jin/installations.R")


################################################################
# Create File Path and Folders
################################################################
print("Status:  Create File Paths")


dir_name <- "survival_heatmaps"

# create data folder where all analysis will be stored
if(!dir.exists(dir_name)) dir.create(dir_name)

# set it as a main directory 
dir <- file.path(getwd(), dir_name)



################################################################
# Create Survival Heatmap - Mega CLuster 
################################################################

srat_obj <- readRDS("seurat_object_leiden_merged_clusters.rds")


# add annotations 
# Cluster 1 - monocytes
# Cluster 2 - CD4+ T cells
# Cluster 3 - CD8+ T cells
# Cluster 4 - NK cells
# Cluster 5 - B cells
# CLuster 6 - V delta 2 cells(, but this one is weird because it contains MAIT cells as well)
# Cluster 7 - V delta 1

Idents(srat_obj) <- "new_res_0.1"

cluster_ids.new <- c(
  "1" = "Monocytes",
  "2" = "CD4+ T cells",
  "3" = "CD8+ T cells",
  "4" = "NK cells",
  "5" = "B cells",
  "6" = "Vδ2",
  "7" = "Vδ1"
)
names(cluster_ids.new) <- levels(srat_obj)
srat_obj <- RenameIdents(srat_obj, cluster_ids.new)
srat_obj$annotation <- Idents(srat_obj)


#' Creates a list of matrices that contains the overall cluster frequency for each HCC patient
#' Differentiated by timepoint (before and after treatment)
#' 
#' @param seurat_obj A seurat object 
#' @return A list with two elements that consists of matrices (before and after treatment)
create_heatmap_matrix <- function(seurat_obj){
  # create the overall table (where baseline and after treatment are not separated)
  df <- seurat_obj@meta.data %>%
    filter(!str_starts(sample_id, "HD_")) %>%
    group_by(patient_id, timepoint, annotation) %>%
    summarise(cell_count = n()) %>%
    group_by(patient_id, timepoint) %>%
    mutate(frequency_prop = cell_count / sum(cell_count)) %>%
    ungroup()
  
  # correct the naming convention (if needed)
  df$patient_id <- gsub("PatID", "Pat_ID", df$patient_id)
  
  
  heatmap_data_T0 <- df %>%
    filter(timepoint == "T0") %>%
    mutate(patient_time = patient_id) %>%   # or keep as patient_id only
    select(patient_time, annotation, frequency_prop) %>%
    pivot_wider(names_from = annotation, values_from = frequency_prop, values_fill = 0) %>%
    column_to_rownames(var = "patient_time")
  
  heatmap_matrix_T0 <- as.matrix(t(heatmap_data_T0))
  
  ############
  
  heatmap_data_T1 <- df %>%
    filter(timepoint == "T1") %>%
    mutate(patient_time = patient_id) %>%
    select(patient_time, annotation, frequency_prop) %>%
    pivot_wider(names_from = annotation, values_from = frequency_prop, values_fill = 0) %>%
    column_to_rownames(var = "patient_time")
  
  heatmap_matrix_T1 <- as.matrix(t(heatmap_data_T1))
  
  return(list(heatmap_matrix_T0, heatmap_matrix_T1))
}


heatmap_matrix_mega_cluster <- create_heatmap_matrix(srat_obj)
heatmap_matrix_T0 <- heatmap_matrix_mega_cluster[1]
heatmap_matrix_T1 <- heatmap_matrix_mega_cluster[2]



# make the annotation column with responder and non- responder, progression-free date or overall survival
clinical_data <- readRDS("/mnt/volumeb/jin/survival_analysis/new_clinical_data_filtered_July.rds")




# each square symbolizes the cluster frequency overall 

# |-|-|-|-|-|-|-|-|-|-|-| OS or PFS (in Days)
# |-|-|-|-|-|-|-|-|-|-|-| Responder vs Non-Responder
# ┌─────┬─────┬─────┬─────┬─────┬─────┐
# Cluster A│  30 │  40 │ ... │     │     │     │ 
# Cluster B│  20 │  10 │     │     │     │     │
# Cluster C│  50 │  50 │     │     │     │     │
# ─────────┼─────┼─────┼─────┼─────┼─────┼─────┤
# Cluster A│  80 │  60 │     │     │     │     │ 
# Cluster B│  10 │  30 │     │     │     │     │
# Cluster C│  10 │  10 │     │     │     │     │
# └─────┴─────┴─────┴─────┴─────┴─────┘
# P1    P2   P3   ...           P10

#' Creates a survival heatmap plot (see outline, above) 
#' 
#' @param matrix A matrix object that contains the overall cluster frequency 
#' @param clinical_data A file with the clinical information from the HCC patients 
#' @param file_name Output file nameing 
#' @param timepoint Can be either T0 or T1 (baseline vs after treatment)
#' @param survival_time Can be either PFS or OS 
#' @param hm_title Title for the heatmap
#' @param width
#' @param height
#' 
#' @return Not an object, but saves the survival heatmap in the directory 
create_heatmap <- function(matrix, 
                           clinical_data, 
                           file_name = "heatmap_mega_cluster_baseline.pdf", 
                           timepoint = "T0", 
                           survival_time = "PFS", 
                           hm_title = "Mega Cluster",
                           width = 10,
                           height = 3){
  
  days <- if (survival_time == "PFS") "PFS_11_24" else "OS_11_24"
  
  annotation_col <- clinical_data %>%
    select(patient_id, !!sym(days), survival_status)
  rownames(annotation_col) <- annotation_col$patient_id
  
  pfs_col_fun <- colorRamp2(
    breaks = seq(0, max(clinical_data[[days]], na.rm = TRUE), length.out = 3),
    colors = c("#2c7bb6", "#abd9e9", "orange")  
  )
  
  # os_col_fun <- colorRamp2(
  #   c(0, 750, 1500),
  #   c("white", "#f4a6c6", "#8b0000")
  # )
  
  freq_col_fun <- colorRamp2(
    breaks = seq(0, max(matrix), length.out = 3),
    #colors = viridis::plasma(3)
    colors = c("white", "blue", "navy")
  )
  
  # staging_colors <- c(
  #   "Response" = "#73e67a",
  #   "Progress" = "#ff3333"
  # )
  
  surviva_col_fun <- c(
    "Short term" = "#ff3333",
    "Mid term" = "yellow",
    "Long term" = "#73e67a"
  )
  
  annotation_col_ordered <- annotation_col[order(-annotation_col[[days]], na.last = TRUE), ]
  heatmap_matrix_ordered <- matrix[, rownames(annotation_col_ordered)]
  
  
  ha_top <- HeatmapAnnotation(
    "{survival_time}" = annotation_col_ordered[[days]],
    "Survival status" = annotation_col_ordered$survival_status,
    gp = gpar(col = "black", lwd = 0.5),
    col = list(
      PFS = pfs_col_fun,
      #OS_11_24 = os_col_fun,
      "Survival status" = surviva_col_fun
    ),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 10)
  )
  
  # check between T0 and T1 
  if(timepoint == "T0"){
    timepoint_name = "at Baseline (T0)"
  } else {
    timepoint_name = "after Treatment (T1)"
  }
  
  hm <- Heatmap(
    heatmap_matrix_ordered,
    name = "Cluster Frequency",
    column_title = str_interp("Cluster Frequencies ${timepoint_name} - ${hm_title}"),
    column_title_gp = gpar(fontsize = 15),
    col = freq_col_fun,
    top_annotation = ha_top,
    rect_gp = gpar(col = "black", lwd = 0.5),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 10),  # Font size for x-axis (patient labels)
    row_names_gp = gpar(fontsize = 10),     # Font size for y-axis (cluster labels)
    show_row_names = TRUE
  )
  
  cairo_pdf(file.path(dir, file_name), width = width, height = height)
  draw(hm, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  
}

create_heatmap(heatmap_matrix_T0, clinical_data)
create_heatmap(heatmap_matrix_T1, 
               clinical_data, 
               file_name = "heatmap_mega_cluster_after_treatment.pdf",
               timepoint = "T1",
               hm_title = "Mega Cluster")


################################################################
# Create Heatmap - NK Cluster 
################################################################

nk_cells <- readRDS("/mnt/volumeb/jin/baseline_pipeline/analysis/analysis_only_NK_cells_with_CD3/Analysis_remove_CD3_cluster_6/Analysis_remove_CD3_cluster_6_rerun_0.2/annotation_abundances/seurat_object_annotated.rds")
#nk_cells$annotation <- Idents(nk_cells)


heatmap_matrix_nk_cells <- create_heatmap_matrix(nk_cells)
heatmap_matrix_T0 <- heatmap_matrix_nk_cells[[1]]
heatmap_matrix_T1 <- heatmap_matrix_nk_cells[[2]]

create_heatmap(heatmap_matrix_T0,
               clinical_data,
               file_name = "heatmap_nk_cells_baseline.pdf",
               hm_title = "NK cells",
               width = 10,
               height = 3) 

create_heatmap(heatmap_matrix_T1,
               clinical_data,
               file_name = "heatmap_nk_cells_after_treatment.pdf",
               timepoint = "T1",
               hm_title = "NK cells",
               width = 10,
               height = 3) 


################################################################
# Create Heatmap - T cell cluster
################################################################

t_cells <- readRDS("/mnt/volumeb/jin/baseline_pipeline/analysis/analysis_only_Tcells_test_run/Merged_clusters/annotation_abundances/seurat_object_annotated.rds")
#t_cells$annotation <- Idents(t_cells)


heatmap_matrix_t_cells <- create_heatmap_matrix(t_cells)
heatmap_matrix_T0 <- heatmap_matrix_t_cells[[1]]
heatmap_matrix_T1 <- heatmap_matrix_t_cells[[2]]

create_heatmap(heatmap_matrix_T0,
               clinical_data,
               file_name = "heatmap_t_cells_baseline.pdf",
               hm_title = "T cells",
               width = 10,
               height = 4) 

create_heatmap(heatmap_matrix_T1,
               clinical_data,
               file_name = "heatmap_t_cells_after_treatment.pdf",
               timepoint = "T1",
               hm_title = "T cells",
               width = 10,
               height = 4) 


################################################################
# Create Heatmap - CD8+ T cell cluster
################################################################

cd8_t_cells <- readRDS("/mnt/volumeb/jin/baseline_pipeline/analysis/analysis_only_Tcells_test_run/Analysis_CD8+_T_cells/Merged_clusters/Merged_clusters/annotation_abundances/seurat_object_annotated.rds")
#t_cells$annotation <- Idents(t_cells)


heatmap_matrix_cd8_t_cells <- create_heatmap_matrix(cd8_t_cells)
heatmap_matrix_T0 <- heatmap_matrix_cd8_t_cells[[1]]
heatmap_matrix_T1 <- heatmap_matrix_cd8_t_cells[[2]]

create_heatmap(heatmap_matrix_T0,
               clinical_data,
               file_name = "heatmap_cd8_t_cells_baseline.pdf",
               hm_title = "CD8+ T cells",
               width = 10,
               height = 3) 

create_heatmap(heatmap_matrix_T1,
               clinical_data,
               file_name = "heatmap_cd8_t_cells_after_treatment.pdf",
               timepoint = "T1",
               hm_title = "T cells",
               width = 10,
               height = 3) 

