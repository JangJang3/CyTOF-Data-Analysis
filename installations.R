################################################################
# Install and Load libraries 
################################################################

# (needed maybe for installation, use in terminal if missing)
# sudo apt install cmake
# sudo apt install libcairo2-dev
# sudo apt install libfontconfig1-dev
# sudo apt install libssl-dev
# sudo apt install libxml2-dev
# sudo apt install libharfbuzz-dev libfribidi-dev 
# sudo apt install libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev 
# sudo apt install libcurl4-openssl-dev
# sudo apt-get install libhdf5-dev 


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if(!require("devtools")){
  install.packages("devtools")
  library(devtools)
}

if(!require("flowCore")){
  BiocManager::install("flowCore")
  library(flowCore)
}

if(!require("SummarizedExperiment")){
  BiocManager::install("SummarizedExperiment")
  library(SummarizedExperiment)
}

if(!require("SingleCellExperiment")){
  BiocManager::install("SingleCellExperiment")
  library(SingleCellExperiment)
}

if(!require("openCyto")){
  BiocManager::install("openCyto")
  library(openCyto)
}

if (!require("CytoML")) {
  BiocManager::install("CytoML")
  library(CytoML)
}

if(!require("flowDensity")){
  BiocManager::install("flowDensity")
  library(flowDensity)
}

if(!require("FlowSOM")){
  BiocManager::install("FlowSOM")
  library(FlowSOM)
}

if(!require("CATALYST")){
  BiocManager::install("CATALYST")
  library(CATALYST)
}

if(!require("CytoNorm")){
  remotes::install_github("saeyslab/CytoNorm")
  library(CytoNorm)
}


if(!require("cytutils")){
  install_github("ismmshimc/cytutils")
  library(cytutils)
}


if(!require(cytofclean)){
  devtools::install_github("JimboMahoney/cytofclean", 
                           dependencies = TRUE)
  library("cytofclean")
}

if(!require("BPCells")){
  install.packages('BPCells', repos = c('https://bnprks.r-universe.dev', 'https://cloud.r-project.org'))
  library(BPCells)
}

packages_needed <- c("Seurat", "ggplot2","tcltk2", "dittoSeq", "ComplexHeatmap", "ggalluvial",
  "clustree", "dplyr", "DT", "openxlsx", "tidyverse", "future", "viridis", "scico", "circlize",
  "future.apply", "patchwork", "cowplot", "glue", "gridExtra", "scCustomize",
  "RColorBrewer", "remotes", "CytoNorm", "flowAI",
  "flowCut", "stringr", "cytutils", "pheatmap", "uwot",
  "ggpubr", "flowWorkspace", "pals", "scales"
)

for (pkg in packages_needed) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

