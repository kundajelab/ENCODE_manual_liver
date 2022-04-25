Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

params = snakemake@params 
input_paths = snakemake@input
output_paths = snakemake@output
log_paths = snakemake@log


proj <- readRDS(file = input_paths[["project_in"]])

proj <- RunPCA(proj, features = VariableFeatures(object = proj))

proj <- FindNeighbors(proj, dims = 1:10)
proj <- FindClusters(proj, resolution = 0.5)

proj <- RunUMAP(proj, dims = 1:10)

plt <- DimPlot(proj, reduction = "umap")
ggsave(output_paths[["umap"]], plt, device = "pdf")

saveRDS(proj, file = output_paths[["project_out"]])