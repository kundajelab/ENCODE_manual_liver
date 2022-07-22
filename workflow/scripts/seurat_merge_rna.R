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

set.seed(params[["seed"]])


projs <- lapply(input_paths[["projects_in"]], readRDS)

proj_merged <- merge(projs[[1]], projs[-1], project = "merged_rna")

proj_merged <- NormalizeData(proj_merged)
proj_merged <- FindVariableFeatures(proj_merged, selection.method = "vst", nfeatures = 2000)

saveRDS(proj_merged, file = output_paths[["project_out"]])

sink(type = "message")
sink()