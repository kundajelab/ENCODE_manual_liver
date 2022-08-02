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

proj <- readRDS(file = input_paths[["project_in"]])

new.cluster.ids <- unlist(params[["cluster_names"]])
names(new.cluster.ids) <- levels(proj)
proj <- RenameIdents(proj, new.cluster.ids)

plt <- DimPlot(proj, reduction = "umap", label = TRUE)
ggsave(output_paths[["umap"]], plt)

proj$log10_count_RNA <- log10(proj$nCount_RNA)

plt <- FeaturePlot(object = proj, features = "log10_count_RNA")
ggsave(output_paths[["umap_qc"]], plt)

write.table(proj@meta.data, file=output_paths[["metadata"]], quote=FALSE, sep='\t', col.names = NA)

saveRDS(proj, file = output_paths[["project_out"]])