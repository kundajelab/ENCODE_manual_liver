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

proj <- readRDS(file = input_paths[["project_rna"]])

ref <- readRDS(file = input_paths[["project_ref"]])
# print(head(ref@meta.data)) ####
# print(ref) ####
# ref <- subset(x = ref, subset = `Factor.Value.inferred.cell.type...authors.labels.` != "")
# print(ref) ####
# print(proj) ####

anchors <- FindTransferAnchors(
  reference = ref,
  query = proj,
  normalization.method = "SCT",
  dims = 1:30, 
#   reference.reduction = "pca",
  reduction = "cca"
)
proj <- MapQuery(
  anchorset = anchors,
  query = proj,
  reference = ref,
  refdata = "cell_type",
  reference.reduction = "pca"
)
proj$cell_type <- proj$predicted.id

proj <- FindNeighbors(proj, dims = 1:30, reduction = "harmony")
proj <- RunUMAP(proj, dims = 1:30, reduction = "harmony")

plt <- DimPlot(proj, reduction = "umap", group.by = "cell_type")
ggsave(output_paths[["umap"]], plt, device = "pdf")

saveRDS(proj, file = output_paths[["project_out"]])