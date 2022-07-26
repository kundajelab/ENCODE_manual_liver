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
print(head(ref@meta.data)) ####
ref <- subset(x = ref, subset = `Factor.Value.inferred.cell.type...authors.labels.` != NA)

anchors <- FindTransferAnchors(
  reference = ref,
  query = proj,
  dims = 1:30, 
  reference.reduction = "pca"
)
proj <- MapQuery(
  anchorset = anchors,
  query = proj,
  reference = ref,
  refdata = "Factor.Value.inferred.cell.type...authors.labels.",
  reference.reduction = "pca"
  # reduction.model = "umap"
)
print(head(proj@meta.data)) ####
print(head(ref@meta.data)) ####
proj$cell_type <- proj$predicted.id

# plt_ref <- DimPlot(proj, reduction = "umap.ref", group.by = "cell_type")
# ggsave(output_paths[["umap_ref"]], plt, device = "pdf")

all.genes <- rownames(proj)
proj <- ScaleData(proj, features = all.genes)

proj <- RunPCA(proj, features = VariableFeatures(object = proj))

proj <- FindNeighbors(proj, dims = 1:10)

proj <- RunUMAP(proj, dims = 1:10)

plt <- DimPlot(proj, reduction = "umap", group.by = "cell_type")
ggsave(output_paths[["umap"]], plt, device = "pdf")

saveRDS(proj, file = output_paths[["project_out"]])