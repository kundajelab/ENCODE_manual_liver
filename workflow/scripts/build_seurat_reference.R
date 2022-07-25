Sys.setenv(CONDA_BUILD_SYSROOT = "/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

params <- snakemake@params
input_paths <- snakemake@input
output_paths <- snakemake@output
log_paths <- snakemake@log

set.seed(params[["seed"]])

# Load the proj dataset
expression_matrix <- ReadMtx(
  mtx = input_paths[["mat"]], 
  features = input_paths[["features"]],
  cells = input_paths[["cells"]]
)
# print(expression_matrix) ####
metadata <- read.table(file = input_paths[["metadata"]], sep = '\t', header = TRUE)
rownames(metadata) <- metadata$Assay
print(head(metadata)) ####

# Initialize the Seurat object with the raw (non-normalized data).
proj <- CreateSeuratObject(
    counts = expression_matrix, 
    project = "reference", 
    min.cells = 3, 
    min.features = 200, 
    meta.data = metadata
)
proj[["percent.mt"]] <- PercentageFeatureSet(proj, pattern = "^MT-")
print(head(proj@meta.data)) ####

plt <- VlnPlot(proj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(output_paths[["qc_violin"]], plt, device = "pdf")

plot1 <- FeatureScatter(proj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(proj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plt <- plot1 + plot2
ggsave(output_paths[["qc_scatter"]], plt, device = "pdf")

proj <- subset(proj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

proj <- NormalizeData(proj, normalization.method = "LogNormalize", scale.factor = 10000)

proj <- FindVariableFeatures(proj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(proj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(proj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plt <- plot1 + plot2
ggsave(output_paths[["var_features"]], plt, device = "pdf")

all.genes <- rownames(proj)
proj <- ScaleData(proj, features = all.genes)

proj <- RunPCA(proj, features = VariableFeatures(object = proj))

proj <- FindNeighbors(proj, dims = 1:10)

proj <- RunUMAP(proj, dims = 1:10, return.model = TRUE)

plt <- DimPlot(proj, reduction = "umap", group.by = "Factor.Value.cell.type.")
ggsave(output_paths[["umap"]], plt, device = "pdf")

saveRDS(proj, file = output_paths[["project_out"]])

sink(type = "message")
sink()