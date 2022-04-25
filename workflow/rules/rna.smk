rule seurat_build:
    """
    Build Seurat project
    """
    input:
        mat = "results/{sample}/fetch/mat.mtx",
        features = "results/{sample}/fetch/features.tsv",
        cells =  "results/{sample}/fetch/barcodes.tsv"
    output:
        project_out = "results/{sample}/rna/seurat_init/proj.rds",
        qc_violin = "results/{sample}/rna/seurat_init/qc_violin.pdf",
        qc_scatter = "results/{sample}/rna/seurat_init/qc_scatter.pdf",
        var_features = "results/{sample}/rna/seurat_init/var_features.pdf",
    params:
        sample_name = lambda w: w.sample,
        seed = config["seurat_seed"]
    log:
        console = "logs/{sample}/rna/seurat_init/console.log"
    conda:
        "../envs/rna.yaml"
    script:
        "../scripts/build_seurat_project.R"

rule seurat_cluster:
    """
    Seurat RNA clustering
    """
    input:
        project_in = "results/{sample}/rna/seurat_init/proj.rds"
    output:
        project_out = "results/{sample}/rna/seurat_clustered/proj.rds",
        umap = "results/{sample}/rna/seurat_clustered/umap.pdf"
    params:
        sample_name = lambda w: w.sample,
        seed = config["seurat_seed"],
    log:
        console = "logs/{sample}/rna/seurat_clustered/console.log"
    conda:
        "../envs/rna.yaml"
    script:
        "../scripts/seurat_cluster.R"