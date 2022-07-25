rule seurat_build_reference:
    """
    Build Seurat reference dataset
    """
    input:
        mat = "reference/fetch/matrix.mtx",
        features = "reference/fetch/features.tsv",
        cells = "reference/fetch/barcodes.tsv",
        metadata = "reference/fetch/metadata.tsv"
    output:
        project_out = "reference/seurat_build_reference/proj.rds",
        qc_violin = "reference/seurat_build_reference/qc_violin.pdf",
        qc_scatter = "reference/seurat_build_reference/qc_scatter.pdf",
        var_features = "reference/seurat_build_reference/var_features.pdf",
    params:
        seed = config["seurat_seed"]
    log:
        console = "logs/reference/seurat_build_reference/seurat_build_rna/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/build_seurat_reference.R"

rule seurat_build_rna:
    """
    Build Seurat project
    """
    input:
        mat = "results/{sample}/fetch/matrix.mtx",
        features = "results/{sample}/fetch/features.tsv",
        cells =  "results/{sample}/fetch/barcodes.tsv"
    output:
        project_out = "results/{sample}/rna/seurat_build_rna/proj.rds",
        qc_violin = "results/{sample}/rna/seurat_build_rna/qc_violin.pdf",
        qc_scatter = "results/{sample}/rna/seurat_build_rna/qc_scatter.pdf",
        var_features = "results/{sample}/rna/seurat_build_rna/var_features.pdf",
    params:
        sample_name = lambda w: w.sample,
        seed = config["seurat_seed"]
    log:
        console = "logs/{sample}/rna/seurat_build_rna/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/build_seurat_rna.R"

rule seurat_merge_rna:
    """
    Merge RNA samples
    """
    input:
        projects_in = expand("results/{sample}/rna/seurat_build_rna/proj.rds", sample=samples_rna)
    output:
        project_out = "results_merged/rna/seurat_merge_rna/proj.rds",
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/merged/rna/seurat_merge_rna/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_merge_rna.R"

rule seurat_cluster_rna:
    """
    Seurat RNA clustering
    """
    input:
        project_rna = "results_merged/rna/seurat_merge_rna/proj.rds",
        project_ref = "reference/seurat_build_reference/proj.rds"
    output:
        project_out = "results_merged/rna/seurat_cluster_rna/proj.rds",
        umap = "results_merged/rna/seurat_cluster_rna/umap.pdf",
        umap_ref = "results_merged/rna/seurat_cluster_rna/umap_ref.pdf"
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/merged/rna/seurat_cluster_rna/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_cluster_rna.R"