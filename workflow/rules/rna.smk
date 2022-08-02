rule seurat_build_reference:
    """
    Build Seurat reference dataset
    """
    input:
        mat = "reference/fetch/matrix.mtx",
        features = "reference/fetch/features.tsv",
        cells = "reference/fetch/barcodes.tsv",
        metadata = "reference/fetch/metadata.csv"
    output:
        project_out = "reference/seurat_build_reference/proj.rds",
        qc_violin = "reference/seurat_build_reference/qc_violin.pdf",
        qc_scatter = "reference/seurat_build_reference/qc_scatter.pdf",
        umap =  "reference/seurat_build_reference/umap.pdf"
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
        cells = "results/{sample}/fetch/barcodes.tsv"
    output:
        project_out = "results/{sample}/rna/seurat_build_rna/proj.rds",
        metadata = "results/{sample}/rna/seurat_build_rna/metadata.tsv",
        qc_violin = "results/{sample}/rna/seurat_build_rna/qc_violin.pdf",
        qc_scatter = "results/{sample}/rna/seurat_build_rna/qc_scatter.pdf",
    params:
        sample_name = lambda w: w.sample,
        seed = config["seurat_seed"]
    log:
        console = "logs/{sample}/rna/seurat_build_rna/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/build_seurat_rna.R"

rule seurat_doublets_rna:
    """
    Filter doublets
    """
    input:
        project_in = "results/{sample}/rna/seurat_build_rna/proj.rds"
    output:
        project_out_all = "results/{sample}/rna/seurat_doublets_rna/proj_all.rds",
        project_out_filtered = "results/{sample}/rna/seurat_doublets_rna/proj_filtered.rds"
    params:
        seed = config["seurat_seed"],
        doublet_rate = config["doublet_formation_rate"]
    log:
        console = "logs/{sample}/rna/seurat_doublets_rna/console.log"
    threads:
        config["max_threads_per_rule"]
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_doublets_rna.R"

rule seurat_label_rna:
    """
    Seurat RNA reference labeling
    """
    input:
        project_rna = "results/{sample}/rna/seurat_doublets_rna/proj_filtered.rds",
        project_ref = "reference/seurat_build_reference/proj.rds"
    output:
        project_out = "results/{sample}/rna/seurat_label_rna/proj.rds",
        umap = "results/{sample}/rna/seurat_label_rna/umap.pdf"
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/{sample}/rna/seurat_label_rna/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_label_rna.R"

rule seurat_merge_rna:
    """
    Merge RNA samples
    """
    input:
        projects_in = expand("results/{sample}/rna/seurat_label_rna/proj.rds", sample=samples_rna+samples_multiome)
    output:
        project_out = "results_merged/rna/seurat_merge_rna/proj.rds",
        umap_dataset = "results_merged/rna/seurat_merge_rna/umap_dataset.pdf",
        umap = "results_merged/rna/seurat_merge_rna/umap.pdf",
    params:
        seed = config["seurat_seed"],
        samples = samples_rna + samples_multiome,
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
        project_in = "results_merged/rna/seurat_merge_rna/proj.rds",
        # project_ref = "reference/seurat_build_reference/proj.rds"
    output:
        project_out = "results_merged/rna/seurat_cluster_rna/proj.rds",
        umap = "results_merged/rna/seurat_cluster_rna/umap_clusters.pdf",
    params:
        seed = config["seurat_seed"],
    log:
        console = "logs/merged/rna/seurat_cluster_rna/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_cluster_rna.R"

rule seurat_name_rna:
    """
    Seurat RNA cluster naming
    """
    input:
        project_in = "results_merged/rna/seurat_cluster_rna/proj.rds",
    output:
        project_out = "results_merged/rna/seurat_name_rna/proj.rds",
        umap = "results_merged/rna/seurat_name_rna/umap_clusters.pdf",
        umap_qc = "results_merged/rna/seurat_name_rna/umap_qc.pdf",
        metadata = "results_merged/rna/seurat_name_rna/metadata.tsv",
    params:
        seed = config["seurat_seed"],
        cluster_names = rna_cluster_names
    log:
        console = "logs/merged/rna/seurat_name_rna/console.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_name_rna.R"