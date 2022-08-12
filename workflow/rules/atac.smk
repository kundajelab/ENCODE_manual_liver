
rule archr_build:
    """
    Build ArchR project
    """
    input:
        frags = expand("results/{sample}/fetch/fragments.tsv.gz", sample=samples_atac+samples_multiome),
        frag_ind = expand("results/{sample}/fetch/fragments.tsv.gz.tbi", sample=samples_atac+samples_multiome),
    output:
        qc_dir = directory("results_merged/atac/archr_qc"),
        project_dir = directory("results_merged/atac/archr_init"),
        arrow_dir = directory("results_merged/atac/archr_arrows")
    params:
        sample_names = samples_atac + samples_multiome,
        seed = config["archr_seed"],
        min_frags = config["archr_min_frags"],
        min_tss_enr = config["archr_min_tss_enr"]
    log:
        console = "logs/merged/atac/archr_build/console.log",
        arrow_create = "logs/merged/atac/archr_build/arrow_create.log",
        doublets = "logs/merged/atac/archr_build/doublets.log",
        save = "logs/merged/atac/archr_build/save.log"
    threads:
        max_threads
    conda:
        "../envs/archr.yaml"
    script:
        "../scripts/build_archr_project.R"

rule archr_lsi:
    """
    ArchR dimensionality reduction
    """
    input:
        project_in = "results_merged/atac/archr_init"
    output:
        project_out = directory("results_merged/atac/archr_lsi")
    params:
        seed = config["archr_seed"]
    log:
        console = "logs/merged/atac/archr_lsi/console.log",
        move = "logs/merged/atac/archr_lsi/move.log",
        lsi_atac = "logs/merged/atac/archr_lsi/lsi_atac.log",
        save = "logs/merged/atac/archr_lsi/save.log"
    threads:
        max_threads
    conda:
        "../envs/archr.yaml"
    script:
        "../scripts/archr_lsi.R"

rule archr_harmony:
    """
    ArchR dimensionality reduction
    """
    input:
        project_in = "results_merged/atac/archr_lsi"
    output:
        project_out = directory("results_merged/atac/archr_harmony")
    params:
        seed = config["archr_seed"]
    log:
        console = "logs/merged/atac/archr_harmony/console.log",
        move = "logs/merged/atac/archr_harmony/move.log",
        umap_plot = "logs/merged/atac/archr_harmony/umap_plot.log",
        save = "logs/merged/atac/archr_harmony/save.log"
    threads:
        max_threads
    conda:
        "../envs/archr.yaml"
    script:
        "../scripts/archr_harmony.R"

rule archr_cluster:
    """
    ArchR clustering
    """
    input:
        project_in = "results_merged/atac/archr_harmony"
    output:
        project_out = directory("results_merged/atac/archr_clustered")
    params:
        seed = config["archr_seed"]
    log:
        console = "logs/merged/atac/archr_cluster/console.log",
        move = "logs/merged/atac/archr_cluster/move.log",
        cluster_atac = "logs/merged/atac/archr_cluster/cluster_atac.log",
        umap_plot = "logs/merged/atac/archr_cluster/umap_plot.log",
        save = "logs/merged/atac/archr_cluster/save.log"
    threads:
        max_threads
    conda:
        "../envs/archr.yaml"
    script:
        "../scripts/archr_cluster.R"

rule transport_rna_labels:
    """
    Parse labels from seurat metadata
    """
    input:
        seurat_data = "results_merged/rna/seurat_name_rna/metadata.tsv",
        wl_atac = "whitelists/atac.txt",
        wl_rna = "whitelists/rna.txt"
    output:
        "results_merged/atac/labels_import.tsv"
    conda:
        "../envs/fetch.yaml"
    script:
        "../scripts/transport_rna_clusters.py"

rule archr_label:
    """
    ArchR cluster labeling
    """
    input:
        project_in = "results_merged/atac/archr_clustered",
        label_data = "results_merged/atac/labels_import.tsv"
    output:
        project_out = directory("results_merged/atac/archr_label")
    params:
        seed = config["archr_seed"]
    log:
        console = "logs/merged/atac/archr_label/console.log",
        move = "logs/merged/atac/archr_label/move.log",
        umap_plot = "logs/merged/atac/archr_label/umap_plot.log",
        save = "logs/merged/atac/archr_label/save.log"
    threads:
        max_threads
    conda:
        "../envs/archr.yaml"
    script:
        "../scripts/archr_label.R"

rule archr_linkage:
    """
    ArchR unconstrained cluster linkage using reference
    """
    input:
        project_in = "results_merged/atac/archr_clustered",
        seurat_data = "reference/seurat_build_reference_log1p/proj.rds"
    output:
        project_out = directory("results_merged/atac/archr_linkage")
    params:
        seed = config["archr_seed"]
    log:
        console = "logs/merged/atac/archr_linkage/console.log",
        move = "logs/merged/atac/archr_linkage/move.log",
        umap_plot = "logs/merged/atac/archr_linkage/umap_plot.log",
        linkage = "logs/merged/atac/archr_linkage/linkage.log",
        save = "logs/merged/atac/archr_linkage/save.log"
    threads:
        max_threads
    conda:
        "../envs/archr.yaml"
    script:
        "../scripts/archr_linkage.R"

rule archr_write_qc:
    """
    ArchR write sample QC data
    """
    input:
        qc_dir = "results_merged/atac/archr_qc"
    output:
        out_dir = directory("results_merged/atac/archr_qc_parsed")
    log:
        console = "logs/merged/atac/archr_write_qc/console.log"
    conda:
        "../envs/archr.yaml"
    script:
        "../scripts/archr_write_qc.R"
