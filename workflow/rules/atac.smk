
rule archr_build:
    """
    Build ArchR project
    """
    input:
        frag = "results/{sample}/fetch/fragments.tsv.gz",
        frag_ind = "results/{sample}/fetch/fragments.tsv.gz.tbi",
    output:
        qc_dir = directory("results/{sample}/atac/archr_qc"),
        project_dir = directory("results/{sample}/atac/archr_init"),
        arrow_dir = directory("results/{sample}/atac/archr_arrows")
    params:
        sample_name = lambda w: w.sample,
        seed = config["archr_seed"]
    log:
        console = "logs/{sample}/atac/archr_build/console.log",
        arrow_create = "logs/{sample}/atac/archr_build/arrow_create.log",
        doublets = "logs/{sample}/atac/archr_build/doublets.log",
        save = "logs/{sample}/atac/archr_build/save.log"
    threads:
        max_threads
    conda:
        "../envs/atac.yaml"
    script:
        "../scripts/build_archr_project.R"

rule archr_lsi:
    """
    ArchR dimensionality reduction
    """
    input:
        project_in = "results/{sample}/atac/archr_init"
    output:
        project_out = directory("results/{sample}/atac/archr_lsi")
    params:
        seed = config["archr_seed"]
    log:
        console = "logs/{sample}/atac/archr_lsi/console.log",
        move = "logs/{sample}/atac/archr_lsi/move.log",
        lsi_atac = "logs/{sample}/atac/archr_lsi/lsi_atac.log",
        save = "logs/{sample}/atac/archr_lsi/save.log"
    threads:
        max_threads
    conda:
        "../envs/atac.yaml"
    script:
        "../scripts/archr_lsi.R"

rule archr_cluster:
    """
    ArchR clustering
    """
    input:
        project_in = "results/{sample}/atac/archr_lsi"
    output:
        project_out = directory("results/{sample}/atac/archr_clustered")
    params:
        seed = config["archr_seed"]
    log:
        console = "logs/{sample}/atac/archr_atac/console.log",
        move = "logs/{sample}/atac/archr_atac/move.log",
        cluster_atac = "logs/{sample}/atac/archr_atac/cluster_atac.log",
        umap_plot = "logs/{sample}/atac/archr_atac/umap_plot.log",
        save = "logs/{sample}/atac/archr_atac/save.log"
    threads:
        max_threads
    conda:
        "../envs/atac.yaml"
    script:
        "../scripts/archr_cluster.R"


