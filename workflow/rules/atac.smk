
rule archr_build:
    """
    Build ArchR project
    """
    input:
        frag = "results/{sample}/fetch/fragments.tsv.gz",
        frag_ind = "results/{sample}/fetch/fragments.tsv.gz.tbi",
    output:
        qc_dir = directory("results/{sample}/cluster/archr_qc"),
        project_dir = directory("results/{sample}/cluster/archr_init"),
        arrow_dir = directory("results/{sample}/cluster/archr_arrows")
    params:
        sample_name = lambda w: w.sample,
        seed = config["archr_seed"]
    log:
        console = "logs/{sample}/cluster/archr_build/console.log",
        arrow_create = "logs/{sample}/cluster/archr_build/arrow_create.log",
        doublets = "logs/{sample}/cluster/archr_build/doublets.log",
        save = "logs/{sample}/cluster/archr_build/save.log"
    threads:
        max_threads
    conda:
        "../envs/cluster.yaml"
    script:
        "../scripts/build_archr_project.R"

rule archr_lsi:
    """
    ArchR dimensionality reduction
    """
    input:
        project_in = "results/{sample}/cluster/archr_init"
    output:
        project_out = directory("results/{sample}/cluster/archr_lsi")
    params:
        seed = config["archr_seed"]
    log:
        console = "logs/{sample}/cluster/archr_lsi/console.log",
        move = "logs/{sample}/cluster/archr_lsi/move.log",
        lsi_atac = "logs/{sample}/cluster/archr_lsi/lsi_atac.log",
        save = "logs/{sample}/cluster/archr_lsi/save.log"
    threads:
        max_threads
    conda:
        "../envs/cluster.yaml"
    script:
        "../scripts/archr_lsi.R"

rule archr_cluster:
    """
    ArchR clustering
    """
    input:
        project_in = "results/{sample}/cluster/archr_lsi"
    output:
        project_out = directory("results/{sample}/cluster/archr_clustered")
    params:
        seed = config["archr_seed"]
    log:
        console = "logs/{sample}/cluster/archr_cluster/console.log",
        move = "logs/{sample}/cluster/archr_cluster/move.log",
        cluster_atac = "logs/{sample}/cluster/archr_cluster/cluster_atac.log",
        umap_plot = "logs/{sample}/cluster/archr_cluster/umap_plot.log",
        save = "logs/{sample}/cluster/archr_cluster/save.log"
    threads:
        max_threads
    conda:
        "../envs/cluster.yaml"
    script:
        "../scripts/archr_cluster.R"


