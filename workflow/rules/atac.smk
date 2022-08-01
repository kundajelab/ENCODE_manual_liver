
rule archr_build:
    """
    Build ArchR project
    """
    input:
        frag = expand("results/{sample}/fetch/fragments.tsv.gz", sample=samples_atac+samples_multiome),
        frag_ind = expand("results/{sample}/fetch/fragments.tsv.gz.tbi", sample=samples_atac+samples_multiome),
    output:
        qc_dir = directory("results_merged/atac/archr_qc"),
        project_dir = directory("results_merged/atac/archr_init"),
        arrow_dir = directory("results_merged/atac/archr_arrows")
    params:
        sample_names = samples_atac + samples_multiome,
        seed = config["archr_seed"],
        min_frags = config["archr_min_frags"]
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

# rule archr_peakmatrix:
#     """
#     ArchR clustering
#     """
#     input:
#         project_in = "results/{sample}/atac/archr_clustered"
#     output:
#         project_out = directory("results/{sample}/atac/archr_peakmatrix"),
#         mat_out = "results/{sample}/atac/archr_peakmatrix_mat/mat.mtx",
#         barcodes_out = "results/{sample}/atac/archr_peakmatrix_mat/barcode_metadata.tsv",
#         peaks_out = "results/{sample}/atac/archr_peakmatrix_mat/peaks.tsv",
#     params:
#         seed = config["archr_seed"]
#     log:
#         console = "logs/{sample}/atac/archr_peakmatrix/console.log",
#         move = "logs/{sample}/atac/archr_peakmatrix/move.log",
#         pseudobulks = "logs/{sample}/atac/archr_peakmatrix/pseudobulks.log",
#         call_peaks = "logs/{sample}/atac/archr_peakmatrix/call_peaks.log",
#         add_peak_mat = "logs/{sample}/atac/archr_peakmatrix/add_peak_mat.log",
#         save_peak_mat = "logs/{sample}/atac/archr_peakmatrix/save_peak_mat.log",
#         save = "logs/{sample}/atac/archr_peakmatrix/save.log"
#     threads:
#         max_threads
#     conda:
#         "../envs/archr.yaml"
#     script:
#         "../scripts/archr_peakmatrix.R"
