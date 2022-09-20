# ENCODE Manual Analyses: Liver ATAC Datasets

Analyses for integrating ENCODE single-nucleus ATAC-Seq liver datasets.

Contact information: 
Austin Wang 
atwang@stanford.edu

## Running the analyses

### Requirements

- A Linux-based OS
- A conda-based Python 3 installation
- [Snakemake v6.6.1+](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) (full installation)
- An ENCODE DCC account with access to the necessary datasets

### Execution

1. Install any necessary requirements above
2. Download the pipeline
    ```
    git clone https://github.com/kundajelab/ENCODE_scatac
    ```
3. Activate the `snakemake` conda environment:
    ```
    conda activate snakemake
    ```
4. Run the pipeline:
    ```
    snakemake -k --use-conda --cores $NCORES 
    ```
    Here, `$NCORES` is the number of cores to utilize

When run for the first time, the pipeline will take some time to install conda packages. 

## Analysis Steps

- Cell filtering based on minimum fragment count and TSS enrichment
- Iterative LSI dimensionality reduction (ArchR)
- Louvain Clustering
- Dataset integration with Harmony
- Cell type labeling guided by [Guilliams et al. *Cell* 2022](https://www.cell.com/cell/fulltext/S0092-8674(21)01481-1), in addition to the ENCODE Liver RNA datasets

All parameters can be found in `config/config.yaml`.

## Data Inputs and Outputs

### Inputs

The pipeline will automatically pull the required datasets from the ENCODE portal. The relevant ENCODE ID's can be found in `config/samples_atac.tsv` and `config/samples_multiome.tsv`.

The pipeline will additionally download relevant non-ENCODE data files. URLs to these can be found in `config/config.yaml`. 

### Outputs

This pipeline contains code for both the RNA and ATAC liver data analyses. Running the pipeline will generate results for both the RNA and ATAC analyses. ATAC outputs are organized in the `export/atac` directory, relative to the current working directory.

```
atac
├── metadata.tsv.gz # Cell metadata file
├── embeddings # Embedding coordinates files
│   └── $EMBEDDING.tsv.gz
├── labels # Cell type labels files 
│   └── $CELL_LABELS_SET.tsv.gz
├── markers # Directory for cell type object marker gene data
│   └── $LABEL_TEMP_ID.tsv.gz # Cell type marker gene file 
├── markers_aux # Directory for cell type object auxiliary data
│   └── $LABEL_TEMP_ID.tar.gz # Tarball for cell type auxiliary data 
├── figures.tar.gz # Tarball for figure data
├── auxiliary_data.tar.gz # Tarball for auxiliary data
└── datasets.txt # A list of datasets used (ENCODE IDs) 
```

Note: intermediate pipeline outputs will also be placed in the current working directory.