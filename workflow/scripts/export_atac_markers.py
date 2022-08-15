import os
import mygene as mg

HEADER = """
# Marker genes for {label_name}
gene_id: ENSEMBL Gene ID
gene_name: Common name of the gene
avg_log2FC: Log2 fold-chage of the average expression between the given cell type and the remaining cells (min 0.5)
FDR: p-value-derived FDR value for gene 
is_enriched: A binary (0/1) value indicating whether the given gene is enriched for the cell type (FDR < 0.1)
"""

COLUMNS = "gene_id\tgene_name\tavg_log2FC\tavg_log2FC\tFDR\tis_enriched\n"

def query_gene(gene_name):
    query, = mg.querymany([gene_name], scopes='symbol', fields='ensembl.gene', species='human')
    if 'ensembl.gene' not in query:
        return gene_name
    if isinstance(query['ensembl.gene'], str): 
        return query['ensembl.gene']

    return query['ensembl.gene'][0]


def export_label(in_path, out_path, gene_cache, label_name):
    with open(in_path) as f, open(out_path, "w") as fo:
        fo.write(HEADER.format(label_name=label_name))
        fo.write(COLUMNS)

        h = f.readline().rstrip('\n').split('\t')
        gene_name_ind = h.index("name")
        lfc_ind = h.index("Log2FC")
        fdr_ind = h.index("FDR")

        for line in f:
            entries = line.rstrip('\n').split('\t')

            gene_name = entries[gene_name_ind]
            avg_logFC = float(entries[lfc_ind])
            fdr = float(entries[fdr_ind])

            gene_id = gene_cache.setdefault(gene_name, query_gene(gene_name))

            is_enriched = int((fdr < 0.1) and (avg_logFC > 0))
            
            line = f"{gene_id}\t{gene_name}\t{avg_logFC}\t{fdr}\t{is_enriched}\n"
            fo.write(line)

def main(markers_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    labels = [i for i in os.listdir(markers_dir) if not i.startswith(".")]
    gene_cache = {}
    for l in labels:
        markers_path = os.path.join(markers_dir, l)
        out_path = os.path.join(out_dir, l)
        label_name = l.split(".")[0].replace("_", " ")
        export_label(markers_path, out_path, gene_cache, label_name)

markers_dir, = snakemake.input

out_dir, = snakemake.output

main(markers_dir, out_dir)