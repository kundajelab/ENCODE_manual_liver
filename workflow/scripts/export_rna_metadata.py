import gzip
import os

HEADER = """
# cell_id: Cell ID used in integrated analysis
# rna_dataset: ENCODE snRNA-Seq dataset ID
# rna_barcode: snRNA-Seq barcode
# atac_dataset: ENCODE snATAC-Seq dataset ID
# atac_barcode: snATAC-Seq barcode
# rna_umi_count: snRNA UMI's per cell
# atac_fragment_count: snATAC-Seq fragments per cell
# rna_frac_mito: Fraction of mitochondrial RNA reads
# rna_frac_ribo: Fraction of ribosomal RNA reads
# atac_tss_enrichment: snATAC-Seq transcription start site enrichment
# passed_filtering: A binary (0/1) value indicating whether the cell passed manual filtering
"""

COLUMNS = (
    "cell_id\trna_dataset\trna_barcode\tatac_dataset\tatac_barcode\t"
    "rna_umi_count\tatac_fragment_count\trna_frac_mito\trna_frac_ribo\tatac_tss_enrichment\tpassed_filtering\n"
)

REV_COMP = str.maketrans("ATGC", "TACG")
def reverse_complement(seq):
    return str.translate(seq, REV_COMP)[::-1]

def load_whitelists(wl_atac_path, wl_rna_path):
    bc_map = {}
    with open(wl_atac_path) as a, open(wl_rna_path) as r:
        for line_a, line_r in zip(a, r):
            bc_a = line_a.rstrip("\n")
            bc_r = line_r.rstrip("\n")

            bc_map[bc_r] = bc_a

    return bc_map

def load_atac_bcs(metadata_atac_path):
    atac_bcs = {}
    with open(metadata_atac_path) as f:
        h = f.readline().rstrip('\n').split('\t')
        cell_id_ind = h.index("cellNames")

        for line in f:
            entries = line.rstrip('\n').split('\t')

            cell_id = entries[cell_id_ind]
            dataset, barcode = cell_id.split("#")
            dataset_parsed_atac = dataset.split("_")[-1].split("-")[0]
            atac_bcs.setdefault(dataset_parsed_atac, set()).add(barcode)

    return atac_bcs

def load_records(metadata_path, atac_bcs, bc_map):
    records = {}
    with open(metadata_path) as f:
        h = f.readline().rstrip('\n').split('\t')
        barcode_ind = 0
        dataset_ind = h.index("orig.ident")
        count_ind = h.index("nCount_RNA")
        mito_ind = h.index("percent.mt")
        ribo_ind = h.index("percent.ribo")

        for line in f:
            entries = line.rstrip('\n').split('\t')

            dataset = entries[dataset_ind]
            barcode = entries[barcode_ind]
            cell_id = f"{dataset}_{barcode}"
            dataset_parsed = dataset.split("_")[0].split("-")[0]

            dataset_parsed_atac = dataset.split("_")[-1].split("-")[0]
            if dataset_parsed_atac in atac_bcs:
                barcode_trans = bc_map[barcode]
                barcode_trans_rc = reverse_complement(bc_map[barcode])
                if barcode_trans in atac_bcs[dataset_parsed_atac]:
                    barcode_atac = barcode_trans
                elif barcode_trans_rc in atac_bcs[dataset_parsed_atac]:
                    barcode_atac = barcode_trans_rc
                else:
                    barcode_atac = "NA"
            else:
                dataset_parsed_atac = "NA"
                barcode_atac = "NA"
            
            umi_count = int(entries[count_ind])
            frac_mito = float(entries[mito_ind]) / 100
            frac_ribo = float(entries[ribo_ind]) / 100

            record = {
                "dataset": dataset_parsed,
                "dataset_atac": dataset_parsed_atac,
                "barcode": barcode,
                "barcode_atac": barcode_atac,
                "umi_count": umi_count,
                "frac_mito": frac_mito,
                "frac_ribo": frac_ribo,
                "pass_filter": False
            }
            records[cell_id] = record

    return records

def load_final_data(final_data_path):
    ids = []
    with open(final_data_path) as f:
        next(f)
        for line in f:
            cell_id = line.rstrip('\n').split('\t')[0]
            ids.append(cell_id)

    return ids

def main(wl_atac_path, wl_rna_path, metadata_paths, metadata_atac_dir, final_data_path, out_path):
    bc_map = load_whitelists(wl_atac_path, wl_rna_path)

    atac_bcs = {}
    metadata_atac_paths = [os.path.join(metadata_atac_dir, i, "metadata.tsv") for i in os.listdir(metadata_atac_dir) if not i.startswith(".")]
    for i in metadata_atac_paths:
        atac_bcs |= load_atac_bcs(i)

    records = {}
    for i in metadata_paths:
        records |= load_records(i, atac_bcs, bc_map)

    final_ids = load_final_data(final_data_path)
    for i in final_ids:
        records[i]["pass_filter"] = True

    with gzip.open(out_path, "wt") as f:
        f.write(HEADER)
        f.write(COLUMNS)

        for cell_id, r in records.items():
            line = f"{cell_id}\t{r['dataset']}\t{r['barcode']}\t{r['dataset_atac']}\t{r['barcode_atac']}\t{r['umi_count']}\tNA\t{r['frac_mito']}\t{r['frac_ribo']}\tNA\t{int(r['pass_filter'])}\n"
            f.write(line)

wl_atac_path = snakemake.input["wl_atac"]
wl_rna_path = snakemake.input["wl_rna"]
metadata_paths = snakemake.input["metadata"]
metadata_atac_dir = snakemake.input["metadata_atac"]
final_data_path = snakemake.input["final_data"]

out_path, = snakemake.output

main(wl_atac_path, wl_rna_path, metadata_paths, metadata_atac_dir, final_data_path, out_path)