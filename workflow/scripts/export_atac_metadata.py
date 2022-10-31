import os
import gzip

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
            bc_a_rc = reverse_complement(bc_a)

            bc_map[bc_a] = bc_r
            bc_map[bc_a_rc] = bc_r

    return bc_map


def load_rna_bcs(metadata_rna_path):
    rna_bcs = {}
    with open(metadata_rna_path) as f:
        h = f.readline().rstrip('\n').split('\t')
        barcode_ind = 0
        dataset_ind = h.index("orig.ident")

        for line in f:
            entries = line.rstrip('\n').split('\t')

            dataset = entries[dataset_ind]
            barcode = entries[barcode_ind]

            rna_bcs.setdefault(dataset, set()).add(barcode)
    
    return rna_bcs

def load_records(metadata_path, rna_bcs, bc_map):
    records = {}
    with open(metadata_path) as f:
        h = f.readline().rstrip('\n').split('\t')
        cell_id_ind = h.index("cellNames")
        frag_ind = h.index("nFrags")
        tss_ind = h.index("TSSEnrichment")

        for line in f:
            entries = line.rstrip('\n').split('\t')

            cell_id = entries[cell_id_ind]
            dataset, barcode = cell_id.split("#")
            dataset_parsed = dataset.split("_")[-1].split("-")[0]

            dataset_parsed_rna = dataset.split("_")[0].split("-")[0]
            if dataset_parsed_rna in rna_bcs:
                barcode_trans = bc_map[barcode]
                if barcode_trans in rna_bcs[dataset_parsed_rna]:
                    barcode_rna = barcode_trans
                else:
                    barcode_rna = "NA"
            else:
                dataset_parsed_rna = "NA"
                barcode_rna = "NA"
            
            frag_count = int(entries[frag_ind])
            tss_enr = float(entries[tss_ind])

            record = {
                "dataset": dataset_parsed,
                "dataset_rna": dataset_parsed_rna,
                "barcode": barcode,
                "barcode_rna": barcode_rna,
                "frag_count": frag_count,
                "tss_enr": tss_enr,
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

def main(wl_atac_path, wl_rna_path, metadata_dir, metadata_rna_paths, final_data_path, out_path):
    bc_map = load_whitelists(wl_atac_path, wl_rna_path)

    rna_bcs = {}
    for i in metadata_rna_paths:
        rna_bcs |= load_rna_bcs(i)

    metadata_paths = [os.path.join(metadata_dir, i, "metadata.tsv") for i in os.listdir(metadata_dir) if not i.startswith(".")]
    records = {}
    for i in metadata_paths:
        records |= load_records(i)

    final_ids = load_final_data(final_data_path)
    for i in final_ids:
        records[i]["pass_filter"] = True

    with gzip.open(out_path, "wt") as f:
        f.write(HEADER)
        f.write(COLUMNS)

        for cell_id, r in records.items():
            line = f"{cell_id}\t{r['dataset_rna']}\t{r['barcode_rna']}\t{r['dataset']}\t{r['barcode']}\tNA\t{r['frag_count']}\tNA\tNA\t{r['tss_enr']}\t{int(r['pass_filter'])}\n"
            f.write(line)

wl_atac_path = snakemake.input["wl_atac"]
wl_rna_path = snakemake.input["wl_rna"]
metadata_dir = snakemake.input["metadata"]
metadata_rna_paths = snakemake.input["metadata_rna"]
final_data_path = snakemake.input["final_data"]

out_path, = snakemake.output

main(wl_atac_path, wl_rna_path, metadata_dir, metadata_rna_paths, final_data_path, out_path)