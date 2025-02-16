# remove 0 byte gene files prior
# atac gene file matched to respectivegenomic signal gene mutation files 

import os
import pandas as pd
import re

DIR1 = "/temp_work/ch262369/BRCA/gene_files"
DIR2 = "/temp_work/ch262369/genomic_signal/extracted_files/out/Breast"
OUTPUT_DIR = "specific_genes_mutation"
os.makedirs(OUTPUT_DIR, exist_ok=True)
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "combined_genes_output.csv")
LOG_MATCHED = os.path.join(OUTPUT_DIR, "matched_genes.log")
LOG_UNMATCHED = os.path.join(OUTPUT_DIR, "unmatched_genes.log")

with open(LOG_MATCHED, "w") as matched, open(LOG_UNMATCHED, "w") as unmatched:
    matched.write("Matched Genes:\n")
    unmatched.write("Unmatched Genes:\n")

def get_exact_match_file(directory, gene):
    for f in os.listdir(directory):
        if re.match(rf"^{gene}(_|$)", f, re.IGNORECASE):
            return os.path.join(directory, f)
    return None

def load_gene_files(gene):
    csv_file = get_exact_match_file(DIR1, gene)
    bed_file = get_exact_match_file(DIR2, gene)
    return csv_file, bed_file

def clean_column_types(data, numeric_columns):
    for col in numeric_columns:
        data[col] = pd.to_numeric(data[col], errors="coerce").astype("Int64")
    return data.dropna()

def merge_data(csv_file, bed_file):
    csv_data = pd.read_csv(csv_file)
    csv_data = clean_column_types(csv_data, ["chromStart", "chromEnd"])
    bed_data = pd.read_csv(bed_file, sep="\t", names=[
        "chrom", "chromStart", "chromEnd", "CX", "CR", "CS", "F5", "U5", "I5", "IX", "U3", "F3", "R", "J"
    ])
    bed_data = clean_column_types(bed_data, ["chromStart", "chromEnd"])
    csv_data["chromStart"] = csv_data["chromStart"].astype("int64")
    bed_data["chromStart"] = bed_data["chromStart"].astype("int64")
    merged = pd.merge_asof(
        csv_data.sort_values("chromStart"),
        bed_data.sort_values("chromStart"),
        left_on="chromStart",
        right_on="chromStart",
        direction="nearest",
        tolerance=5
    )
    return merged.dropna()

all_results = []
for gene in os.listdir(DIR1):
    gene_name = gene.split('_')[0]
    csv_file, bed_file = load_gene_files(gene_name)
    if csv_file and bed_file:
        try:
            merged_data = merge_data(csv_file, bed_file)
            if not merged_data.empty:
                all_results.append(merged_data.assign(gene=gene_name))
                with open(LOG_MATCHED, "a") as matched:
                    matched.write(f"{gene_name}\n")
            else:
                with open(LOG_UNMATCHED, "a") as unmatched:
                    unmatched.write(f"{gene_name}: No matching rows in merged data\n")
        except Exception as e:
            with open(LOG_UNMATCHED, "a") as unmatched:
                unmatched.write(f"{gene_name}: Error during merge ({e})\n")
    else:
        with open(LOG_UNMATCHED, "a") as unmatched:
            unmatched.write(f"{gene_name}: Files not found ({'CSV missing' if not csv_file else ''}{', ' if not csv_file and not bed_file else ''}{'BED missing' if not bed_file else ''})\n")

if all_results:
    pd.concat(all_results).to_csv(OUTPUT_FILE, index=False)
