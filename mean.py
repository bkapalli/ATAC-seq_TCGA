# bigwig to wig (for TCGA)
# calculate mean score of all samples for all genes 

import os
import pandas as pd
import numpy as np
import time

unfiltered_wig_folder = "/temp_work/ch262369/BRCA/wig_files"
gene_file_path = "/home/ch262369/Desktop/ncbi_gene_hg38.tsv"
output_directory = "/temp_work/ch262369/BRCA/gene_files"
log_file = "/temp_work/ch262369/BRCA/processed_genes.log"

os.makedirs(output_directory, exist_ok=True)

bin_size = 100

gene_df = pd.read_csv(gene_file_path, sep='\t')

if not os.path.exists(log_file):
    with open(log_file, "w") as f:
        f.write("Processed Genes Log\n")

def log_progress(message):
    with open(log_file, "a") as f:
        f.write(message + "\n")

def process_wig_file(file):
    wig = pd.read_csv(file, sep='\t', comment='#', header=None, 
                      names=['chrom', 'chromStart', 'chromEnd', 'score'], dtype={'chrom': str})
    wig = wig.dropna(subset=['chromStart', 'chromEnd', 'score'])
    
    relevant_genes = gene_df[gene_df['chrom'].isin(wig['chrom'].unique())]
    
    for _, gene_row in relevant_genes.iterrows():
        chrom, flank_start, flank_end, gene_name = gene_row['chrom'], gene_row['flankStart'], gene_row['flankEnd'], gene_row['gene']
        output_file = os.path.join(output_directory, f"{gene_name}_{chrom}_{flank_start}_{flank_end}_mean_signal.csv")
        
        if os.path.exists(output_file):
            log_progress(f"Skipping {gene_name} ({chrom}:{flank_start}-{flank_end}), file already exists.")
            continue
        
        region_range = pd.DataFrame({
            'chromStart': np.arange(flank_start, flank_end, bin_size),
            'chromEnd': np.arange(flank_start + bin_size, flank_end + bin_size, bin_size)
        })
        region_range['score_sum'] = 0
        region_range['count'] = 0
        
        gene_wig_data = wig[(wig['chrom'] == chrom) & (wig['chromEnd'] > flank_start) & (wig['chromStart'] < flank_end)]
        
        for _, row in gene_wig_data.iterrows():
            start, end, score = max(row['chromStart'], flank_start), min(row['chromEnd'], flank_end), row['score']
            bins = np.arange(start, end, bin_size)
            for bin_start in bins:
                bin_end = min(bin_start + bin_size, end)
                overlap = bin_end - bin_start
                idx = region_range['chromStart'] == bin_start
                region_range.loc[idx, 'score_sum'] += score * overlap
                region_range.loc[idx, 'count'] += overlap
        
        region_range['mean_score'] = region_range['score_sum'] / region_range['count']
        region_range['mean_score'] = region_range['mean_score'].fillna(0)
        region_range[['chromStart', 'chromEnd', 'mean_score']].to_csv(output_file, index=False)
        log_progress(f"Finished processing {gene_name} ({chrom}:{flank_start}-{flank_end})")

start_time = time.time()
for file in os.listdir(unfiltered_wig_folder):
    if file.endswith('.wig'):
        process_wig_file(os.path.join(unfiltered_wig_folder, file))

print(f"All gene regions processed and saved in {output_directory} in {time.time() - start_time:.4f} seconds")
