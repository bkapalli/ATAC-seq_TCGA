import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

fold_change_dir = '/Volumes/bkapalli/atac_final_figures/BRCA/fold_change'
output_dir = 'ind_plots_fold_change'
os.makedirs(output_dir, exist_ok=True)

# add mutation gene file name + region with desired mutation shaded region and hotspot
# example: 'signal_matrix_gene.csv', ('chr', , , "chr:-")
mutation_files = [
]

dpi_value = 600
left_margin = 0.15

for file_name, region_info in mutation_files:
    chrom, filter_start, filter_end, highlight_region = region_info
    gene_name = file_name.split('_fold_change')[0]
    fold_change_file = os.path.join(fold_change_dir, f"{file_name}")
    if not os.path.exists(fold_change_file):
        print(f"Fold change file for {gene_name} not found. Skipping...")
        continue

    fold_change_df = pd.read_csv(fold_change_file)
    fold_change_df[['chrom', 'coords']] = fold_change_df['region'].str.split(':', expand=True)
    fold_change_df[['chromStart', 'chromEnd']] = fold_change_df['coords'].str.split('-', expand=True).astype(int)
    highlight_start, highlight_end = map(int, highlight_region.split(':')[1].split('-'))
    overlap_regions = fold_change_df[
        (fold_change_df['chromStart'] <= highlight_end) &
        (fold_change_df['chromEnd'] >= highlight_start)
    ]

    plt.figure(figsize=(12, 6), dpi=dpi_value)
    sns.lineplot(data=fold_change_df, x=fold_change_df.index, y='fold_change', color='black', linewidth=3, errorbar=None)

    if not overlap_regions.empty:
        min_overlap_start = overlap_regions.index.min()
        max_overlap_end = overlap_regions.index.max()
        plt.axvspan(min_overlap_start, max_overlap_end, color='#F09D4F', alpha=0.3)

    plt.ylabel('Fold Change', labelpad=15, fontsize=40)  
    ax = plt.gca()
    
    ax.yaxis.set_major_locator(plt.MaxNLocator(3))  
    
    plt.yticks(fontsize=32) 
    
    ax.spines['left'].set_linewidth(2) 
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.spines['bottom'].set_visible(True)
    ax.spines['bottom'].set_linewidth(2) 

    plt.xlabel('Genomic Position', fontsize=40, labelpad=15)
    plt.xticks([]) 

    plt.tight_layout()
    plt.subplots_adjust(left=left_margin)
    plot_filename = os.path.join(output_dir, f"{gene_name}_fold_change_plot.png")
    plt.savefig(plot_filename, dpi=dpi_value, transparent=True)
    plt.close()
