import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

input_dir = '/Volumes/bkapalli/atac_final_figures/BRCA'
output_dir = 'ind_plots_mean_signal'
os.makedirs(output_dir, exist_ok=True)

# add mutation gene file name + region with desired mutation shaded region and hotspot
# example: 'signal_matrix_gene.csv', ('chr', , , "chr:-")
mutation_files = [
]

subtypes = ['Basal', 'HER2E', 'LumA', 'LumB', 'normal-like']
dpi_value = 600
left_margin = 0.15

for file_name, region_info in mutation_files:
    chrom, filter_start, filter_end, highlight_region = region_info

    basal_df = pd.DataFrame()
    combined_non_basal_df = pd.DataFrame()

    for subtype in subtypes:
        file_path = os.path.join(input_dir, subtype, file_name)
        if not os.path.exists(file_path):
            continue

        df = pd.read_csv(file_path)
        df = df[(df['chrom'] == chrom) & (df['chromStart'] >= filter_start) & (df['chromEnd'] <= filter_end)]

        if subtype == 'Basal':
            basal_df = df
        else:
            combined_non_basal_df = pd.concat([combined_non_basal_df, df])

    if basal_df.empty or combined_non_basal_df.empty:
        print(f"Skipping {file_name} due to missing signal score data.")
        continue

    basal_df = basal_df.drop_duplicates(subset=['chrom', 'chromStart', 'chromEnd'])
    combined_non_basal_df = combined_non_basal_df.drop_duplicates(subset=['chrom', 'chromStart', 'chromEnd'])

    combined_non_basal_df['sqrt_mean_signal_score'] = np.sqrt(combined_non_basal_df['mean_signal_score'])
    basal_df['sqrt_mean_signal_score'] = np.sqrt(basal_df['mean_signal_score'])

    combined_non_basal_df.reset_index(inplace=True)
    basal_df.reset_index(inplace=True)

    highlight_start, highlight_end = map(int, highlight_region.split(':')[1].split('-'))
    overlap_regions = combined_non_basal_df[
        (combined_non_basal_df['chromStart'] <= highlight_end) &
        (combined_non_basal_df['chromEnd'] >= highlight_start)
    ]

    plt.figure(figsize=(12, 6), dpi=dpi_value)

    if not overlap_regions.empty:
        min_overlap_start = overlap_regions['index'].min()
        max_overlap_end = overlap_regions['index'].max()
        plt.axvspan(min_overlap_start, max_overlap_end, color='#F09D4F', alpha=0.3)

    sns.lineplot(data=combined_non_basal_df, x='index', y='sqrt_mean_signal_score', color='blue', linewidth=3, errorbar=None, label='Non-Basal')
    sns.lineplot(data=basal_df, x='index', y='sqrt_mean_signal_score', color='green', linewidth=3, errorbar=None, label='Basal')

    plt.legend(loc='upper right', bbox_to_anchor=(1, 1), fontsize=35, ncol=1, frameon=False)

    ax = plt.gca()
    ax.spines['left'].set_linewidth(2)
    ax.yaxis.set_major_locator(plt.MaxNLocator(3))

    min_signal_score = min(basal_df['sqrt_mean_signal_score'].min(), combined_non_basal_df['sqrt_mean_signal_score'].min())
    max_signal_score = max(basal_df['sqrt_mean_signal_score'].max(), combined_non_basal_df['sqrt_mean_signal_score'].max())
    ax.set_ylim(min_signal_score, max_signal_score)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)  
    ax.spines['bottom'].set_linewidth(2)  
    
    plt.xlabel('Genomic Position', fontsize=40, labelpad=15)  
    plt.xticks([])  
    ax.set_xlim(left=combined_non_basal_df['index'].min(), right=combined_non_basal_df['index'].max())

    plt.ylabel('ATAC-seq Mean', labelpad=15, fontsize=40)
    plt.yticks(fontsize=32)


    plt.tight_layout()
    plt.subplots_adjust(left=left_margin)

    plot_filename = os.path.join(output_dir, f"{file_name.split('.')[0]}_mean_signal_plot.png")
    plt.savefig(plot_filename, dpi=dpi_value, transparent=True)
    plt.close()
