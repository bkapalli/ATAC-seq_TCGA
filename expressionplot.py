import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def add_jitter(values, scale=0.1):
    return values + np.random.uniform(-scale, scale, size=len(values))

# Define gene file location with atac and expression data 
basal_file = '/Volumes/bkapalli/atac_final_figures/BRCA/expression/updated_basal_signal_matrix_with_gene.txt'
nonbasal_file = '/Volumes/bkapalli/atac_final_figures/BRCA/expression/updated_nonbasal_signal_matrix_with_gene.txt'

basal_df = pd.read_csv(basal_file, sep='\t')
nonbasal_df = pd.read_csv(nonbasal_file, sep='\t')

# add mutation region
chrom = 
start = 
end = 

chr10_region_basal = basal_df[
    (basal_df['chrom'] == chrom) &
    (basal_df['chromStart'] < end) &
    (basal_df['chromEnd'] > start)
]

if chr10_region_basal.empty:
    raise ValueError(f"No overlapping rows found for {chrom}:{start}-{end} in the basal file.")

chr10_region_basal_mean = chr10_region_basal.iloc[:, 3:].mean(axis=0)

chr10_region_nonbasal = nonbasal_df[
    (nonbasal_df['chrom'] == chrom) &
    (nonbasal_df['chromStart'] < end) &
    (nonbasal_df['chromEnd'] > start)
]

if chr10_region_nonbasal.empty:
    raise ValueError(f"No overlapping rows found for {chrom}:{start}-{end} in the non-basal file.")

chr10_region_nonbasal_mean = chr10_region_nonbasal.iloc[:, 3:].mean(axis=0)

x_basal_atac = np.sqrt(chr10_region_basal_mean.values)
x_nonbasal_atac = np.sqrt(chr10_region_nonbasal_mean.values)

y_basal_expression = pd.to_numeric(basal_df.iloc[-1, 3:].values, errors='coerce')
y_nonbasal_expression = pd.to_numeric(nonbasal_df.iloc[-1, 3:].values, errors='coerce')

combined_expression = np.concatenate([y_basal_expression, y_nonbasal_expression])

normalized_expression = 100 * (combined_expression - combined_expression.min()) / (combined_expression.max() - combined_expression.min())

y_basal_normalized = normalized_expression[:len(y_basal_expression)]
y_nonbasal_normalized = normalized_expression[len(y_basal_expression):]

y_basal_normalized_jittered = add_jitter(y_basal_normalized)
y_nonbasal_normalized_jittered = add_jitter(y_nonbasal_normalized)

plt.figure(figsize=(12, 6))
plt.scatter(x_basal_atac, y_basal_normalized_jittered, color='green', label='Basal', alpha=0.7, s=100)
plt.scatter(x_nonbasal_atac, y_nonbasal_normalized_jittered, color='blue', label='Non-Basal', alpha=0.7, s=100)
plt.xlabel('ATAC-seq at hotspot', fontsize=40, labelpad=30)
plt.ylabel('Expression %', fontsize=40)

x_min = min(np.min(x_basal_atac), np.min(x_nonbasal_atac))
x_max = max(np.max(x_basal_atac), np.max(x_nonbasal_atac))

x_min = int(np.floor(x_min))
x_max = int(np.ceil(x_max))

plt.xticks([0, x_max], fontsize=32)
plt.yticks([0, 100], fontsize=32)

plt.legend().remove()

ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)

plt.gca().patch.set_alpha(0)
plt.tight_layout()

output_file = "expression_plot_gene.png"
plt.savefig(output_file, format='png', dpi=300, transparent=True)
plt.show()
