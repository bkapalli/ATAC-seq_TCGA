import os
import pandas as pd

# Depending on cancer type change the name
# Directory containing the WIG files
dir_path = "/Volumes/Elements/BCHinternship/Bulk_ATAC/lung_cancer/LUSC/wig_files"
filtered_wig_folder = "Filtered_WIG"
processed_folder = "Processed_CSV"
os.makedirs(filtered_wig_folder, exist_ok=True)
os.makedirs(processed_folder, exist_ok=True)

# Define genomic region of interest
chrom = 
start_region = 
end_region = 

files = [os.path.join(dir_path, f) for f in os.listdir(dir_path) if f.endswith('.wig')]

# Process each WIG file
def read_and_filter_wig(file, chrom, start_region, end_region):
    print(f"Processing file: {file} for region {chrom}:{start_region}-{end_region}")
    wig = pd.read_csv(file, sep='\t', header=None, comment='#', 
                      names=['chrom', 'chromStart', 'chromEnd', 'score'], dtype={'chrom': str})
    
    wig['chromStart'] = pd.to_numeric(wig['chromStart'], errors='coerce')
    wig['chromEnd'] = pd.to_numeric(wig['chromEnd'], errors='coerce')
    wig = wig.dropna(subset=['chromStart', 'chromEnd'])
    
    filtered_wig = wig[(wig['chrom'] == chrom) & (wig['chromStart'] >= start_region) & (wig['chromEnd'] <= end_region)]
    
    output_file = os.path.join(filtered_wig_folder, os.path.basename(file).replace(".wig", "_filtered.wig"))
    filtered_wig.to_csv(output_file, sep='\t', index=False, header=False)
    
    return filtered_wig

# Initialize list to store data
all_data = []

for file in files:
    sample_name = os.path.basename(file).replace(".wig", "")
    filtered_wig = read_and_filter_wig(file, chrom, start_region, end_region)
    
    all_data.append(filtered_wig[['chrom', 'chromStart', 'chromEnd', 'score']])

# Concatenate all the data and compute the mean of the signal scores
merged_df = pd.concat(all_data, axis=0)
merged_df = merged_df.groupby(['chrom', 'chromStart', 'chromEnd'], as_index=False).mean()

# Rename the columns to include the mean signal score
merged_df = merged_df.rename(columns={'score': 'mean_signal_score'})

# Save the merged dataframe
output_file = os.path.join(processed_folder, "gene_ATAC_mean_signal.csv")
merged_df.to_csv(output_file, index=False)
print(f"Merged ATAC mean signal file saved: {output_file}")

