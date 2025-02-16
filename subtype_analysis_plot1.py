import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

file_path = 'combined_genes_output_mutation.csv'
data = pd.read_csv(file_path)

y_axis_columns = ['CX', 'CR', 'CS', 'F5', 'U5', 'I5', 'IX', 'U3', 'F3', 'R', 'J']

for y_col in y_axis_columns:
    score_hits = data[['mean_score', y_col]]
    
    score_hits_sorted = score_hits.sort_values(by='mean_score', ascending=False).reset_index(drop=True)
    
    score_hits_sorted[y_col] = score_hits_sorted[y_col].cumsum()
    
    if score_hits_sorted[y_col].iloc[-1] == 0:
        print(f"Skipping column {y_col} due to no valid data.")
        continue
    
    x = np.arange(0, len(score_hits_sorted) + 1) / len(score_hits_sorted)
    y = np.concatenate(([0], score_hits_sorted[y_col] / score_hits_sorted[y_col].iloc[-1]))
    
    plt.figure()
    plt.plot(x, y, label='significant')
    plt.plot([0, 1], [0, 1], label='random')
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.xlabel('% of peaks sorted by mean_score')
    plt.ylabel(f'% of {y_col}')
    plt.legend()
    plt.show()
