promoter_columns = ['F5', 'U5', 'I5']      # Extended Promoter Region
intronic_columns = ['I5', 'IX']            # Full Intronic Sequence
coding_columns = ['CX', 'CR', 'CS']        # Coding Region

data['Promoter'] = data[promoter_columns].sum(axis=1)
data['Intronic'] = data[intronic_columns].sum(axis=1)
data['Coding'] = data[coding_columns].sum(axis=1)

categories = {'Promoter': 'Promoter Region', 
              'Intronic': 'Intronic Region', 
              'Coding': 'Coding Region'}

plt.figure(figsize=(8, 6))

for category, label in categories.items():
    score_hits_sorted = data[['mean_score', category]].sort_values(by='mean_score', ascending=False).reset_index(drop=True)
    
    score_hits_sorted[category] = score_hits_sorted[category].cumsum()
    
    if score_hits_sorted[category].iloc[-1] == 0:
        print(f"Skipping {label} due to no valid data.")
        continue
    
    x = np.arange(0, len(score_hits_sorted) + 1) / len(score_hits_sorted)
    y = np.concatenate(([0], score_hits_sorted[category] / score_hits_sorted[category].iloc[-1]))
    
    plt.plot(x, y, label=label)

plt.plot([0, 1], [0, 1], 'k--', label='Random')

plt.xlim(0, 1)
plt.ylim(0, 1)
plt.xlabel('% of peaks sorted by mean_score')
plt.ylabel('% of mutations')
plt.title('Cumulative Distribution: Functional Groupings vs. Random')
plt.legend()
plt.tight_layout()
plt.show()
