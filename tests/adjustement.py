from statsmodels.stats.multitest import multipletests

# Same p-values as in C++
p_values = [0.01, 0.04, 0.03, 0.002, 0.05]

# Apply Holm-Bonferroni
_, p_adjusted, _, _ = multipletests(p_values, method='holm')

# Print
print("Original p-values:")
print([round(p, 4) for p in p_values])
print("\nHolm-adjusted p-values:")
print([round(p, 4) for p in p_adjusted])
