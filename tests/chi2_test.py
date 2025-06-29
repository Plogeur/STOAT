from scipy.stats import chi2_contingency

# Define your datasets
datasets = [
    {"g0": [79, 18], "g1": [96, 23]},
    {"g0": [1, 0], "g1": [0, 1]},
    {"g0": [0, 0], "g1": [0, 1]},
    {"g0": [30, 5], "g1": [2, 25]},
    {"g0": [10, 20], "g1": [20, 10]},
]

# Perform Chi2 test for each dataset
for idx, data in enumerate(datasets, 1):
    table = [data["g0"], data["g1"]]
    print(f"\nDataset {idx}: {table}")

    try:
        chi2, p, dof, expected = chi2_contingency(table, correction=False)
        print(f"  Chi2 Statistic: {chi2:.4f}")
        print(f"  p-value: {p:.4f}")
        print(f"  Degrees of Freedom: {dof}")
        print(f"  Expected Frequencies:\n{expected}")
    except ValueError as e:
        print(f"  Cannot compute Chi2: {e}")
