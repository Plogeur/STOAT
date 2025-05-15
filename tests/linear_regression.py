# import sys
# import pandas as pd
# import statsmodels.api as sm
# from statsmodels.stats.multitest import multipletests

# # --- Load input files from command line ---
# feature_file = sys.argv[1]
# pheno_file = sys.argv[2]

# # --- Load features ---
# features = pd.read_csv(feature_file, sep='\t')
# features = features.rename(columns={features.columns[0]: "IID"})

# # --- Load phenotype ---
# pheno = pd.read_csv(pheno_file, sep='\t', usecols=["IID", "PHENO"])

# # --- Merge on sample ID ---
# df = pd.merge(features, pheno, on="IID")

# # --- Prepare data for model ---
# X = df.drop(columns=["IID", "PHENO"])
# X = sm.add_constant(X)  # Adds intercept term
# y = df["PHENO"]

# # --- Fit linear regression model ---
# model = sm.OLS(y, X)
# result = model.fit()

# # --- Output coefficients and p-values ---
# print(result.summary())  # Shows coef, std err, t, p-value, conf int

# p_adjusted = multipletests(result.pvalues[1:], method='holm')[1]

# for i in range(len(result.pvalues)-1) :
#     print("p_values : ", result.pvalues[i+1])
#     print("p_adjusted : ", p_adjusted[i])

import numpy as np
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

# --- Define synthetic data (3 features, 5 samples) ---
X = np.array([
    [1, 5, 8],
    [2, 6, 7],
    [3, 5, 6],
    [4, 5, 5],
    [5, 6, 4],
    [6, 5, 3],
    [7, 4, 2],
    [8, 6, 1],
    [9, 5, 0],
    [10, 4, -1],
    [11, 5, -2],
    [12, 4, -3],
    [13, 5, -4],
    [14, 6, -5],
    [15, 5, -6],
    [16, 4, -7],
    [17, 6, -8],
    [18, 5, -9],
    [19, 4, -10],
    [20, 5, -11]
])

y = np.array([
    2.1, 3.2, 4.1, 5.3, 6.0,
    7.2, 8.1, 9.4, 10.5, 11.7,
    12.9, 14.0, 14.8, 16.1, 17.4,
    18.2, 19.6, 20.5, 21.9, 22.7
])

# --- Add intercept ---
X = sm.add_constant(X)

# --- Fit linear regression model ---
model = sm.OLS(y, X)
result = model.fit()

# --- Output full summary ---
print(result.summary())

# --- p-value correction (Holm) ---
raw_pvals = result.pvalues[1:]  # Skip intercept
p_adjusted = multipletests(raw_pvals, method='holm')[1]

# --- Show betas and adjusted p-values ---
for i in range(len(raw_pvals)):
    print(f"Feature x{i+1}: beta = {result.params[i+1]:.18f}, p = {raw_pvals[i]:.18f}, p_adj = {p_adjusted[i]:.18f}")

# python3 linear_regression.py regression_test/quantitative_15_18.tsv regression_test/quantitative_phenotype.tsv