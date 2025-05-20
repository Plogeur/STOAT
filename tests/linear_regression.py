import sys
import pandas as pd
import statsmodels.api as sm

# Predictor matrix X (20 rows, 1 column)
X = np.array([
    [0], [0], [0], [0], [0], [0],
    [1], [1],
    [0.5], [0.5], [0.5], [0.5],
    [1], [1], [1], [1], [1], [1], [1], [1]
])

# --- Merge on sample ID ---
df = pd.merge(features, pheno, on="IID")

# --- Prepare data for model ---
X = df.drop(columns=["IID", "PHENO"])
# X = sm.add_constant(X)  # Adds intercept term
y = df["PHENO"]

# Add intercept
X = sm.add_constant(X)

# Fit OLS model
model = sm.OLS(y, X).fit()

# Print full summary
print(model.summary())

# python3 linear_regression.py ../output/regression/4220_4223.tsv ../data/quantitative/phenotype.tsv

# python3 linear_regression.py ../output/regression/4220_4223.tsv ../data/quantitative/phenotype.tsv
