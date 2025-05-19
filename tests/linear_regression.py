import pandas as pd
import sys
import statsmodels.api as sm

# --- Load input files from command line ---
feature_file = sys.argv[1]
pheno_file = sys.argv[2]

# --- Load features ---
features = pd.read_csv(feature_file, sep='\t')
features = features.rename(columns={features.columns[0]: "IID"})

# Drop specific columns if they exist
columns_to_drop = [">4220>4222>4223"]
features = features.drop(columns=columns_to_drop, errors='ignore')

# --- Load phenotype ---
pheno = pd.read_csv(pheno_file, sep='\t', usecols=["IID", "PHENO"])

# --- Merge on sample ID ---
df = pd.merge(features, pheno, on="IID")

# Separate target and features
y = df["PHENO"]
X = df.drop(columns=["IID", "PHENO"])

print(X)

# Add intercept
X = sm.add_constant(X)

# Fit OLS model
model = sm.OLS(y, X).fit()

# Print full summary
print(model.summary())

# python3 linear_regression.py ../output/regression/4220_4223.tsv ../data/quantitative/phenotype.tsv
