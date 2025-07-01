import sys
import pandas as pd
import statsmodels.api as sm

# --- Load input files from command line ---
feature_file = sys.argv[1]
pheno_file = sys.argv[2]

# --- Load features ---
features = pd.read_csv(feature_file, sep='\t')
features = features.rename(columns={features.columns[0]: "IID"})

# --- Load phenotype ---
pheno = pd.read_csv(pheno_file, sep='\t', usecols=["IID", "PHENO"])

# --- Merge on sample ID ---
df = pd.merge(features, pheno, on="IID")

# --- Prepare data for model ---
X = df.drop(columns=["IID", "PHENO"])

# Apply transformation rule to create new column 2 based on column 1
col1 = X.iloc[:, 1]

# Apply your rule
col2 = col1.apply(lambda v: 1 if v == 0 else (0 if v == 1 else 0.5))

# Add col2 as a new feature
X['>7690843>7690845>7690846'] = col2

# Add intercept
X = sm.add_constant(X)

# Target variable
y = df["PHENO"]

print("X : ", X)
print("y : ", y)

# Fit OLS model
model = sm.OLS(y, X).fit()

# Print full summary
print(model.summary())

# python3 linear_regression.py ../output/regression/4220_4223.tsv ../data/quantitative/phenotype.tsv

# python3 linear_regression.py ../output/regression/4220_4223.tsv ../data/quantitative/phenotype.tsv
