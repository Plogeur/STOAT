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
X = sm.add_constant(X)  # Adds intercept term
y = df["PHENO"]

# --- Fit linear regression model ---
model = sm.OLS(y, X)
result = model.fit()

# --- Output coefficients and p-values ---
print(result.summary())  # Shows coef, std err, t, p-value, conf int
