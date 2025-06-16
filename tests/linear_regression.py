import sys
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

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
X = X.drop(df.columns[2], axis=1)

# X = sm.add_constant(X)  # Adds intercept term
y = df["PHENO"]

# Add intercept
X = sm.add_constant(X)

# Fit OLS model
model = sm.OLS(y, X).fit()
# Traceback (most recent call last):
#   File "/home/mbagarre/Bureau/STOAT_CXX/tests/linear_regression.py", line 34, in <module>
#     # Fit OLS model
#             ^^^^^^^^
#   File "/home/mbagarre/miniconda3/lib/python3.12/site-packages/statsmodels/regression/linear_model.py", line 921, in __init__
#     super().__init__(endog, exog, missing=missing,
#   File "/home/mbagarre/miniconda3/lib/python3.12/site-packages/statsmodels/regression/linear_model.py", line 746, in __init__
#     super().__init__(endog, exog, missing=missing,
#   File "/home/mbagarre/miniconda3/lib/python3.12/site-packages/statsmodels/regression/linear_model.py", line 200, in __init__
#     super().__init__(endog, exog, **kwargs)
#   File "/home/mbagarre/miniconda3/lib/python3.12/site-packages/statsmodels/base/model.py", line 270, in __init__
#     super().__init__(endog, exog, **kwargs)
#   File "/home/mbagarre/miniconda3/lib/python3.12/site-packages/statsmodels/base/model.py", line 95, in __init__
#     self.data = self._handle_data(endog, exog, missing, hasconst,
#                 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#   File "/home/mbagarre/miniconda3/lib/python3.12/site-packages/statsmodels/base/model.py", line 135, in _handle_data
#     data = handle_data(endog, exog, missing, hasconst, **kwargs)
#            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#   File "/home/mbagarre/miniconda3/lib/python3.12/site-packages/statsmodels/base/data.py", line 675, in handle_data
#     return klass(endog, exog=exog, missing=missing, hasconst=hasconst,
#            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#   File "/home/mbagarre/miniconda3/lib/python3.12/site-packages/statsmodels/base/data.py", line 84, in __init__
#     self.endog, self.exog = self._convert_endog_exog(endog, exog)
#                             ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#   File "/home/mbagarre/miniconda3/lib/python3.12/site-packages/statsmodels/base/data.py", line 509, in _convert_endog_exog
#     raise ValueError("Pandas data cast to numpy dtype of object. "
# ValueError: Pandas data cast to numpy dtype of object. Check input data with np.asarray(data).

# Print full summary
print(model.summary())

# python3 linear_regression.py ../output/regression/4220_4223.tsv ../data/quantitative/phenotype.tsv

# python3 linear_regression.py ../output/regression/4220_4223.tsv ../data/quantitative/phenotype.tsv
