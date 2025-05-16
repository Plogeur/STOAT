import numpy as np
import statsmodels.api as sm

# --- Synthetic Data ---
X = np.array([
    [1,  0],
    [1,  0],
    [1,  0],
    [1,  0],
    [1,  0],
    [1,  0],
    [0, 1],
    [0, 1],
    [0.5,  0.5],
    [0.5,  0.5],
    [0.5,  0.5],
    [0.5,  0.5],
    [0, 1],
    [0, 1],
    [0, 1],
    [0, 1],
    [0, 1],
    [0, 1],
    [0, 1],
    [0, 1],
])

phenotype = np.array([0,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,1,1,0,0])

# Add intercept
# X = sm.add_constant(X)

# Fit OLS model
model = sm.OLS(phenotype, X).fit()

# Print full summary
print(model.summary())
