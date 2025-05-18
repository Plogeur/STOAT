import numpy as np
import statsmodels.api as sm

# Predictor matrix X (20 rows, 1 column)
X = np.array([
    [0], [0], [0], [0], [0], [0],
    [1], [1],
    [0.5], [0.5], [0.5], [0.5],
    [1], [1], [1], [1], [1], [1], [1], [1]
])

# Response vector y (20 values)
y = np.array([0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0])

# Add intercept
X = sm.add_constant(X)

# Fit OLS model
model = sm.OLS(y, X).fit()

# Print full summary
print(model.summary())

# python3 linear_regression.py