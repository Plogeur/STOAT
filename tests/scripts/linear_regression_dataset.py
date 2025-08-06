import numpy as np
import statsmodels.api as sm

def linear_regression_summary(df, quantitative_phenotype):
    X = np.array(df)
    y = np.array(quantitative_phenotype)

    # Add intercept
    X = sm.add_constant(X)

    try:
        model = sm.OLS(y, X)
        results = model.fit()

        # Return full summary text
        return results.summary()

    except Exception as e:
        return f"Error during regression: {e}"

# Test cases data
test_cases = [
    {
        "name": "Linear Regression 1 - Perfect Linear Relationship",
        "df": [
            [0],
            [1],
            [0]
        ],
        "quantitative_phenotype": [2.0, 4.0, 6.0]
    },
    {
        "name": "Linear Regression 2 - Moderate",
        "df": [
            [0.5, 0],
            [0, 0.5],
            [1, 0],
            [0, 1],
            [0, 0.5]
        ],
        "quantitative_phenotype": [10.5, 13.0, 15.8, 19.7, 21.5]
    },
    {
        "name": "Linear Regression 3 - Weaker Correlation",
        "df": [
                [1, 0],
                [1, 0],
                [1, 0],
                [1, 0],
                [1, 0],
                [1, 0],
                [1, 0],
                [0, 1],
                [0, 0]
        ],
        "quantitative_phenotype": [4.5, 7.0, 9.2, 10.9, 13.0, 14.0, 11.0, 15.0, 16.0]
    }
]

# Run and print summary for each test
for test in test_cases:
    print(f"Test: {test['name']}")
    summary = linear_regression_summary(test["df"], test["quantitative_phenotype"])
    print(summary)
    print("="*80)

# python3 linear_regression_dataset.py


#                             OLS Regression Results                            
# ==============================================================================
# Dep. Variable:                      y   R-squared:                       0.000
# Model:                            OLS   Adj. R-squared:                 -1.000
# Method:                 Least Squares   F-statistic:                     0.000
# Date:                Wed, 06 Aug 2025   Prob (F-statistic):               1.00
# Time:                        16:15:46   Log-Likelihood:                -5.7281
# No. Observations:                   3   AIC:                             15.46
# Df Residuals:                       1   BIC:                             13.65
# Df Model:                           1                                         
# Covariance Type:            nonrobust                                         
# ==============================================================================
#                  coef    std err          t      P>|t|      [0.025      0.975]
# ------------------------------------------------------------------------------
# const          4.0000      2.000      2.000      0.295     -21.412      29.412
# x1          1.332e-15      3.464   3.85e-16      1.000     -44.016      44.016
# ==============================================================================
# Omnibus:                          nan   Durbin-Watson:                   1.000
# Prob(Omnibus):                    nan   Jarque-Bera (JB):                0.281
# Skew:                           0.000   Prob(JB):                        0.869
# Kurtosis:                       1.500   Cond. No.                         2.41