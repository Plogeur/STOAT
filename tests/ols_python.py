import pandas as pd # type: ignore
import statsmodels.api as sm # type: ignore

def linear_regression(df:pd.DataFrame, pheno:dict) :

    df = df.astype(int)
    df['Target'] = df.index.map(pheno)
    x = df.drop('Target', axis=1)
    y = df['Target']

    result = sm.OLS(y, x).fit()
    print(result.summary())

    print("R-squared:", result.rsquared)  # R-squared value
    print("F-statistic:", result.fvalue)  # F-statistic
    print("p-value of F-statistic:", result.f_pvalue)  # p-value of the F-statistic

def main():
    # Using same data as in C++ example
    data = {
        'sample1': [1, 0, 1],
        'sample2': [0, 1, 1],
        'sample3': [1, 1, 0],
        'sample4': [0, 0, 1],
        'sample5': [1, 1, 1],
    }

    df = pd.DataFrame(data).T  # Transpose to have samples as index
    df.index.name = 'Sample'  # Name the index for clarity

    # Same phenotype values
    pheno = {
        'sample1': 2.5,
        'sample2': 1.8,
        'sample3': 3.0,
        'sample4': 1.2,
        'sample5': 3.5
    }

    # Run regression
    linear_regression(df, pheno)

if __name__ == "__main__":
    main()