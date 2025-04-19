import scipy.stats as stats
import numpy as np

def run_tests(table):
    """
    Perform Chi-squared test on 2xN table.
    Perform Fisher's exact test only if table is 2x2.

    Parameters:
    table (list of lists): A 2xN contingency table.

    Returns:
    dict: Dictionary containing p-values for chi2 and optionally fisher tests
    """
    table_array = np.array(table)

    if table_array.shape[0] != 2:
        raise ValueError("Table must have exactly 2 rows (2xN format)")

    if np.any(table_array.sum(axis=0) == 0) or np.any(table_array.sum(axis=1) == 0):
        raise ValueError("Table has a full zero row or column, which is not allowed")

    # Chi-squared test
    chi2_stat, chi2_p, _, _ = stats.chi2_contingency(table_array)

    # Fisher's exact test only for 2x2
    if table_array.shape == (2, 2):
        fisher_oddsratio, fisher_p = stats.fisher_exact(table_array)
    else:
        fisher_p = None

    return {
        "chi2_pvalue": chi2_p,
        "fisher_pvalue": fisher_p
    }

def print_results(name, table):
    print(f"\n{name}")
    print("Contingency Table:")
    for row in table:
        print(row)
    try:
        results = run_tests(table)
        print("Chi-squared test p-value:", results['chi2_pvalue'])
        if results['fisher_pvalue'] is not None:
            print("Fisher's exact test p-value:", results['fisher_pvalue'])
        else:
            print("Fisher's exact test not applicable (table not 2x2)")
    except Exception as e:
        print("Error:", e)

if __name__ == "__main__":
    examples = {
        "Example 1: 2x2 Balanced": [[10, 20], [20, 10]],
        "Example 2: 2x2 Strong effect": [[30, 5], [2, 25]],
        "Example 3: 2x3 Table": [[10, 15, 5], [20, 10, 10]],
        "Example 4: 2x4 Table": [[5, 10, 15, 20], [20, 15, 10, 5]],
        "Example 5: 2x5 Uniform": [[10, 10, 10, 10, 10], [10, 10, 10, 10, 10]],
        "Example 6: All Zeros": [[0, 0], [0, 0]],
        "Example 6: All Zeros": [[0, 0], [0, 1]],
        "Example 7: Full Zero Row": [[0, 0, 0], [10, 20, 30]],
        "Example 8: Full Zero Column": [[0, 10, 5], [0, 20, 15]],
        "Example 9: One Cell Non-Zero": [[1, 0], [0, 1]],
    }

    for name, table in examples.items():
        print_results(name, table)
