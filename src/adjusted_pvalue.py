import pandas as pd
import numpy as np
import argparse

# Function to read the GWAS file
def read_gwas_file(file_path):
    df = pd.read_csv(file_path, sep='\t')
    return df

# Function to combine P_FISHER and P_CHI2
def combine_pvalues(df):
    return df[['P_FISHER', 'P_CHI2']].mean(axis=1)

# Function to apply Benjamini-Hochberg correction
def apply_bh_correction(df, combined_pvalues):
    p_adjusted = combined_pvalues * len(df) / (df.index + 1)
    return p_adjusted

# Function to format p-values (smaller memory usage & better visibility)
def format_pvalue(p_value):
    if np.isnan(p_value):  
        return "NA"  # Keep empty if NaN
    return f"{p_value:.4e}" if p_value < 0.0001 else f"{p_value:.4f}"

# Function to modify the GWAS DataFrame and add the adjusted p-value
def modify_gwas_file(df, type_):
    if type_ == "binary":
        # sort df p-value
        df = df.sort_values('P_CHI2')
        df = df.reset_index(drop = True)
        combined_pvalues = combine_pvalues(df)
        p_adjusted = apply_bh_correction(df, combined_pvalues)
        df.insert(df.columns.get_loc('P_CHI2') + 1, 'P_ADJUSTED', p_adjusted)  # Insert column
        df['P_FISHER'] = df['P_FISHER'].apply(format_pvalue)
        df['P_CHI2'] = df['P_CHI2'].apply(format_pvalue)

    else:
        df = df.sort_values('P')
        df = df.reset_index(drop = True)
        combined_pvalues = df['P']
        p_adjusted = apply_bh_correction(df, combined_pvalues)
        df.insert(df.columns.get_loc('P') + 1, 'P_ADJUSTED', p_adjusted)  
        df['P'] = df['P'].apply(format_pvalue)

    # Apply formatting to reduce memory usage
    df['P_ADJUSTED'] = df['P_ADJUSTED'].apply(format_pvalue)
    return df

# Function to save the modified GWAS file (in-place update)
def save_modified_file(df, file_path):
    df.to_csv(file_path, sep='\t', index=False)  # Save the file in place

# Function to create a file with top variants based on adjusted p-value < 10^-5
def save_top_variants(df, file_path):
    top_variants_df = df[df['P_ADJUSTED'].apply(lambda x: x != "NA" and float(x) < 1e-5)]  # Filter top variants
    top_variants_df.to_csv(file_path, sep='\t', index=False)
    print(f"Top variants saved to: {file_path}")

# Main function to execute the process
def main(input_file_path, type_):
    df = read_gwas_file(input_file_path)
    modified_df = modify_gwas_file(df, type_)
    save_modified_file(modified_df, input_file_path)  
    
    # Save top variants (p-value < 10^-5)
    top_variants_file = input_file_path.replace(".tsv", "_top_variants.tsv")
    save_top_variants(modified_df, top_variants_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the Pvalue Stoat GWAS analysis")
    parser.add_argument("--pvalue", type=str, required=True, help="Path to the GWAS pvalue file")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-b", "--binary", action="store_true", help="Specify binary type analysis")
    group.add_argument("-q", "--quantitative", action="store_true", help="Specify quantitative type analysis")
    
    args = parser.parse_args()

    type_ = "binary" if args.binary else "quantitative"
    main(args.pvalue, type_)
