import pandas as pd
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create a covariate file based on a phenotype file.")
    
    parser.add_argument("--pheno", required=True, help="Path to the phenotype file (must contain 'IID' column).")
    parser.add_argument("--out", required=True, help="Path to save the output covariate file.")
    
    # Covariate values (constant for all samples)
    parser.add_argument("--sex", type=int, required=True, help="Sex covariate value.")
    parser.add_argument("--cp1", type=float, required=True, help="CP1 covariate value.")
    parser.add_argument("--cp2", type=float, required=True, help="CP2 covariate value.")
    parser.add_argument("--cp3", type=float, required=True, help="CP42 covariate value.")
    
    args = parser.parse_args()

    # Read phenotype file
    phenos = pd.read_csv(args.pheno, sep="\t")

    if "IID" not in phenos.columns:
        raise ValueError("Phenotype file must contain an 'IID' column.")

    # Create covariate dataframe
    covariates = phenos[["IID"]].copy()
    covariates["SEX"] = args.sex
    covariates["CP1"] = args.cp1
    covariates["CP2"] = args.cp2
    covariates["CP3"] = args.cp3

    # Save covariate file
    covariates.to_csv(args.out, sep="\t", index=False)
    print(f"Covariate file created successfully: {args.out}")

# python3 covar.py --pheno ../data/binary/phenotype.tsv --sex 0 --cp1 25.215 --cp2 75.84 --cp3 45 --out ../data/binary/covariate.tsv
# python3 covar.py --pheno ../data/quantitative/phenotype.tsv --sex 0 --cp1 1.25 --cp2 0.0045 --cp3 0.5 --out ../data/quantitative/covariate.tsv
