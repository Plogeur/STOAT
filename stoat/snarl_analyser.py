import argparse
from stoat import utils
from cyvcf2 import VCF # type: ignore
import numpy as np # type: ignore
import pandas as pd # type: ignore
import statsmodels.api as sm # type: ignore
from scipy.stats import chi2_contingency # type: ignore
from scipy.stats import fisher_exact # type: ignore
from typing import Optional, List
# from limix.stats import logisticMixedModel # type: ignore
# from limix.stats import scan # type: ignore
import subprocess
import time
import os

class Matrix :
    def __init__(self, default_row_number:int=1_000_000, column_number:int=2):
        self.default_row_number = default_row_number 
        self.matrix = np.zeros((default_row_number, column_number), dtype=bool)
        self.row_header = None

    def get_matrix(self) -> np.zeros:
        return self.matrix

    def set_matrix(self, expended_matrix:int) -> None:
        self.matrix = expended_matrix

    def get_row_header(self) -> Optional[dict]:
        return self.row_header

    def get_default_row_number(self) -> int:
        return self.default_row_number

    def set_row_header(self, row_header:dict) -> None:
        self.row_header = row_header  

    def add_data(self, idx_snarl:int, idx_geno:int) -> None:
        self.matrix[idx_snarl, idx_geno] = 1

class SnarlProcessor:
    def __init__(self, vcf_path:str, list_samples:list):
        self.list_samples = list_samples
        self.matrix = Matrix(self._count_lines_with_wc(vcf_path)*4, len(self.list_samples)*2)
        self.vcf_path = vcf_path

    def _count_lines_with_wc(self, file_path):
        result = subprocess.run(['wc', '-l', file_path], stdout=subprocess.PIPE, text=True)
        line_count = int(result.stdout.split()[0])
        return line_count

    def expand_matrix(self):
        """Expands a given numpy matrix by doubling the number of rows."""
        data_matrix = self.matrix.get_matrix()
        current_rows, current_cols = data_matrix.shape
        new_rows = current_rows + self.matrix.get_default_row_number() # add + default_row_number row  

        # Create a new matrix of zeros with the expanded size
        expanded_matrix = np.zeros((new_rows, current_cols), dtype=data_matrix.dtype)
        expanded_matrix[:current_rows, :] = data_matrix
        self.matrix.set_matrix(expanded_matrix)

    def determine_str(self, s:str, length_s:int, i:int) -> tuple[int, int]:
        """Extract an integer from a string starting at index i."""
        start_idx = i
        while i < length_s and s[i] not in ['>', '<']:
            i += 1
        return i, s[start_idx:i]

    def decompose_string(self, s:str) -> List[str]:
        """Decompose a string with snarl information."""
        result = []
        i = 0
        length_s = len(s)
        prev_int = None
        prev_sym = None
        
        while i < length_s:
            start_sym = s[i]
            i += 1
            i, current_int = self.determine_str(s, length_s, i)

            if prev_int is not None and prev_sym is not None:
                result.append(f"{prev_int}")
                result.append(f"{current_int}")
            
            prev_int = current_int
            prev_sym = start_sym
        
        return result

    def decompose_snarl(self, lst:List[str]) -> List[List[str]]:
        """Decompose a list of snarl strings."""
        return [self.decompose_string(s) for s in lst]

    def get_or_add_index(self, row_header_dict:dict, key:str, length_ordered_dict:int) -> int:
        """ 
        Retrieve the index of the key if it exists in the OrderedDict.
        If the key does not exist, add it and return the new index.
        """
        if key in row_header_dict:
            return row_header_dict[key]
        else:
            new_index = length_ordered_dict
            row_header_dict[key] = new_index
            return new_index
    
    def push_matrix(self, idx_snarl:int, decomposed_snarl:str, row_header_dict:dict, index_column:int) -> None:
        """Add True to the matrix if snarl is found"""

        # Retrieve or add the index in one step and calculate the length once
        length_ordered_dict = len(row_header_dict)
        idx_snarl = self.get_or_add_index(row_header_dict, decomposed_snarl, length_ordered_dict)

        # Check if a new matrix chunk is needed (only if length > 1)
        current_rows_number = self.matrix.get_matrix().shape[0]
        if length_ordered_dict > current_rows_number -1 :
            self.expand_matrix()

        # Add data to the matrix
        self.matrix.add_data(idx_snarl, index_column)

    def fill_matrix(self) -> None:
        """Parse VCF file (main function)"""
        row_header_dict = dict()

        # Parse variant line by line
        for variant in VCF(self.vcf_path):
            genotypes = variant.genotypes  # Extract genotypes once per variant
            
            if variant.INFO.get('LV') > 0 :
                continue

            path_list = variant.INFO.get('AT', '').split(',')  # Extract and split snarl list once per variant
            list_list_decomposed_node = self.decompose_snarl(path_list)  # Decompose snarls once per variant

            for index_column, genotype in enumerate(genotypes) :

                allele_1, allele_2 = genotype[:2]  # assume there are only 2 allele
                col_idx = index_column * 2

                if allele_1 == -1 or allele_2 == -1 :# case where we got ./.
                    continue

                for decompose_allele_1 in list_list_decomposed_node[allele_1] :
                    self.push_matrix(allele_1, decompose_allele_1, row_header_dict, col_idx)

                for decompose_allele_2 in list_list_decomposed_node[allele_2] :
                    self.push_matrix(allele_2, decompose_allele_2, row_header_dict, col_idx + 1)

        self.matrix.set_row_header(row_header_dict)

    def binary_table(self, snarls:dict, binary_groups:tuple[dict, dict], kinship_matrix:pd.DataFrame=None, covar:Optional[dict]=None, gaf:bool=False, output:str="output/binary_output.tsv"):
        """
        Generate a binary table with statistical results and write to a file.
        """
        
        common_headers = (
            "CHR\tPOS\tNODE\tREF\tALT\tP_FISHER\tP_CHI2\tALLELE_NUM\tMIN_ROW_INDEX\tNUM_COLUM\tINTER_GROUP\tAVERAGE"
        )
        headers = f"{common_headers}\tGROUP_PATHS\n" if gaf else f"{common_headers}\n"

        with open(output, 'wb') as outf:
            outf.write(headers.encode('utf-8'))

            for node, node_info in snarls.items() :
                chromosome, position = node_info
                # Create the binary table, considering covariates if provided
                df = self.create_binary_table(binary_groups, node)
                # if kinship_matrix and covar :
                #     p_value, beta, vcomp = lmm_pvalue = self.LMM_binary(df, kinship_matrix, covar)
                ref = alt = 'NA'
                # Perform statistical tests and compute descriptive statistics
                fisher_p_value, chi2_p_value, allele_number, min_sample, numb_colum, inter_group, average, group_paths = self.binary_stat_test(df, gaf)
                common_data = (
                f"{chromosome}\t{position}\t{node}\t{ref}\t{alt}\t"
                f"{fisher_p_value}\t{chi2_p_value}\t{allele_number}\t{min_sample}\t"
                f"{numb_colum}\t{inter_group}\t{average}")
                data = f"{common_data}\t{group_paths}\n" if gaf else f"{common_data}\n"
            
                outf.write(data.encode('utf-8'))

    def quantitative_table(self, snarls:dict, quantitative_dict:dict, kinship_matrix:pd.DataFrame=None, covar:Optional[dict]=None, output:str="output/quantitative_output.tsv") :

        with open(output, 'wb') as outf:
            headers = 'CHR\tPOS\tNODE\tREF\tALT\tRSQUARED\tBETA\tSE\tP\tALLELE_NUM\n'
            outf.write(headers.encode('utf-8'))
            for node, node_info in snarls.items() :
                chromosome, position = node_info
                df, allele_number = self.create_quantitative_table(node)
                rsquared, beta, se, pvalue = self.linear_regression(df, quantitative_dict)
                ref = alt = 'NA'
                data = (
                f"{chromosome}\t{position}\t{node}\t{ref}\t{alt}\t"
                f"{rsquared}\t{beta}\t{se}\t"
                f"{pvalue}\t{allele_number}\n")
                outf.write(data.encode('utf-8'))

    def identify_correct_path(self, node:str) -> list:
        """
        Return a list of column indices where all specific elements of this column in the matrix are 1.
        """

        try :
            row_index = self.matrix.get_row_header()[node]
        except :
            return []

        # Extract the rows from the matrix using rows_to_check
        extracted_rows = self.matrix.get_matrix()[row_index, :]
        return extracted_rows

    def create_binary_table(self, binary_groups:dict, node:str) -> pd.DataFrame:
        """Generates a binary table DataFrame indicating the presence of snarl paths in given groups based on matrix data"""

        # Initialize g0 and g1 with zeros, corresponding to the length of column_headers
        g0 = [0, 0]
        g1 = [0, 0]

        # Count occurrences in g0 and g1
        extracted_rows = self.identify_correct_path(node)

        for idx, allele in enumerate(extracted_rows):
            if allele :
                srr = self.list_samples[idx // 2]
                if binary_groups[srr] == 0:
                    g0[0] += 1
                else :
                    g1[0] += 1

        g0[1] = len(self.list_samples)*2 - g0[0]
        g1[1] = len(self.list_samples)*2 - g1[0]

        # Create and return the DataFrame
        df = pd.DataFrame([g0, g1], index=['G0', 'G1'], columns=[node, f'not_passing_{node}'])
        return df

    def create_quantitative_table(self, node: str) -> pd.DataFrame:
        """
        Creates a quantitative table for a given node, showing passing and not passing alleles,
        and calculates the total allele count.
        """

        length_sample = len(self.list_samples)
        extracted_rows = self.identify_correct_path(node)

        # Initialize genotype counts for 'passing' and 'not passing' alleles
        passing = [0] * length_sample

        # Count passing alleles
        for idx, allele in enumerate(extracted_rows):
            if allele:
                passing[idx // 2] += 1

        # Calculate 'not passing' alleles
        # not_passing = [2 - count for count in passing]

        # Create a DataFrame
        df = pd.DataFrame(
            {'Passing': passing},
            index=self.list_samples)

        # Calculate total allele count
        allele_number = sum(passing) # + sum(not_passing)
        return df, allele_number

    def linear_regression(self, df:pd.DataFrame, pheno:dict) -> tuple :

        df = df.astype(int)
        df['Target'] = df.index.map(pheno)
        x = df.drop('Target', axis=1)
        y = df['Target']
        result = sm.OLS(y, x).fit()

        rsquared = f"{result.rsquared:.4e}" if result.rsquared < 0.0001 else f"{result.rsquared:.4f}"

        # Mean of beta coefficients
        beta_mean = f"{result.params.mean():.4e}" if result.params.mean() < 0.0001 else f"{result.params.mean():.4f}"

        # Mean of standard errors
        se_mean = f"{result.bse.mean():.4e}" if result.bse.mean() < 0.0001 else f"{result.bse.mean():.4f}"
        formatted_p_value = f"{result.f_pvalue:.4e}" if result.f_pvalue < 0.0001 else f"{result.f_pvalue:.4f}"

        return rsquared, beta_mean, se_mean, formatted_p_value

    # def LMM_quantitatif(self, df:pd.DataFrame, kinship_matrix:pd.DataFrame, covar:dict, pheno:dict) -> tuple:
    #     """
    #     Perform Linear Mixed Model (LMM) for quantitative phenotype data.
    #     """

    #     # Ensure the covariate matrix is in a DataFrame and map the covariates and phenotype correctly
          # df = df.astype(int)
          # df['Target'] = df.index.map(pheno)
          # x = df.drop('Target', axis=1)
          # y = df['Target']

          # x_with_const = sm.add_constant(x)

    #     # Perform Linear Mixed Model scan (assuming a function like scan)
    #     results = scan(y=y, K=kinship_matrix, covariates=x)
        
    #     # Extract metrics from the results object (p-value, beta, beta_se, log-likelihood, heritability)
    #     p_value = round(results.stats["pv"], 4) # P-values for each covariate
    #     beta = results.stats["beta"]            # Effect sizes (coefficients for covariates)
    #     beta_se = results.stats["beta_se"]      # Standard errors for effect sizes
    #     ll = results.stats["ll"]                # Log-likelihood of the model
    #     heritability = results.stats["h2"]      # Heritability estimate (proportion of variance explained by GRM)
    #     return p_value, beta, beta_se, ll, heritability

    def chi2_test(self, df:pd.DataFrame) -> str:
        """Calculate p_value using chi-2 test"""

        # Check if dataframe has at least 2 columns and more than 0 counts in every cell
        if df.shape[1] >= 2 and np.all(df.sum(axis=0)) and np.all(df.sum(axis=1)):
            # Perform Chi-Square test
            p_value = chi2_contingency(df)[1] # from scipy.stats import chi2_contingency
            p_value = f"{p_value:.4e}" if p_value < 0.0001 else f"{p_value:.4f}"

        else:
            p_value = "NA"

        return p_value

    def fisher_test(self, df:pd.DataFrame) -> str:
        """Calcul p_value using fisher exact test"""

        try:
            p_value = fisher_exact(df)[1] # from scipy.stats import fisher_exact
            p_value = f"{p_value:.4e}" if p_value < 0.0001 else f"{p_value:.4f}"

        except ValueError as e:
            p_value = 'NA'
        
        return p_value

    # # Logistic Mixed Model
    # def LMM_binary(self, df, kinship_matrix, covar):
    #     """
    #     Perform Logistic Mixed Model on a binary phenotype.
    #     """
    #     # Map phenotype to df
    #     df['Target'] = df.index.map(pheno)
    #     y = df['Target'].values
    #     X = covar[df.index].values  # Covariates should match the index of genotype data
    #     lmm = logisticMixedModel(y=y, K=kinship_matrix, X=X)
        
    #     # Fit the model with covariates
    #     lmm.fit(X)
    #     beta = lmm.beta                 # Effect sizes (log-odds for covariates)
    #     p_value = round(lmm.pv, 4)      # P-values for fixed effects
    #     vcomp = lmm.vcomp               # Variance components (relatedness)

    #     return p_value, beta, vcomp
    
    def format_group_paths(self, df:pd.DataFrame) :
        """Format group paths as a string for adding gaf column information."""
        result = []
        for column in df.columns:
            column_values = [f"{df.loc[group, column]}" for group in df.index]
            result.append(":".join(column_values))

        final_str = ",".join(result)
        return final_str
    
    def binary_stat_test(self, df:pd.DataFrame, gaf:bool=False) -> tuple:
        """ Perform statistical tests and calculate descriptive statistics on the binary analysis."""
        fisher_p_value = self.fisher_test(df)
        chi2_p_value = self.chi2_test(df)

        allele_number = int(df.iloc[:, 0].sum()) # Calculate the number of samples in the DataFrame
        inter_group = int(df.iloc[:, 0].min()) # Calculate the sum of the minimum values from each column
        numb_colum = df.shape[1] # Get the total number of columns in the DataFrame
        average = float(allele_number / numb_colum) # Calculate the average of the total sum divided by the number of columns
        row_sums = df.sum(axis=1) # Calculate the sum of values for each row in the DataFrame
        min_row_index = row_sums.min() # Find the minimum sum among the row sums
        group_paths = self.format_group_paths(df) if gaf else ""

        return fisher_p_value, chi2_p_value, allele_number, min_row_index, numb_colum, inter_group, average, group_paths
    
if __name__ == "__main__" :
    parser = argparse.ArgumentParser(description="Parse and analyse snarl from vcf file")
    parser.add_argument("vcf_path", type=utils.check_format_vcf_file, help="Path to the vcf file (.vcf or .vcf.gz)")
    parser.add_argument("snarl", type=utils.check_format_list_path, help="Path to the snarl file that containt snarl and aT (.txt or .tsv)")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-b", "--binary", type=utils.check_format_pheno, help="Path to the binary phenotype file (.txt or .tsv)")
    group.add_argument("-q", "--quantitative", type=utils.check_format_pheno, help="Path to the quantitative phenotype file (.txt or .tsv)")
    parser.add_argument("-c", "--covariate", type=utils.check_covariate_file, required=False, help="Path to the covariate file (.txt or .tsv)")
    parser.add_argument("-o", "--output", type=str, required=False, help="Path to the output file")
    args = parser.parse_args()

    start = time.time()
    list_samples = utils.parsing_samples_vcf(args.vcf_path)
    vcf_object = SnarlProcessor(args.vcf_path, list_samples)
    vcf_object.fill_matrix()
    print(f"Time Matrix : {time.time() - start} s")

    start = time.time()
    snarl = utils.parse_snarl_path_file_dict(args.snarl)
    covar = utils.parse_covariate_file(args.covariate) if args.covariate else None

    output_dir = args.output or "output"    
    os.makedirs(output_dir, exist_ok=True)
    output = os.path.join(output_dir, "stoat2bin.assoc.tsv")

    if args.binary:
        binary_group = utils.parse_pheno_binary_file(args.binary)
        vcf_object.binary_table(snarl, binary_group, covar, output=output)

    # python3 stoat/snarl_analyser.py tests/simulation/binary_data/merged_output.vcf tests/simulation/binary_data/snarl_paths.tsv -b tests/simulation/binary_data/phenotype.tsv -o tests/binary_tests_output

    if args.quantitative:
        quantitative_dict = utils.parse_pheno_quantitatif_file(args.quantitative)
        vcf_object.quantitative_table(snarl, quantitative_dict, covar, output=output)

    # python3 stoat/snarl_analyser.py tests/simulation/quantitative_data/merged_output.vcf tests/simulation/quantitative_data/snarl_paths.tsv -q tests/simulation/quantitative_data/phenotype.tsv -o tests/quantitative_tests_output

    print(f"Time P-value : {time.time() - start} s")


