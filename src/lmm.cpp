#include "lmm.hpp"

void lmm_binary(
    const std::vector<std::vector<size_t>>& df,                  
    const std::vector<bool>& phenotype_binary,      
    const KinshipMatrix& kinship,                                              
    const std::vector<std::vector<double>>& covariates,
    std::string& p_value_str, std::string& beta_str, std::string& se_str, std::string& r2_str) {
    
}

void lmm_quantitative(
    const std::vector<std::vector<size_t>>& df,                  
    const vector<double>& phenotype_table,      
    const KinshipMatrix& kinship,                                              
    const std::vector<std::vector<double>>& covariates,
    std::string& p_value_str, std::string& beta_str, 
    std::string& se_str, std::string& r2_str) {

}

