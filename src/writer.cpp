#include "writer.hpp"

//#define DEBUG_WRITER

namespace stoat_vcf {

void write_binary_covar_header(std::ofstream& outstream) {
    outstream << "CHR\tPOS\tSNARL\tTYPE\tP\tP_ADJUSTED\tBETA\tSE\tALLELE_NUM\tALLELE_PATHS" << endl;
}

void write_binary_header(std::ofstream& outstream) {
    outstream << "CHR\tPOS\tSNARL\tTYPE\tP_FISHER\tP_CHI2\tP_ADJUSTED\tALLELE_NUM\tMIN_ROW_INDEX\tNUM_COLUM\tINTER_GROUP\tAVERAGE\tGROUP_PATHS" << endl;
}

void write_quantitative_header(std::ofstream& outstream) {
    outstream << "CHR\tPOS\tSNARL\tTYPE\tP\tP_ADJUSTED\tRSQUARE\tBETA\tSE\tALLELE_NUM\tALLELE_PATHS" << endl;
}

void write_eqtl_header(std::ofstream& outstream) {
    outstream <<  "CHR\tPOS\tSNARL\tTYPE\tGENE\tP\tP_ADJUSTED\tRSQUARE\tBETA\tSE\tALLELE_NUM\tALLELE_PATHS" << endl;
}

}//end namespace

