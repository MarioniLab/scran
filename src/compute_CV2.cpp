#include "scran.h"

/* This computes the mean and CV2 for every gene, while also dividing through by the size factors. 
 * Unless log_prior is specified, in which case the values in 'ptr' are assumed to be log-expressed.
 * These are converted back to the raw scale without any size factor calculations.
 */
template <typename T>
SEXP compute_CV2_internal(const T* ptr, const matrix_info& MAT, SEXP subset_row, SEXP size_factors, SEXP log_prior) {
    subset_values rsubout=check_subset_vector(subset_row, MAT.nrow);
    const int rslen=rsubout.first;
    const int* rsptr=rsubout.second;

    const size_t& ncells=MAT.ncol;
    if (ncells < 2) {
        throw std::runtime_error("need two or more cells to compute variances");
    }

    const bool to_unlog=(log_prior!=R_NilValue);
    const double* sfptr=NULL;
    double lp=0;
    if (to_unlog) { 
        if (!isReal(log_prior) || LENGTH(log_prior)!=1) {
            throw std::runtime_error("prior count should be a double-precision scalar");
        }
        lp=asReal(log_prior);
        if (size_factors!=R_NilValue){ 
            throw std::runtime_error("size factors cannot be specified for log-expression input");
        }
    } else {
        if (!isReal(size_factors)) {
            throw std::runtime_error("size factors should be double-precision");
        } else  if (LENGTH(size_factors)!=int(ncells)) {
            throw std::runtime_error("number of size factors is not equal to number of cells");
        }
        sfptr=REAL(size_factors);
    } 

    // Temporaries.
    double* tmp=(double*)R_alloc(MAT.ncol, sizeof(double));
    size_t col;
    const T* ptr_cpy;
    double diff;
    
    SEXP output=PROTECT(allocVector(VECSXP, 2));
    try {
        SET_VECTOR_ELT(output, 0, allocVector(REALSXP, rslen));
        double* omptr=REAL(VECTOR_ELT(output, 0));
        SET_VECTOR_ELT(output, 1, allocVector(REALSXP, rslen));
        double* ovptr=REAL(VECTOR_ELT(output, 1));
   
        for (int ri=0; ri<rslen; ++ri) {
            ptr_cpy=ptr+rsptr[ri];

            if (!to_unlog) { 
                for (col=0; col<ncells; ++col) {
                    tmp[col]=ptr_cpy[col*MAT.nrow]/sfptr[col];
                }
            } else {
                for (col=0; col<ncells; ++col) {
                    double& current=(tmp[col]=std::pow(2, ptr_cpy[col*MAT.nrow]));
                    current-=lp;
                    if (current < 0) { current=0; }
                }
            }

            // Computing the mean.
            double& curmean=(omptr[ri]=0);
            for (col=0; col<ncells; ++col) {
                curmean+=tmp[col];                    
            }
            curmean/=ncells;
            
            // Computing the variance.
            double& curvar=(ovptr[ri]=0);
            for (col=0; col<ncells; ++col) {
                diff=tmp[col]-curmean;
                curvar+=diff*diff;
            }
            curvar/=(ncells-1);
        }
    } catch (std::exception& e) {
        UNPROTECT(1);
        throw;
    }

    UNPROTECT(1);
    return output;    
}

SEXP compute_CV2(SEXP exprs, SEXP subset_row, SEXP size_factors, SEXP log_prior) try {
    const matrix_info MAT=check_matrix(exprs);
    if (MAT.is_integer) {
        return compute_CV2_internal<int>(MAT.iptr, MAT, subset_row, size_factors, log_prior);
    } else {
        return compute_CV2_internal<double>(MAT.dptr, MAT, subset_row, size_factors, log_prior);
    }
} catch (std::exception& e) {
    return mkString(e.what());
}



