#include "scran.h"

/* This computes the mean and CV2 for every gene, while also dividing through by the size factors. */
template <typename T>
SEXP compute_CV2_internal(const T* ptr, const matrix_info& MAT, SEXP subset_row, SEXP size_factors) {
    subset_values rsubout=check_subset_vector(subset_row, MAT.nrow);
    const int rslen=rsubout.first;
    const int* rsptr=rsubout.second;

    const size_t& ncells=MAT.ncol;
    if (ncells < 2) {
        throw std::runtime_error("need two or more cells to compute variances");
    }
    if (!isReal(size_factors)) {
        throw std::runtime_error("size factors should be double-precision");
    } else  if (LENGTH(size_factors)!=int(ncells)) {
        throw std::runtime_error("number of size factors is not equal to number of cells");
    }
    const double* sfptr=REAL(size_factors);

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
            for (col=0; col<ncells; ++col) {
                tmp[col]=ptr_cpy[col*MAT.nrow]/sfptr[col];
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

SEXP compute_CV2(SEXP exprs, SEXP subset_row, SEXP size_factors) try {
    const matrix_info MAT=check_matrix(exprs);
    if (MAT.is_integer) {
        return compute_CV2_internal<int>(MAT.iptr, MAT, subset_row, size_factors);
    } else {
        return compute_CV2_internal<double>(MAT.dptr, MAT, subset_row, size_factors);
    }
} catch (std::exception& e) {
    return mkString(e.what());
}



