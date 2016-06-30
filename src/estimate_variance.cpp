#include "scran.h"

SEXP estimate_variance (SEXP qr, SEXP qraux, SEXP exprs, SEXP subset) try {
    matrix_info QR=check_matrix(qr);
    if (QR.is_integer) { 
        throw std::runtime_error("Q matrix should be double-precision");
    }
    matrix_info expr_vals=check_matrix(exprs);
    if (expr_vals.is_integer){ 
        throw std::runtime_error("expression matrix should be double-precision"); 
    }
    if (QR.nrow!=expr_vals.ncol) {
        throw std::runtime_error("'ncol(exprs)' does not equal 'nrow(QR)'");
    }
    if (!isReal(qraux) || size_t(LENGTH(qraux))!=QR.ncol) {
        throw std::runtime_error("QR auxiliary vector should be double-precision and of length 'ncol(Q)'");
    }
    const double* qrxptr=REAL(qraux);
    subset_values subout=check_subset_vector(subset, expr_vals.nrow);
    const int slen=subout.first;
    const int* sptr=subout.second;

    // Setting up for Q-based multiplication.
    run_dormqr multQ(QR.nrow, QR.ncol, QR.dptr, qrxptr, 'T');

    // Setting up pointers to each cell in expr_vals.
    const double** expr_ptrs=(const double**)R_alloc(expr_vals.ncol, sizeof(const double*));
    if (expr_vals.ncol) {
        expr_ptrs[0]=expr_vals.dptr;
        for (size_t c=1; c<expr_vals.ncol; ++c) {
            expr_ptrs[c]=expr_ptrs[c-1]+expr_vals.nrow;
        }
    }

    SEXP output=PROTECT(allocVector(VECSXP, 2));
    try {
        SET_VECTOR_ELT(output, 0, allocVector(REALSXP, slen));
        double* omptr=REAL(VECTOR_ELT(output, 0));
        SET_VECTOR_ELT(output, 1, allocVector(REALSXP, slen));
        double* ovptr=REAL(VECTOR_ELT(output, 1));

        // Running through each gene and reporting its variance and mean.
        size_t c;
        for (int s=0; s<slen; ++s) {
            const int& curgene=sptr[s];
            double& curmean=(omptr[s]=0);
            for (c=0; c<expr_vals.ncol; ++c) {
                const double& curx=expr_ptrs[c][curgene];
                curmean += curx;
                multQ.rhs[c]=curx;
            }
            curmean /= expr_vals.ncol;

            multQ.run();
            double& curvar=(ovptr[s]=0);
            for (c=QR.ncol; c<QR.nrow; ++c) { // only using the residual effects.
                curvar += multQ.rhs[c]*multQ.rhs[c];
            }
            curvar /= QR.nrow - QR.ncol;         
        }
    } catch (std::exception& e){
        UNPROTECT(1);
        throw;
    }

    UNPROTECT(1);
    return output;
} catch (std::exception&e ){
    return mkString(e.what());
}
