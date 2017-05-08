#include "scran.h"

/* This function computes residuals in a nice and quick manner.
 * It takes a pre-computed QR matrix from qr(), and optionally
 * a subset of integer indices, and returns a matrix of residuals.
 */

SEXP get_residuals(SEXP exprs, SEXP qr, SEXP qraux, SEXP subset, SEXP lower_bound) try {
    const matrix_info emat=check_matrix(exprs);
    if (emat.is_integer) {
        throw std::runtime_error("expression matrix must be double-precision"); 
    }
    const size_t& ncells=emat.ncol;
    std::vector<const double*> eptrs(ncells);
    if (ncells) {
        eptrs[0]=emat.dptr;
        for (size_t c=1; c<ncells; ++c) {
            eptrs[c]=eptrs[c-1]+emat.nrow;
        }
    }

    // Checking the subset vector.
    subset_values subout=check_subset_vector(subset, emat.nrow);
    const int slen=subout.first;
    const int* sptr=subout.second;
    
    // Checking the QR matrix.
    const matrix_info QR=check_matrix(qr);
    if (QR.is_integer){ 
        throw std::runtime_error("QR matrix must be double-precision");
    }
    if (!isReal(qraux) || size_t(LENGTH(qraux))!=QR.ncol) {
        throw std::runtime_error("QR auxiliary vector should be double-precision and of length 'ncol(Q)'");
    }
    const double* qrxptr=REAL(qraux);
    run_dormqr multQ1(QR.nrow, QR.ncol, QR.dptr, qrxptr, 'T');
    run_dormqr multQ2(QR.nrow, QR.ncol, QR.dptr, qrxptr, 'N');

    // Checking the lower bound.
    if (!isReal(lower_bound) || LENGTH(lower_bound)!=1) {
        throw std::runtime_error("lower bound should be a numeric scalar");
    }
    const double lbound=asReal(lower_bound);
    const bool check_lower=R_FINITE(lbound);
    
    SEXP output=PROTECT(allocMatrix(REALSXP, slen, ncells));
    try {
        std::vector<double*> optrs(ncells);
        if (ncells) {
            optrs[0]=REAL(output);
            for (size_t c=1; c<ncells; ++c) {
                optrs[c]=optrs[c-1]+slen;
            }
        }

        size_t c;
        double* temporary=(double*)R_alloc(ncells, sizeof(double));
        std::deque<int> below_bound;
        double lowest;

        for (int s=0; s<slen; ++s) {
            for (c=0; c<ncells; ++c) {
                temporary[c]=eptrs[c][sptr[s]];
            }
            
            // Identifying elements below the lower bound.
            if (check_lower) { 
                for (c=0; c<ncells; ++c) {
                    if (temporary[c] <= lbound) {
                        below_bound.push_back(c);                        
                    }
                }
            }

            multQ1.run(temporary); // Getting main+residual effects.
            for (c=0; c<QR.ncol; ++c) {
                temporary[c]=0; // setting main effects to zero.
            }
            multQ2.run(temporary); // Getting residuals.

            // Forcing the values below the boundary to a value below the smallest residual.
            if (check_lower && !below_bound.empty()) {
                lowest=(*std::min_element(temporary, temporary+ncells)) - 1;
                for (c=0; c<below_bound.size(); ++c) {
                    temporary[below_bound[c]]=lowest;
                }
                below_bound.clear();
            }

            for (c=0; c<ncells; ++c) {
                optrs[c][s]=temporary[c];
            }
        }
    } catch (std::exception& e) {
        UNPROTECT(1);
        throw;
    }

    UNPROTECT(1);
    return output;
} catch (std::exception& e) {
    return mkString(e.what());
}

