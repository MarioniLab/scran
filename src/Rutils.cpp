#include "scran.h"

/* Various utility functions that don't really warrant their own page, 
 * but are helpful for avoiding construction of large temporaries. */

/* A function to (a) subset by row, (b) subset by column, and (c) divide through by the library sizes. 
 * The output is equivalent to t(t(MAT[row_subset,col_subset])/lib.sizes) where lib.sizes itself is
 * computed as colSums(MAT[row_subset,col_subset]).
 */

template <typename T>
SEXP subset_and_divide_internal(const T* ptr, const matrix_info& MAT, SEXP row_subset, SEXP col_subset) {
    // Checking row subset vector
    subset_values rsubout=check_subset_vector(row_subset, int(MAT.nrow));
    const int rslen=rsubout.first;
    const int* rsptr=rsubout.second;

    // Checking column subset vector
    subset_values csubout=check_subset_vector(col_subset, int(MAT.ncol));
    const int cslen=csubout.first;
    const int* csptr=csubout.second;

    SEXP output=PROTECT(allocVector(VECSXP, 2));
    try {
        SET_VECTOR_ELT(output, 0, allocVector(REALSXP, cslen));
        double* olptr=REAL(VECTOR_ELT(output, 0));

        SET_VECTOR_ELT(output, 1, allocMatrix(REALSXP, rslen, cslen));
        double* onptr=REAL(VECTOR_ELT(output, 1));
        const T* curptr;

        for (int cs=0; cs<cslen; ++cs) {
            curptr=ptr + MAT.nrow*csptr[cs];

            double& curlib=(olptr[cs]=0);
            for (int rs=0; rs<rslen; ++rs) { 
                curlib+=curptr[rsptr[rs]];
            }
            if (curlib < 0.00000001) {
                throw std::runtime_error("cells should have non-zero library sizes");
            }

            for (int rs=0; rs<rslen; ++rs) { 
                onptr[rs] = curptr[rsptr[rs]]/curlib;
            }
            onptr+=rslen;
        }

    } catch (std::exception& e) {
        UNPROTECT(1);
        throw;
    }

    UNPROTECT(1);
    return output;
}

SEXP subset_and_divide(SEXP matrix, SEXP row_subset, SEXP col_subset) try {
    matrix_info MAT=check_matrix(matrix);
    if (MAT.is_integer){
        return subset_and_divide_internal<int>(MAT.iptr, MAT, row_subset, col_subset);
    } else {
        return subset_and_divide_internal<double>(MAT.dptr, MAT, row_subset, col_subset);
    }
} catch (std::exception& e) {
    return mkString(e.what());
}

/* This function computes error-tolerant ranks for a subset of genes in a subset of cells. */

struct data_holder {
    data_holder(int nvals) : nobs(nvals) {
        holder1=PROTECT(allocVector(REALSXP, nobs));
        holder2=PROTECT(allocVector(REALSXP, nobs));
        arg1=REAL(holder1);
        arg2=REAL(holder2);
        index=(int*)R_alloc(nobs, sizeof(int));
        return;
    };

    ~data_holder() {
        UNPROTECT(2);
        return;
    }

    const int nobs;
    int* index; 
    SEXP holder1, holder2;
    double* arg1, *arg2;
};


template <typename T>
SEXP rank_subset_internal (const T* ptr, const matrix_info& MAT, SEXP subset_row, SEXP subset_col, SEXP tol) {
    if (!isReal(tol) || LENGTH(tol)!=1) {
        throw std::runtime_error("tolerance must be a double-precision scalar");
    }
    const T tolerance=asReal(tol);
    subset_values rsubout=check_subset_vector(subset_row, MAT.nrow);
    const int rslen=rsubout.first;
    const int* rsptr=rsubout.second;
    subset_values csubout=check_subset_vector(subset_col, MAT.ncol);
    const int cslen=csubout.first;
    const int* csptr=csubout.second;
    
    // Setting up some pointers to the matrix.
    const T** ptrs=(const T**)R_alloc(MAT.ncol, sizeof(const T*));
    if (MAT.ncol) {
        ptrs[0]=ptr;
        for (size_t c=1; c<MAT.ncol; ++c) {
            ptrs[c]=ptrs[c-1]+MAT.nrow;
        }
    }

    SEXP output=PROTECT(allocMatrix(INTSXP, cslen, rslen));
    try {
        if (!cslen || !rslen) {
            UNPROTECT(1);
            return output;
        }
        int* optr=INTEGER(output);
        data_holder dh(cslen);
        Rx_random_seed myseed;

        int cs;
        double last_unique;
        for (int rs=0; rs<rslen; ++rs) {
            for (cs=0; cs<cslen; ++cs) {
                dh.index[cs]=cs;
                dh.arg1[cs]=ptrs[csptr[cs]][rsptr[rs]];
            }

            // First stage sorting and equalization of effective ties.
            R_orderVector1(dh.index, dh.nobs, dh.holder1, FALSE, FALSE);
            last_unique=dh.arg1[dh.index[0]]; // Should be okay, as we remove cases where cslen=0.
            for (cs=1; cs<cslen; ++cs) {
                double& val=dh.arg1[dh.index[cs]];
                if (val - last_unique <= tolerance) {
                    val=last_unique;
                } else {
                    last_unique=val;
                }
            }

            // Second stage sorting with broken ties.
            for (cs=0; cs<cslen; ++cs) {
                dh.index[cs]=cs;
                dh.arg2[cs]=unif_rand();
            }
            R_orderVector(dh.index, dh.nobs, Rf_lang2(dh.holder1, dh.holder2), FALSE, FALSE);

            // Filling the output matrix.
            for (cs=0; cs<cslen; ++cs){ 
                optr[dh.index[cs]]=cs+1;
            }
            optr+=cslen;
        }
    } catch (std::exception& e) {
        UNPROTECT(1);
        throw;
    }

    UNPROTECT(1);
    return output;
}

SEXP rank_subset(SEXP exprs, SEXP subset_row, SEXP subset_col, SEXP tol) try {
    const matrix_info MAT=check_matrix(exprs);
    if (MAT.is_integer) {
        return rank_subset_internal<int>(MAT.iptr, MAT, subset_row, subset_col, tol);
    } else {
        return rank_subset_internal<double>(MAT.dptr, MAT, subset_row, subset_col, tol);
    }
} catch (std::exception& e) {
    return mkString(e.what());
}

/* This function computes residuals in a nice and quick manner. */

SEXP get_residuals(SEXP exprs, SEXP qr, SEXP qraux, SEXP subset) try {
    const matrix_info emat=check_matrix(exprs);
    if (emat.is_integer) {
        throw std::runtime_error("expression matrix must be double-precision"); 
    }
    const double** eptrs=(const double**)R_alloc(emat.ncol, sizeof(const double*));
    if (emat.ncol) {
        eptrs[0]=emat.dptr;
        for (size_t c=1; c<emat.ncol; ++c) {
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
    
    SEXP output=PROTECT(allocMatrix(REALSXP, slen, emat.ncol));
    try {
        double** optrs=(double**)R_alloc(emat.ncol, sizeof(double*));
        if (emat.ncol) {
            optrs[0]=REAL(output);
            for (size_t c=1; c<emat.ncol; ++c) {
                optrs[c]=optrs[c-1]+slen;
            }
        }

        size_t c;
        double* temporary=(double*)R_alloc(emat.ncol, sizeof(double));
        for (int s=0; s<slen; ++s) {
            for (c=0; c<emat.ncol; ++c) {
                temporary[c]=eptrs[c][sptr[s]];
            }

            multQ1.run(temporary); // Getting main+residual effects.
            for (c=0; c<QR.ncol; ++c) {
                temporary[c]=0; // setting main effects to zero.
            }

            multQ2.run(temporary); // Getting residuals.
            for (c=0; c<emat.ncol; ++c) {
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

