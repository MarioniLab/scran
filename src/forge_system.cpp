#include "scran.h"

/*** A function to (a) subset by row, (b) subset by column, and (c) divide through by the library sizes. 
 *** The output is equivalent to t(t(MAT[row_subset,col_subset])/lib.sizes) where lib.sizes itself is
 *** computed as colSums(MAT[row_subset,col_subset]).
 ***/
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

/*** A function to estimate the pooled size factors and construct the linear equations. ***/
SEXP forge_system (SEXP exprs, SEXP ordering, SEXP size, SEXP ref) try {
    // Checking inputs
    const matrix_info emat=check_matrix(exprs);
    const int& ncells=emat.ncol;
    const int& ngenes=emat.nrow;
    if (ncells==0) { throw std::runtime_error("at least one cell required for normalization"); }
    const double** eptrs=(const double**)R_alloc(ncells, sizeof(const double*));
    eptrs[0]=emat.dptr;
    int cell=0;
    for (cell=1; cell<ncells; ++cell) { eptrs[cell]=eptrs[cell-1]+ngenes; }

    if (!isInteger(ordering)) { throw std::runtime_error("ordering vector should be integer"); }
    if (LENGTH(ordering)<ncells*2-1)  { throw std::runtime_error("ordering vector is too short for number of cells"); }
    const int* orptr=INTEGER(ordering);

    if (!isInteger(size) || LENGTH(size) > 1) { throw std::runtime_error("size should be an integer scalar"); }
    const int SIZE=asInteger(size);
    if (SIZE < 1 || SIZE > ncells) { throw std::runtime_error("size should be within [1, number of cells]"); }

    if (!isNumeric(ref)) { throw std::runtime_error("reference expression vector should be double-precision"); }
    const double* rptr=REAL(ref);
    if (ngenes!=LENGTH(ref)) { throw std::runtime_error("length of reference vector is inconsistent with number of cells"); }

    // Setting up the output matrix.
    SEXP output=PROTECT(allocVector(VECSXP, 3));
try { 
    SET_VECTOR_ELT(output, 0, allocVector(INTSXP, SIZE * ncells));
    int* row_optr=INTEGER(VECTOR_ELT(output, 0));
    SET_VECTOR_ELT(output, 1, allocVector(INTSXP, SIZE * ncells));
    int* col_optr=INTEGER(VECTOR_ELT(output, 1));
    SET_VECTOR_ELT(output, 2, allocVector(REALSXP, ncells));
    double* ofptr=REAL(VECTOR_ELT(output, 2));

    // Running through the ordering.
    int index=0, gene=0;
    const bool is_even=bool(ngenes%2==0);
    const int halfway=int(ngenes/2);
    double* combined=(double*)R_alloc(ngenes, sizeof(double));
    const double* cur_eptr;
    const int* cur_window;
    
    for (cell=0; cell<ncells; ++cell) {
        std::fill(combined, combined+ngenes, 0);
        std::fill(row_optr, row_optr+SIZE, cell);
        row_optr+=SIZE;
        cur_window=orptr+cell;

        for (index=0; index<SIZE; ++index) {
            const int& curcell=cur_window[index];
            *(col_optr++) = curcell;
            cur_eptr=eptrs[curcell];

            for (gene=0; gene<ngenes; ++gene) { 
                combined[gene]+=cur_eptr[gene];
            }
        }

        // Computing the median ratio against the reference.
        for (gene=0; gene<ngenes; ++gene) { 
            combined[gene]/=rptr[gene];
        }
        std::partial_sort(combined, combined+halfway+1, combined+ngenes);
        if (is_even) {
            ofptr[cell]=(combined[halfway]+combined[halfway-1])/2;
        } else {
            ofptr[cell]=combined[halfway];
        }
    }    
} catch (std::exception& e) {
    UNPROTECT(1);
    return mkString(e.what());
}
    UNPROTECT(1);
    return output;
} catch (std::exception& e) {
    return mkString(e.what());
}

