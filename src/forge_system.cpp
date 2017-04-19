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
SEXP forge_system (SEXP exprs, SEXP ordering, SEXP sizes, SEXP ref) try {
    // Checking input matrix.
    const matrix_info emat=check_matrix(exprs);
    const int ncells=emat.ncol;
    const int ngenes=emat.nrow;
    if (ncells==0) { throw std::runtime_error("at least one cell required for normalization"); }

    std::vector<const double*> eptrs(ncells);
    eptrs[0]=emat.dptr;
    for (int cell=1; cell<ncells; ++cell) { eptrs[cell]=eptrs[cell-1]+ngenes; }

    // Checking the input sizes.
    if (!isInteger(sizes) || LENGTH(sizes)==0) { throw std::runtime_error("sizes should be a non-empty integer vector"); }
    const int nsizes=LENGTH(sizes);
    const int * szptr=INTEGER(sizes);
    int s, total_SIZE=0;
    for (s=0; s<nsizes; ++s) { 
        const int& SIZE=szptr[s];
        if (SIZE < 1 || SIZE > ncells) { throw std::runtime_error("each element of sizes should be within [1, number of cells]"); }
        if (s!=0 && SIZE < szptr[s-1]) { throw std::runtime_error("sizes should be sorted"); }
        total_SIZE+=SIZE;
    }

    // Checking reference and ordering.
    if (!isNumeric(ref)) { throw std::runtime_error("reference expression vector should be double-precision"); }
    const double* rptr=REAL(ref);
    if (ngenes!=LENGTH(ref)) { throw std::runtime_error("length of reference vector is inconsistent with number of cells"); }

    if (!isInteger(ordering)) { throw std::runtime_error("ordering vector should be integer"); }
    if (LENGTH(ordering)<ncells*2-1)  { throw std::runtime_error("ordering vector is too short for number of cells"); }
    const int* orptr=INTEGER(ordering);

    // Setting up the output matrix.
    SEXP output=PROTECT(allocVector(VECSXP, 3));
try { 
    SET_VECTOR_ELT(output, 0, allocVector(INTSXP, total_SIZE * ncells));
    int * row_optr=INTEGER(VECTOR_ELT(output, 0));

    SET_VECTOR_ELT(output, 1, allocVector(INTSXP, total_SIZE * ncells));
    int * col_optr=INTEGER(VECTOR_ELT(output, 1));

    SET_VECTOR_ELT(output, 2, allocVector(REALSXP, nsizes * ncells));
    double* ofptr=REAL(VECTOR_ELT(output, 2));

    int index=0, gene=0;
    std::vector<double> combined(ngenes), ratios(ngenes);
    const double* cur_eptr;
    const int* cur_window;
    int rownum;

    const bool is_even=bool(ngenes%2==0);
    const int halfway=int(ngenes/2);
    double medtmp;
    
    // Running through the sliding windows.
    for (int win=0; win<ncells; ++win) {
        std::fill(combined.begin(), combined.end(), 0);
        cur_window=orptr+win;

        index=0;
        for (s=0; s<nsizes; ++s) {
            const int& SIZE=szptr[s];
            rownum=ncells*s + win; // Setting the row so that all rows with the same SIZE are consecutive.
            std::fill(row_optr, row_optr+SIZE, rownum);
            row_optr+=SIZE;

            while (index<SIZE) {
                const int& curcell=cur_window[index];
                col_optr[index] = curcell;
                cur_eptr=eptrs[curcell];
                for (gene=0; gene<ngenes; ++gene) { 
                    combined[gene]+=cur_eptr[gene];
                }
                ++index;
            }
            
            col_optr+=SIZE;
            if (s+1!=nsizes) {
                // Copying over to the next set of column assignments.
                std::copy(col_optr-SIZE, col_optr, col_optr);
            }

            // Computing the ratio against the reference.
            for (gene=0; gene<ngenes; ++gene) { 
                ratios[gene]=combined[gene]/rptr[gene];
            }

            // Computing the median (faster than partial sort).
            std::nth_element(ratios.begin(), ratios.begin()+halfway, ratios.end());
            if (is_even) {
                medtmp=ratios[halfway];
                std::nth_element(ratios.begin(), ratios.begin()+halfway-1, ratios.end());
                ofptr[rownum]=(medtmp+ratios[halfway-1])/2;
            } else {
                ofptr[rownum]=ratios[halfway];
            }       
        }

/*
        std::partial_sort(combined, combined+halfway+1, combined+ngenes);
        if (is_even) {
            ofptr[cell]=(combined[halfway]+combined[halfway-1])/2;
        } else {
            ofptr[cell]=combined[halfway];
        } 
*/
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

