#include "scran.h"

SEXP forge_system (SEXP ng, SEXP nc, SEXP exprs, SEXP ordering, SEXP size, SEXP ref) try {
    // Checking inputs
    if (!isInteger(ng) || LENGTH(ng) > 1) { throw std::runtime_error("number of genes must be a positive integer"); }
    if (!isInteger(nc) || LENGTH(nc) > 1) { throw std::runtime_error("number of cells must be a positive integer"); }
    const int ncells=asInteger(nc), ngenes=asInteger(ng);
    if (ncells==0) { throw std::runtime_error("at least one cell required for normalization"); }
    
    if (!isInteger(ordering)) { throw std::runtime_error("ordering vector should be integer"); }
    if (LENGTH(ordering)<ncells*2-1)  { throw std::runtime_error("ordering vector is too short for number of cells"); }
    const int* orptr=INTEGER(ordering);

    if (!isInteger(size) || LENGTH(size) > 1) { throw std::runtime_error("size should be an integer scalar"); }
    const int SIZE=asInteger(size);
    if (SIZE < 1 || SIZE > ncells) { throw std::runtime_error("size should be within [1, number of cells]"); }

    if (!isNumeric(ref)) { throw std::runtime_error("reference expression vector should be double-precision"); }
    const double* rptr=REAL(ref);
    if (ngenes!=LENGTH(ref)) { throw std::runtime_error("length of reference vector is inconsistent with number of cells"); }

    if (!isNumeric(exprs)) { throw std::runtime_error("expression matrix should be double-precision"); }
    if (LENGTH(exprs)!=ncells*ngenes) { throw std::runtime_error("matrix dimensions are inconsistent with the number of genes/cells"); }
    const double** eptrs=(const double**)R_alloc(ncells, sizeof(const double*));
    eptrs[0]=REAL(exprs);
    int cell=0;
    for (cell=1; cell<ncells; ++cell) { eptrs[cell]=eptrs[cell-1]+ngenes; }

    // Setting up the output matrix.
    SEXP output=PROTECT(allocVector(VECSXP, 2));
try { 
    SET_VECTOR_ELT(output, 0, allocMatrix(REALSXP, ncells, ncells));
    double** oxptrs=(double**)R_alloc(ncells, sizeof(double*));
    oxptrs[0]=REAL(VECTOR_ELT(output, 0));
    for (cell=1; cell<ncells; ++cell) { oxptrs[cell]=oxptrs[cell-1]+ncells; } 
    std::fill(oxptrs[0], oxptrs[0] + ncells*ncells, 0);

    SET_VECTOR_ELT(output, 1, allocVector(REALSXP, ncells));
    double* ofptr=REAL(VECTOR_ELT(output, 1));

    // Running through the ordering.
    int index=0, gene=0;
    const bool is_even=bool(ngenes%2==0);
    const int halfway=int(ngenes/2);
    double* combined=(double*)R_alloc(ngenes, sizeof(double));
    
    for (cell=0; cell<ncells; ++cell) {
        std::fill(combined, combined+ngenes, 0);

        for (index=0; index<SIZE; ++index) {
            const int& curcell=orptr[index+cell];
            oxptrs[curcell][cell]=1;

            for (gene=0; gene<ngenes; ++gene) { 
                combined[gene]+=eptrs[curcell][gene];
            }
        }

        // Computing the median ratio against the reference.
        for  (gene=0; gene<ngenes; ++gene) { 
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

