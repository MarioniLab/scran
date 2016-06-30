#include "scran.h"

matrix_info::matrix_info (int nr, int nc, bool ii) : nrow(nr), ncol(nc), is_integer(ii), iptr(NULL), dptr(NULL) {}

matrix_info check_matrix(SEXP matrix) {
    int type;
    if (isReal(matrix)) {
        type=0;
    } else if (isInteger(matrix)) {
        type=1;
    } else if (isLogical(matrix)) { 
        type=2;
    } else {
        throw std::runtime_error("matrix must be integer or double-precision");
    }

    SEXP dims=getAttrib(matrix, R_DimSymbol);
    if (!isInteger(dims) || LENGTH(dims)!=2) { 
        throw std::runtime_error("dimensions of the matrix should be an integer vector of length 2");
    }
    int nrow=INTEGER(dims)[0], ncol=INTEGER(dims)[1];
    if (LENGTH(matrix)!=nrow*ncol) {
        throw std::runtime_error("recorded dimensions of the matrix are not consistent with its length"); 
    }

    matrix_info output(nrow, ncol, type>0);
    switch (type) {
        case 0:
            output.dptr=REAL(matrix);
            break;
        case 1:
            output.iptr=INTEGER(matrix);
            break;
        case 2:
            output.iptr=LOGICAL(matrix);
            break;
    }

    return output;
}

bool isNA(int x) {
    return x==NA_INTEGER;
}

bool isNA(double x) {
    return ISNA(x);
}

subset_values check_subset_vector(SEXP subset, int maxdim) {
    if (!isInteger(subset)) { 
        throw std::runtime_error("subset vector must be an integer vector");
    }
    const int slen=LENGTH(subset);
    const int* sptr=INTEGER(subset);
    for (int s=0; s<slen; ++s) {
        if (sptr[s] < 0 || sptr[s] >= maxdim) {
            throw std::runtime_error("subset indices out of range");
        }
    }
    return std::make_pair(slen, sptr);
}

