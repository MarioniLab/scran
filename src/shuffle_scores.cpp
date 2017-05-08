#include "scran.h"
#include <cstring>

template <typename T>
double get_proportion (const T* expr, const int& npairs, const int& minpairs, const int * m1, const int * m2) {
    int was_first=0, was_total=0;
    for (int p=0; p<npairs; ++p) {
        const T& first=expr[m1[p]];
        const T& second=expr[m2[p]];
        if (first > second) { ++was_first; }
        if (first != second) { ++was_total; }      
    }
    if (was_total < minpairs) { return(NA_REAL); }
    return(double(was_first)/was_total);
}

template <typename T>
SEXP shuffle_scores_internal (SEXP mycells, const T* eptr, const matrix_info& emat, 
        SEXP marker1, SEXP marker2, SEXP used, SEXP iter, SEXP miniter, SEXP minpair) { 

    if (!isInteger(mycells)) { throw std::runtime_error("cell indices must be an integer vector"); }
    const int* cptr=INTEGER(mycells);
    const int nc=LENGTH(mycells);
    const int ng=int(emat.nrow);
    const int totalcells=int(emat.ncol);

    if (!isInteger(marker1) || !isInteger(marker2)) { throw std::runtime_error("vectors of markers must be integer"); }
    const int npairs = LENGTH(marker1);
    if (npairs!=LENGTH(marker2)) { throw std::runtime_error("vectors of markers must be of the same length"); }
    if (!isInteger(used)) { throw std::runtime_error("vector of used gene indices must be integer"); }
    const int nused=LENGTH(used);
    const int* m1_ptr=INTEGER(marker1), * m2_ptr=INTEGER(marker2), * uptr=INTEGER(used);

    if (!isInteger(iter) || LENGTH(iter)!=1) { throw std::runtime_error("number of iterations must be an integer scalar"); }
    const int nit=asInteger(iter);
    if (!isInteger(miniter) || LENGTH(miniter)!=1) { throw std::runtime_error("minimum number of iterations must be an integer scalar"); }
    const int minit=asInteger(miniter);
    if (!isInteger(minpair) || LENGTH(minpair)!=1) { throw std::runtime_error("minimum number of pairs must be an integer scalar"); }
    const int minp=asInteger(minpair);

    // Checking marker sanity.    
    for (int marker=0; marker<npairs; ++marker) {
        const int& m1m=m1_ptr[marker];
        if (m1m >= nused || m1m < 0) { throw std::runtime_error("first marker indices are out of range"); }
        const int& m2m=m2_ptr[marker];
        if (m2m >= nused || m2m < 0) { throw std::runtime_error("second marker indices are out of range"); }
    }

    // Checking gene index sanity.
    for (int u=0; u<nused; ++u) {
        const int& usedex=uptr[u];
        if (usedex >= ng || usedex < 0) { throw std::runtime_error("used gene indices are out of range"); }
    }

    SEXP output=PROTECT(allocVector(REALSXP, nc));
    try {
        // Assorted temporary objects.
        double* optr=REAL(output);
        double curscore, newscore;
        int gene;
        int it, below, total;

        Rx_random_seed my_seed;
        const T* cell_exprs;
        T* current_exprs=(T*)R_alloc(nused, sizeof(T));

        for (int cell=0; cell<nc; ++cell) {
            const int& curcell=cptr[cell];
            if (curcell < 1 || curcell > totalcells) {
                throw std::runtime_error("cell indices are out of range");
            }

            // Storing the values to be shuffled in a separate array.
            cell_exprs=eptr + ng * (curcell - 1);
            for (gene=0; gene<nused; ++gene) { current_exprs[gene]=cell_exprs[uptr[gene]]; }
                
            curscore=get_proportion<T>(current_exprs, npairs, minp, m1_ptr, m2_ptr);
            if (ISNA(curscore)) { 
                optr[cell]=NA_REAL;
                continue;
            }

            // Iterations of shuffling to obtain a null distribution for the score.
            below=total=0;
            for (it=0; it < nit; ++it) {
                Rx_shuffle(current_exprs, current_exprs+nused);
                newscore=get_proportion<T>(current_exprs, npairs, minp, m1_ptr, m2_ptr);
                if (!ISNA(newscore)) { 
                    if (newscore < curscore) { ++below; }
                    ++total;
                }
            }
            
            optr[cell]=(total < minit ?  NA_REAL : double(below)/total);
        }
    } catch (std::exception& e) {
        UNPROTECT(1);
        throw;
    }

    UNPROTECT(1);
    return output;
}

SEXP shuffle_scores(SEXP mycells, SEXP exprs, SEXP marker1, SEXP marker2, SEXP indices, SEXP iter, SEXP miniter, SEXP minpair) try {
    const matrix_info emat=check_matrix(exprs);
    if (emat.is_integer) {
        return shuffle_scores_internal<int>(mycells, emat.iptr, emat, marker1, marker2, indices, iter, miniter, minpair);
    } else {
        return shuffle_scores_internal<double>(mycells, emat.dptr, emat, marker1, marker2, indices, iter, miniter, minpair);
    }
} catch (std::exception& e) {
    return mkString(e.what());
}

/* We could just assign ties random directions; then we'd only have to shuffle
 * once for all cells, and then we could use the same null distribution across
 * multiple cells, without worrying about whether or not one cell has more ties
 * than the other. The problem is that there's no protection from spuriously
 * high scores due to random breaking of ties; normally (for correlations),
 * we'd provide protection by controlling the type I error rate, but we're not
 * generating p-values here so it's harder to do.
 */

SEXP auto_shuffle(SEXP incoming, SEXP nits) {
    const int N=LENGTH(incoming);
    const int niters=asInteger(nits);
    SEXP output=PROTECT(allocMatrix(REALSXP, N, niters));
    {
        Rx_random_seed myseed;
        double* optr=REAL(output);
        const double* source=REAL(incoming);
        for (int i=0; i<niters; ++i) {
            std::memcpy(optr, source, N*sizeof(double));
            Rx_shuffle(optr, optr+N);
            source=optr;
            optr+=N;
        }
    }
    UNPROTECT(1);
    return(output);
}
