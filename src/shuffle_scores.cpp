#include "scran.h"
#include <cstring>

double get_proportion (const double* expr, const int* access, const int * permuted_index, const int& npairs, const int& minpairs, const int * m1, const int * m2) {
    // Assumes m1, m2 are 1-based, and permuted_index is already shifted by 1 to compensate.
    int was_first=0, was_total=0;
    for (int p=0; p<npairs; ++p) {
        const double& first=expr[permuted_index[access[m1[p]]]];
        const double& second=expr[permuted_index[access[m2[p]]]];
        if (first > second) { ++was_first; }
        if (first != second) { ++was_total; }      
    }
    if (was_total < minpairs) { return(NA_REAL); }
    return(double(was_first)/was_total);
}

SEXP shuffle_scores (SEXP mycells, SEXP ngenes, SEXP exprs, SEXP marker1, SEXP marker2, SEXP iter, SEXP miniter, SEXP minpair) try {
    if (!isInteger(mycells) || LENGTH(mycells)!=2) { throw std::runtime_error("cell indices must be an integer vector of length 2"); }
    const int start=INTEGER(mycells)[0] - 1;
    const int last=INTEGER(mycells)[1];
    const int nc=last-start;
    
    if (!isInteger(ngenes) || LENGTH(ngenes)!=1) { throw std::runtime_error("number of genes must be an integer scalar"); }
    const int ng=asInteger(ngenes);
    if (!isReal(exprs)) { throw std::runtime_error("matrix of expression values must be double-precision"); }
    if (nc*ng > LENGTH(exprs)) { throw std::runtime_error("size of expression matrix is not consistent with provided dimensions"); }
    if (!isInteger(marker1) || !isInteger(marker2)) { throw std::runtime_error("vectors of marker pair genes must be integer"); }
    const int npairs = LENGTH(marker1);
    if (npairs!=LENGTH(marker2)) { throw std::runtime_error("vectors of marker pairs must be of the same length"); }

    if (!isInteger(iter) || LENGTH(iter)!=1) { throw std::runtime_error("number of iterations must be an integer scalar"); }
    const int nit=asInteger(iter);
    if (!isInteger(miniter) || LENGTH(miniter)!=1) { throw std::runtime_error("minimum number of iterations must be an integer scalar"); }
    const int minit=asInteger(miniter);
    if (!isInteger(minpair) || LENGTH(minpair)!=1) { throw std::runtime_error("minimum number of pairs must be an integer scalar"); }
    const int minp=asInteger(minpair);

    const double** exp_ptrs=(const double**)R_alloc(nc, sizeof(const double*));
    if (nc) { exp_ptrs[0] = REAL(exprs) + start * ng; }
    int cell;
    for (cell=1; cell<nc; ++cell) { exp_ptrs[cell]=exp_ptrs[cell-1] + ng; }
    const int* m1_ptr=INTEGER(marker1), * m2_ptr=INTEGER(marker2);

    /* The idea is that the indices in each pair point to 'access', which in turn point to an entry of 'ref_permute/repermute'.
     * The permutation vector contains indices pointing to the rows of 'exprs' for each gene, for all and only rows used in any of the pairs.
     * Shuffling the values of the permutation vector is equivalent to permuting only the used rows amongst themselves.
     */
    int* access=(int*)R_alloc(ng, sizeof(int)) - 1;
    for (int g=1; g<=ng; ++g) { access[g]=-1; }
    int* ref_repermute=(int*)R_alloc(ng, sizeof(int));
    int stored=0;

    // Checking marker sanity.    
    for (int marker=0; marker<npairs; ++marker) {
        if (m1_ptr[marker] > ng || m1_ptr[marker] <= 0) { throw std::runtime_error("first marker indices are out of range"); }
        if (access[m1_ptr[marker]] < 0) { 
            access[m1_ptr[marker]]=stored;
            ref_repermute[stored]=m1_ptr[marker]-1;
            ++stored;
        }
        if (m2_ptr[marker] > ng || m2_ptr[marker] <= 0) { throw std::runtime_error("second marker indices are out of range"); }
        if (access[m2_ptr[marker]] < 0) { 
            access[m2_ptr[marker]]=stored;
            ref_repermute[stored]=m2_ptr[marker]-1;
            ++stored;
        }
    }

    int* repermute=(int*)R_alloc(stored, sizeof(int));
    SEXP output=PROTECT(allocVector(REALSXP, nc));
    try {
        double* optr=REAL(output);
        double curscore, newscore;
        int gene;
        int it, below, total;
        Rx_random_seed my_seed;

        for (int cell=0; cell<nc; ++cell) {
            for (gene=0; gene<stored; ++gene) { repermute[gene]=ref_repermute[gene]; }
            curscore=get_proportion(exp_ptrs[cell], access, repermute, npairs, minp, m1_ptr, m2_ptr);
            if (ISNA(curscore)) { 
                optr[cell]=NA_REAL;
                continue;
            }


            below=total=0;
            for (it=0; it < nit; ++it) {
                Rx_shuffle(repermute, repermute+stored);
                newscore=get_proportion(exp_ptrs[cell], access, repermute, npairs, minp, m1_ptr, m2_ptr);
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
