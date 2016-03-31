#include "scran.h"

double rho_mult (double Ncells) {
    return 6/(Ncells*(Ncells*Ncells-1));
}

SEXP get_null_rho (SEXP cells, SEXP iters) try {
    if (!isInteger(cells) || LENGTH(cells)!=1)  {
        throw std::runtime_error("number of cells should be an integer scalar"); 
    }
    const int Ncells=asInteger(cells);
    if (Ncells <= 1) { throw std::runtime_error("number of cells should be greater than 2"); }
    if (!isInteger(iters) || LENGTH(iters)!=1)  {
        throw std::runtime_error("number of iterations should be an integer scalar"); 
    }
    const int Niters=asInteger(iters);
    if (Niters <= 0) { throw std::runtime_error("number of iterations should be positive"); }

    int* rankings=(int*)R_alloc(Ncells, sizeof(int));
    int cell;
    for (cell=0; cell<Ncells; ++cell) { rankings[cell]=cell; }    

    SEXP output=PROTECT(allocVector(REALSXP, Niters));
    try {
        double* optr=REAL(output);
        Rx_random_seed myseed;
        double tmp;
        const double mult=rho_mult(Ncells);

        for (int it=0; it<Niters; ++it) {
            Rx_shuffle(rankings, rankings + Ncells);
            tmp=0;
            for (cell=0; cell<Ncells; ++cell) {
                tmp+=(rankings[cell]-cell)*(rankings[cell]-cell);
            }
            tmp*=mult;
            optr[it]=1-tmp;            
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

SEXP compute_rho(SEXP g1, SEXP g2, SEXP cells, SEXP rankings, SEXP nulls) try {
    if (!isInteger(cells) || LENGTH(cells)!=1)  {
        throw std::runtime_error("number of cells should be an integer scalar"); 
    }
    const int Ncells=asInteger(cells);
    if (Ncells <= 1) { throw std::runtime_error("number of cells should be greater than 2"); }

    if (!isInteger(g1) || !isInteger(g2)) { 
        throw std::runtime_error("gene indices must be integer vectors");
    }
    const int Npairs=LENGTH(g1);
    if (Npairs!=LENGTH(g2)) { 
        throw std::runtime_error("gene index vectors must be of the same length"); 
    }
    const int *g1ptr=INTEGER(g1), *g2ptr=INTEGER(g2);
    
    if (!isInteger(rankings)) { 
        throw std::runtime_error("ranking matrix must be integer");
    }
    const int Ngenes=LENGTH(rankings)/Ncells;
    if (Ngenes*Ncells!=LENGTH(rankings)) { 
        throw std::runtime_error("dimensions of ranking matrix are not consistent with number of cells"); 
    }
    const int* rptr=INTEGER(rankings);

    if (!isReal(nulls)) {
        throw std::runtime_error("null distribution must be double-precision");
    }
    const int Nnulls=LENGTH(nulls);
    double* nptr=REAL(nulls), *after=nptr+Nnulls;

    SEXP output=PROTECT(allocVector(VECSXP, 2));
    try {
        SET_VECTOR_ELT(output, 0, allocVector(REALSXP, Npairs));
        double* orptr=REAL(VECTOR_ELT(output, 0));
        SET_VECTOR_ELT(output, 1, allocVector(REALSXP, Npairs));
        double* opptr=REAL(VECTOR_ELT(output, 1));
        
        int off1, off2, cell;
        double tmp, tmp2;
        const double mult=rho_mult(Ncells); 
        double* finder, *ref_finder;
        int left, right;

        for (int p=0; p<Npairs; ++p) {
            off1=(g1ptr[p]-1)*Ncells;
            off2=(g2ptr[p]-1)*Ncells;
            if (g1ptr[p] < 1 || g1ptr[p] > Ngenes) {
                throw std::runtime_error("first gene index is out of range");
            }
            if (g2ptr[p] < 1 || g2ptr[p] > Ngenes) {
                throw std::runtime_error("second gene index is out of range");
            }

            // Computing the correlation.
            tmp=0;
            for (cell=0; cell<Ncells; ++cell) {
                tmp2=rptr[off1+cell] - rptr[off2+cell];
                tmp+=tmp2*tmp2;
            }
            tmp*=mult;
            orptr[p]=1-tmp;
            
            // Computing the p-value.
            ref_finder=std::lower_bound(nptr, after, orptr[p]);
            finder=ref_finder;
            right=after-finder; // equal to or greater than the current rho.
            while (finder!=nptr) {
                --finder;
                if (orptr[p] - *finder > 0.00000001) { break; }
                ++right;
            }

            finder=ref_finder;
            left=finder-nptr; // don't add 1, as we don't know it's equal to.
            while (finder!=after && *finder - orptr[p] < 0.00000001) { // checking for equal to's.
                ++finder;
                ++left;
            }
            
            opptr[p]=double(std::min(right, left)+1)*2/(Nnulls+1);
            if (opptr[p] > 1) { opptr[p]=1; }
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


