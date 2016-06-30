#include "scran.h"

double rho_mult (double Ncells) {
    return 6/(Ncells*(Ncells*Ncells-1));
}

/*** Null distribution estimation without a design matrix. ***/
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
        double tmp, tmpdiff;
        const double mult=rho_mult(Ncells);

        for (int it=0; it<Niters; ++it) {
            Rx_shuffle(rankings, rankings + Ncells);
            tmp=0;
            for (cell=0; cell<Ncells; ++cell) {
                tmpdiff=rankings[cell]-cell;
                tmp+=tmpdiff*tmpdiff;
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

/*** Null distribution estimation with a design matrix. ***/
SEXP get_null_rho_design(SEXP qr, SEXP qraux, SEXP iters) try {
    matrix_info QR=check_matrix(qr);
    if (QR.is_integer) {
        throw std::runtime_error("Q matrix must be double-precision");
    }
    if (!isReal(qraux) || size_t(LENGTH(qraux))!=QR.ncol) {
        throw std::runtime_error("QR auxiliary vector should be double-precision and of length 'ncol(Q)'");
    }
    const double* qrxptr=REAL(qraux);
    if (!isInteger(iters) || LENGTH(iters)!=1) {
        throw std::runtime_error("number of iterations should be an integer scalar");
    }
    const int Niters=asInteger(iters);
    if (Niters <= 0) { throw std::runtime_error("number of iterations should be positive"); }
    
    // Setting up to multiply by the Q matrix.
    run_dormqr multQ(QR.nrow, QR.ncol, QR.dptr, qrxptr, 'N');
    const int& Nobs=QR.nrow;
    const int& Ncoef=QR.ncol;
    double* effects=multQ.rhs;

    SEXP output=PROTECT(allocVector(REALSXP, Niters));
    try {
        double* optr=REAL(output);
        Rx_random_seed myseed;
    
        std::deque<std::pair<double, int> > collected1(Nobs), collected2(Nobs);
        std::deque<int> rank1(Nobs), rank2(Nobs);
        const double mult=rho_mult(Nobs);

        // Simulating residuals, using the Q-matrix to do it.
        // We set the main effects to zero (hence, starting from "Ncoefs") and simulate normals for the residual effects.
        // We then use this to reconstruct the residuals themselves, and then compute correlations between them.
        int mode, row, col;
        double tmpdiff;

        for (int it=0; it<Niters; ++it) {
            for (mode=0; mode<2; ++mode) {
                for (col=0; col<Ncoef; ++col) {
                    effects[col]=0;
                }
                for (col=Ncoef; col<Nobs; ++col) {
                    effects[col]=norm_rand();
                }

                // Computing the residuals.
                multQ.run();

                // Sorting.
                std::deque<std::pair<double, int> >& current=(mode ? collected1 : collected2);
                for (row=0; row<Nobs; ++row) {
                    current[row].first=effects[row];
                    current[row].second=row;
                }
                std::sort(current.begin(), current.end());
                std::deque<int>& rank=(mode ? rank1 : rank2);
                for (row=0; row<Nobs; ++row) {
                    rank[current[row].second]=row;
                }
            }

            // Computing the squared difference in the ranks.
            double& rho=(optr[it]=0);
            for (row=0; row<Nobs; ++row) {
                tmpdiff=rank1[row]-rank2[row];
                rho+=tmpdiff*tmpdiff;
            }
            rho*=mult;
            rho=1-rho;
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

/*** Estimating correlations (without expanding into a matrix to do so via 'cor'). ***/
SEXP compute_rho(SEXP g1, SEXP g2, SEXP rankings) try {
    const matrix_info rmat=check_matrix(rankings);
    if (!rmat.is_integer) {
        throw std::runtime_error("rankings must be integer");
    }
    const int* rptr=rmat.iptr;
    const int& Ncells=rmat.nrow;
    if (Ncells <= 1) { throw std::runtime_error("number of cells should be greater than 2"); }
    const int& Ngenes=rmat.ncol;

    if (!isInteger(g1) || !isInteger(g2)) { 
        throw std::runtime_error("gene indices must be integer vectors");
    }
    const int Npairs=LENGTH(g1);
    if (Npairs!=LENGTH(g2)) { 
        throw std::runtime_error("gene index vectors must be of the same length"); 
    }
    const int *g1ptr=INTEGER(g1), *g2ptr=INTEGER(g2);
    
    SEXP output=PROTECT(allocVector(REALSXP, Npairs));
    try {
        double* orptr=REAL(output);
        
        const int* r1ptr, * r2ptr;
        int cell, tmp, working;
        const double mult=rho_mult(Ncells); 

        for (int p=0; p<Npairs; ++p) {
            const int& g1x=g1ptr[p];
            const int& g2x=g2ptr[p];
            if (g1x < 1 || g1x > Ngenes) {
                throw std::runtime_error("first gene index is out of range");
            }
            if (g2x < 1 || g2x > Ngenes) {
                throw std::runtime_error("second gene index is out of range");
            }
            r1ptr=rptr+(g1x-1)*Ncells;
            r2ptr=rptr+(g2x-1)*Ncells;

            // Computing the correlation.
            working=0;
            for (cell=0; cell<Ncells; ++cell) {
                tmp=r1ptr[cell] - r2ptr[cell];
                working+=tmp*tmp;
            }
            orptr[p]=1 - working*mult;
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


