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

/*** This function computes error-tolerant ranks for a subset of genes in a subset of cells. ***/
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

/*** This function computes residuals in a nice and quick manner, prior to ranking. ***/
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
        int cell, tmp;
        double working;
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


