#include "scran.h"

template <typename T> 
double** average_ranks(const T* ptr, const matrix_info& MAT, const int* subset, const int subset_len) {
    /* Note that the vectors returned in this function correspond to the ordering in 'subset',
     * not any ordering of 0->MAT.nrow. In other words, it's as if MAT was already subsetted.
     */
    double** optrs=(double**)R_alloc(MAT.ncol, sizeof(double*));
    for (size_t c=0; c<MAT.ncol; ++c) {
        optrs[c]=(double*)R_alloc(subset_len, sizeof(double));
    }

    std::deque<std::pair<T, int> > collected(subset_len);
    for (size_t c=0; c<MAT.ncol; ++c) {
        for (size_t s=0; s<subset_len; ++s) {
            if (isNA(ptr[s])) { 
                throw std::runtime_error("missing values not supported in quickCluster");
            }
            collected[s].first=ptr[subset[s]];
            collected[s].second=s;
        }
        std::sort(collected.begin(), collected.end());
        
        // Need a bit more effort to deal with tied ranks.
        double accumulated_rank=0;
        size_t n_same_rank=0, s_same;
        for (size_t s=0; s<subset_len; ++s) {
            ++n_same_rank;
            accumulated_rank+=s;
            if (s==subset_len-1 || collected[s].first!=collected[s+1].first) {
                accumulated_rank/=n_same_rank;
                ++accumulated_rank; // to get to 1-based ranks.

                s_same=s;
                while (n_same_rank) {
                    optrs[c][collected[s_same].second]=accumulated_rank;
                    --n_same_rank;
                    --s_same;
                }
                accumulated_rank=0;
            }
        }

        ptr+=MAT.nrow;
    }
    return optrs;
}

template <typename T>
SEXP cordist_internal(const T* ptr, const matrix_info& MAT, SEXP subset) {
    if (!isInteger(subset)) { 
        throw std::runtime_error("subset vector must be integer");
    } 
    const int slen=LENGTH(subset);
    if (slen<2) {
        throw std::runtime_error("need at least 2 observations to compute correlations");
    }
    const int* sptr=INTEGER(subset);
    for (size_t s=0; s<slen; ++s) {
        if (sptr[s]<0 || sptr[s]>=MAT.nrow) {
            throw std::runtime_error("subset indices are out of range"); 
        }
    }

    SEXP output=PROTECT(allocMatrix(REALSXP, MAT.ncol, MAT.ncol));
    try {
        double ** optrs=(double**)R_alloc(MAT.ncol, sizeof(double*));
        if (MAT.ncol) {
            optrs[0]=REAL(output);
            for (size_t cell=1; cell<MAT.ncol; ++cell) {
               optrs[cell]=optrs[cell-1]+MAT.ncol;
           }   
        }

        /* Ranking the expression of all genes. Note that the ranks are
         * in terms of the subset, not in terms of 0->MAT.nrow.
         */
        double** ranked=average_ranks(ptr, MAT, sptr, slen);

        // Getting the standard deviation of the ranks for each cell.
        std::deque<double> stdev(MAT.ncol), means(MAT.ncol);
        double tmpdiff;
        for (size_t c=0; c<MAT.ncol; ++c) {
            double& curmean=means[c];
            for (size_t s=0; s<slen; ++s) {
                curmean+=ranked[c][s];
            }
            curmean/=slen;

            double& curdev=stdev[c];
            const double* currank=ranked[c];
            for (size_t s=0; s<slen; ++s) {
                tmpdiff=currank[s]-curmean;
                curdev+=tmpdiff*tmpdiff;
            }
            curdev/=slen-1;
            if (curdev==0) {
                throw std::runtime_error("rank variances of zero detected for a cell");
            }
            curdev=std::sqrt(curdev);
        }

        // Running through every pair of cells and computing the correlation-based distance.
        const double* c1ranks, *c2ranks;
        for (size_t c1=0; c1<MAT.ncol; ++c1) {
            c1ranks=ranked[c1];
            const double& c1mean=means[c1];
            const double& c1dev=stdev[c1];

            for (size_t c2=0; c2<c1; ++c2) {
                c2ranks=ranked[c2];
                const double& c2mean=means[c2];
                const double& c2dev=stdev[c2];
                
                // Converting to a distance metric based on the Spearman correlation.
                double& curdist=(optrs[c1][c2]=0); 
                for (size_t s=0; s<slen; ++s) {
                    curdist+=(c1ranks[s]-c1mean)*(c2ranks[s]-c2mean);
                }
                curdist/=slen-1;
                curdist/=c1dev*c2dev;
                if (curdist >= 1) {
                    curdist=0; // avoid rooting negatives due to numerical imprecision.
                } else {
                    curdist=std::sqrt(0.5*(1-curdist));
                }

                // Filling up the other guy.
                optrs[c2][c1]=curdist;
            }
            optrs[c1][c1]=0;
        }
    } catch (std::exception& e) {
        UNPROTECT(1);
        throw;
    }

    UNPROTECT(1);
    return output;
}

SEXP compute_cordist(SEXP exprs, SEXP subset) try {
    matrix_info MAT=check_matrix(exprs);
    if (MAT.is_integer){ 
        return cordist_internal<int>(MAT.iptr, MAT, subset);
    } else {
        return cordist_internal<double>(MAT.dptr, MAT, subset);
    }
} catch (std::exception& e) {
    return mkString(e.what());
}
