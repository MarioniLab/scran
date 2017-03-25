#include "scran.h"

template <typename T> 
void average_ranks(const T* ptr, const matrix_info& MAT, const int* subset, const int subset_len, double* incoming) {
    std::deque<std::pair<T, int> > collected(subset_len);
    size_t s;
    double accumulated_rank, sum_squares;
    size_t n_same_rank, s_same;
    const double mean_adj=double(subset_len-1)/2;

    for (size_t c=0; c<MAT.ncol; ++c) {
        for (s=0; s<subset_len; ++s) {
            if (isNA(ptr[s])) { 
                throw std::runtime_error("missing values not supported in quickCluster");
            }
            collected[s].first=ptr[subset[s]];
            collected[s].second=s;
        }
        std::sort(collected.begin(), collected.end());

        // Need a bit more effort to deal with tied ranks.
        accumulated_rank=0;
        n_same_rank=0;
        sum_squares=0;
        for (s=0; s<subset_len; ++s) {
            ++n_same_rank;
            accumulated_rank+=s;

            if (s==subset_len-1 || collected[s].first!=collected[s+1].first) {
                accumulated_rank /= n_same_rank;
                accumulated_rank -= mean_adj; // getting to a cosine distance.
                sum_squares += accumulated_rank * accumulated_rank * n_same_rank;
                
                s_same=s;
                while (n_same_rank) {
                    incoming[collected[s_same].second]=accumulated_rank;
                    --n_same_rank;
                    --s_same;
                }
                accumulated_rank=0;
            }
        }

        if (sum_squares==0) {
            throw std::runtime_error("rank variances of zero detected for a cell");
        }

        // Converting to cosine values.
        sum_squares = std::sqrt(sum_squares)*2;
        for (s=0; s<subset_len; ++s) {
            incoming[s]/=sum_squares;
        }

        ptr+=MAT.nrow;
        incoming+=subset_len;
    }
    return;
}

template <typename T>
SEXP cordist_internal(const T* ptr, const matrix_info& MAT, SEXP subset, SEXP return_ranks) {
    /// Checking the subset values.
    subset_values SS=check_subset_vector(subset, MAT.nrow);
    const int slen=SS.first;
    const int* sptr=SS.second;
    if (slen<2) {
        throw std::runtime_error("need at least 2 observations to compute correlations");
    }

    // Determining if ranks should be returned.
    if (!isLogical(return_ranks) || LENGTH(return_ranks)!=1) {
        throw std::runtime_error("return_ranks should be a logical scalar");
    }
    const bool rr=asLogical(return_ranks);

    // Computing the ranks.
    SEXP output1=PROTECT(allocMatrix(REALSXP, slen, MAT.ncol));
    try {
        average_ranks(ptr, MAT, sptr, slen, REAL(output1));
    } catch (std::exception& e) {
        UNPROTECT(1);
        throw;
    }
    if (rr) {
        UNPROTECT(1);
        return output1;
    }
   
    // Computing the distances between ranks (i.e., the correlation-based distance). 
    SEXP output2=PROTECT(allocMatrix(REALSXP, MAT.ncol, MAT.ncol));
    try {
        double ** optrs1=(double**)R_alloc(MAT.ncol, sizeof(double*));
        if (MAT.ncol) {
            optrs1[0]=REAL(output1);
            for (size_t cell=1; cell<MAT.ncol; ++cell) {
               optrs1[cell]=optrs1[cell-1]+slen;
           }   
        }

        double ** optrs2=(double**)R_alloc(MAT.ncol, sizeof(double*));
        if (MAT.ncol) {
            optrs2[0]=REAL(output2);
            for (size_t cell=1; cell<MAT.ncol; ++cell) {
               optrs2[cell]=optrs2[cell-1]+MAT.ncol;
           }   
        }

        // Running through every pair of cells and computing the correlation-based distance.
        const double* c1ranks, *c2ranks;
        size_t c1, c2;
        double curdiff;
        for (c1=0; c1<MAT.ncol; ++c1) {
            c1ranks=optrs1[c1];
            for (c2=0; c2<c1; ++c2) {
                c2ranks=optrs1[c2];
               
                // Converting to a distance metric based on the Spearman correlation.
                double& curdist=(optrs2[c1][c2]=0); 
                for (size_t s=0; s<slen; ++s) {
                    curdiff=c1ranks[s]-c2ranks[s];
                    curdist+=curdiff*curdiff;
                }
                curdist=std::sqrt(curdist);

                // Filling up the other guy.
                optrs2[c2][c1]=curdist;
            }
            optrs2[c1][c1]=0;
        }
    } catch (std::exception& e) {
        UNPROTECT(2);
        throw;
    }

    UNPROTECT(2);
    return output2;
}

SEXP compute_cordist(SEXP exprs, SEXP subset, SEXP return_ranks) try {
    matrix_info MAT=check_matrix(exprs);
    if (MAT.is_integer){ 
        return cordist_internal<int>(MAT.iptr, MAT, subset, return_ranks);
    } else {
        return cordist_internal<double>(MAT.dptr, MAT, subset, return_ranks);
    }
} catch (std::exception& e) {
    return mkString(e.what());
}
