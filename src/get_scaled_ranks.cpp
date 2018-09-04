#include "scran.h"

template <typename T, class V, class M> 
SEXP average_ranks_internal(const M mat, SEXP intype, SEXP subset, SEXP transpose, SEXP as_sparse) { 
    /// Checking the subset values.
    const size_t ncells=mat->get_ncol();
    const size_t ngenes=mat->get_nrow();
    auto SS=check_subset_vector(subset, ngenes);
    const size_t slen=SS.size();

    // Checking if we should transpose or not.
    const bool do_transpose=check_logical_scalar(transpose, "transpose specification");
    int out_nr=(do_transpose ? ncells : slen);
    int out_nc=(do_transpose ? slen : ncells);

    // Creating the output matrix.
    const bool sparsify=check_logical_scalar(as_sparse, "sparse specification");
    auto omat=beachmat::create_numeric_output(out_nr, out_nc, 
        sparsify ? beachmat::SPARSE_PARAM : beachmat::SIMPLE_PARAM);

    // Various other bits and pieces.
    std::vector<std::pair<T, int> > collected(slen);
    const double mean_adj=double(slen-1)/2;
    V incoming(ngenes);
    Rcpp::NumericVector outgoing(slen);

    for (size_t c=0; c<ncells; ++c) {
        mat->get_col(c, incoming.begin());

        // Sorting all subsetted values.
        auto sIt=SS.begin();
        for (size_t s=0; s<slen; ++s, ++sIt) {
            const T& curval=incoming[*sIt];
            if (isNA(curval)) { 
                throw std::runtime_error("missing values not supported in quickCluster");
            }
            collected[s].first=curval;
            collected[s].second=s;
        }
        std::sort(collected.begin(), collected.end());

        // Need a bit more effort to deal with tied ranks.
        double accumulated_rank=0, sum_squares=0;
        size_t n_same_rank=0;
        size_t most_common_rank=0, common_rank_begin=0;
        double common_rank_value=0;

        for (size_t s=0; s<slen; ++s) {
            ++n_same_rank;
            accumulated_rank+=s;

            if (s==slen-1 || collected[s].first!=collected[s+1].first) {
                accumulated_rank /= n_same_rank;
                accumulated_rank -= mean_adj; // getting to a cosine distance.
                sum_squares += accumulated_rank * accumulated_rank * n_same_rank;
                
                if (sparsify && n_same_rank > most_common_rank) {
                    most_common_rank = n_same_rank;
                    common_rank_begin = s;
                    common_rank_value = accumulated_rank;
                }

                size_t s_same=s;
                while (n_same_rank) {
                    outgoing[collected[s_same].second]=accumulated_rank;
                    --n_same_rank;
                    --s_same;
                }
                accumulated_rank=0;
            }
        }

        // Sparsifying by subtracting off the most common rank.
        if (sparsify) {
            for (auto& o : outgoing) {
                o -= common_rank_value;
            }
            for (size_t s=common_rank_begin; s<most_common_rank; ++s) {
                outgoing[collected[s].second]=0; // force to zero, just in case.
            }
        }

        // Converting to cosine values.
        if (sum_squares==0) {
            throw std::runtime_error("rank variances of zero detected for a cell");
        }
        sum_squares = std::sqrt(sum_squares)*2;
        for (auto oIt=outgoing.begin(); oIt!=outgoing.end(); ++oIt) {
            (*oIt)/=sum_squares;
        }

        if (do_transpose) { 
            omat->set_row(c, outgoing.begin());
        } else {
            omat->set_col(c, outgoing.begin());
        }
    }

    return omat->yield();
}

SEXP get_scaled_ranks(SEXP exprs, SEXP subset, SEXP transpose, SEXP as_sparse) {
    BEGIN_RCPP
    int rtype=beachmat::find_sexp_type(exprs);
    if (rtype==INTSXP) { 
        auto mat=beachmat::create_integer_matrix(exprs);
        return average_ranks_internal<int, Rcpp::IntegerVector>(mat.get(), exprs, subset, transpose, as_sparse);
    } else {
        auto mat=beachmat::create_numeric_matrix(exprs);
        return average_ranks_internal<double, Rcpp::NumericVector>(mat.get(), exprs, subset, transpose, as_sparse);
    }
    END_RCPP
}
